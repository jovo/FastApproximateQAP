/***************************************************************************
 *   Copyright (C) 2008 by Mikhail Zaslavskiy   *
 *   mikhail.zaslavskiy@ensmp.fr   *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/
#include "algorithm_sch.h"

match_result algorithm_sch::match(graph& g, graph& h,gsl_matrix* gm_P_i,gsl_matrix* gm_ldh,double dalpha_ldh)
{
	bool bblast_match_end=(get_param_i("blast_match_proj")==1);
	bool bgreedy=(get_param_i("hungarian_greedy")==1);
	match_result mres=algorithm_qcv::match(g,h,gm_P_i,gm_ldh,dalpha_ldh);
	if (bverbose)
		*gout<<"SCH algorithm"<<std::endl;	
	//some duplicate variables
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue); 
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	//memory allocation
	gsl_eigen_symmv_workspace * gesw= gsl_eigen_symmv_alloc (N);
	gsl_vector* eval_g=gsl_vector_alloc(N);
	gsl_vector* eval_h=gsl_vector_alloc(N);
	gsl_matrix* evec_g=gsl_matrix_alloc(N,N);
	gsl_matrix* evec_h=gsl_matrix_alloc(N,N);
	if (bverbose) *gout<<"Memory allocation finished"<<std::endl;
	//eigenvalues and eigenvectors for both matrices
	gsl_eigen_symmv (gm_Ag_d, eval_g,evec_g,gesw);
	if (bverbose) *gout<<"Ag eigen vectors"<<std::endl;
	gsl_eigen_symmv (gm_Ah_d, eval_h,evec_h,gesw);
	
	
	if (bverbose) *gout<<"Ah eigen vectors"<<std::endl;
	gsl_eigen_symmv_sort (eval_g, evec_g, GSL_EIGEN_SORT_VAL_DESC);
	gsl_eigen_symmv_sort (eval_h, evec_h, GSL_EIGEN_SORT_VAL_DESC);
	gsl_matrix_abs(evec_g);
	gsl_matrix_abs(evec_h);

	gsl_matrix *gm_lambda=gsl_matrix_alloc(N,N);
	gsl_matrix_set_all(gm_lambda,0);
	for (int i=0;i<eval_g->size;i++)
	for (int j=0;j<eval_h->size;j++)
		{
			gsl_matrix_set(gm_lambda,i,j,abs(eval_g->data[i]-eval_h->data[j]));
		};
	/*for (int i=0;i<eval_g->size;i++)
		{
			gsl_vector_view gvv_c=gsl_matrix_column(evec_g,i);
			gsl_vector_scale(&gvv_c.vector,pow(eval_g->data[i],0.5));
		};
	for (int j=0;j<eval_h->size;j++)
		{
			gsl_vector_view gvv_c=gsl_matrix_column(evec_h,j);
			gsl_vector_scale(&gvv_c.vector,pow(eval_h->data[j],0.5));
		};*/
	//scale matrix construction
	gsl_matrix* C=gsl_matrix_alloc(N,N);
	gsl_matrix_transpose(evec_h);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_lambda,evec_h,0,C);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,evec_g,C,0,gm_lambda);
	gsl_matrix_memcpy(C,gm_lambda);
	//produce weighted eigen vectors
	/*for (int i=0;i<eval_g->size;i++)
	for (int j=0;j<eval_h->size;j++)
		{
			C->data[i+N*j]=0;
			for (int k=0;k<N;k++)
			for (int m=0;m<N;m++)
				C->data[i+N*j]+=pow(gm_Ag_d->data[i+k*N]-gm_Ah_d->data[j+m*N],2);
		};*/
	
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
		
//	gsl_matrix_transpose(C);
	//gsl_matrix_scale(C,-1);
	gsl_matrix_printout(C,"C",pdebug.strvalue);
	//gsl_matrix_set_all(C,1);

	//multiplication of QCV minimum on C
	gsl_matrix_memcpy(gm_lambda,mres.gm_P_exact);
	gsl_matrix_add_constant(mres.gm_P_exact,-0.5);
	gsl_matrix_mul_elements(mres.gm_P_exact,C);
	gsl_matrix_transpose(mres.gm_P_exact);
	//permuation projection
	gsl_matrix_scale(mres.gm_P_exact,-10000);
	gsl_matrix_hungarian(mres.gm_P_exact,mres.gm_P,NULL,NULL,false,(bblast_match_end?gm_ldh:NULL),bgreedy);
	gsl_matrix_memcpy(mres.gm_P_exact,gm_lambda);
	//memory release
	gsl_matrix_free(C);
	gsl_matrix_free(gm_lambda);
	gsl_matrix_free(evec_g);
	gsl_matrix_free(evec_h);
	gsl_vector_free(eval_g);
	gsl_vector_free(eval_h);
	//initial score
	mres.vd_trace.clear();
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,mres.gm_P,cscore_matrix));
	//other output parameters
	mres.dres=mres.vd_trace.at(1);
	mres.inum_iteration=2;
	return mres;
}

