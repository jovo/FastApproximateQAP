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
#include "algorithm_rank.h"

match_result algorithm_rank::match(graph& g, graph& h,gsl_matrix* gm_P_i, gsl_matrix* _gm_ldh,double dalpha_ldh)
{
	if (bverbose)
		*gout<<"Rank matching"<<std::endl;
        bool bblast_match_end=(get_param_i("blast_match_proj")==1);
	//some duplicate variables
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue); 
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	double d1=gsl_matrix_sum(gm_Ag_d);
	double d2=gsl_matrix_sum(gm_Ah_d);
	gsl_matrix* C=gsl_matrix_alloc(N,N);
	gsl_matrix* C_old=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_temp=gsl_matrix_alloc(N,N);
	gsl_matrix_set_all(C,1.0/(N));
	int max_num_it=1000;
	int num_it=0;
	
	gsl_vector_view gvv_C_old=gsl_vector_view_array(C_old->data,N*N);
	gsl_vector_view gvv_R=gsl_vector_view_array(C->data,N*N);
	
	
	gsl_vector* gv_deg_g=gsl_vector_alloc(N);
	gsl_vector* gv_deg_h=gsl_vector_alloc(N);
	gsl_matrix_sum(gm_Ag_d,1, gv_deg_g);
	gsl_matrix_sum(gm_Ah_d,1, gv_deg_h);
	//adjacency matrix normalization
	for (int i=0;i<N;i++)
		for (int j=0;j<N;j++)
		{
		 if (gv_deg_g->data[j]>0) gm_Ag_d->data[i*N+j]/=gv_deg_g->data[j];
		 if (gv_deg_h->data[j]>0) gm_Ah_d->data[i*N+j]/=gv_deg_h->data[j];
		};
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"gm_Ag_d",pdebug.strvalue); 
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"gm_Ah_d",pdebug.strvalue); 
	gsl_vector_free(gv_deg_g);
	gsl_vector_free(gv_deg_h);
	double ddiff;bool bcontinue=true;double dm_old;
	if (dalpha_ldh>0) { 
		gsl_matrix_scale(gm_ldh,dalpha_ldh);
		gsl_matrix_memcpy(C,gm_ldh);
			  }
		else
		      gsl_matrix_set_all(C,1.0/N);
	
	if (pdebug.ivalue) gsl_matrix_printout(C,"C",pdebug.strvalue); 
	bool biter_algo=true;
	if (biter_algo)
	{
		while (bcontinue)
		{	
			gsl_matrix_memcpy(C_old,C);
			gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, gm_Ag_d,C,0,gm_temp);
			if (pdebug.ivalue) gsl_matrix_printout(gm_temp,"gm_temp",pdebug.strvalue); 
			//gsl_blas_dgemm(CblasNoTrans,CblasTrans, 1-dalpha_ldh, gm_temp,gm_Ah_d,0,C);
			gsl_matrix_transpose(gm_temp);
			gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1-dalpha_ldh,gm_Ah_d, gm_temp,0,C);
			gsl_matrix_transpose(C);
			if (pdebug.ivalue) gsl_matrix_printout(C,"C",pdebug.strvalue); 
			if (dalpha_ldh>0) gsl_matrix_add(C,gm_ldh);
			//normalization
			dm_old=gsl_blas_dnrm2(&gvv_R.vector);
			gsl_matrix_scale(C,1.0/dm_old);
			gsl_matrix_sub(C_old,C);
			if (pdebug.ivalue) gsl_matrix_printout(C_old,"C_old",pdebug.strvalue); 
			ddiff=gsl_blas_dnrm2(&gvv_C_old.vector);
			bcontinue=!(ddiff<1e-2);
			bcontinue=(bcontinue && (num_it++<max_num_it));
			if (bverbose)
				*gout<<"Dif="<<ddiff<<std::endl;
		};
	}
	else //eigen decomposition approach
	{
		//memory allocation
		gsl_eigen_symmv_workspace * gesw= gsl_eigen_symmv_alloc (N);
		gsl_vector* eval_g=gsl_vector_alloc(N);
		gsl_vector* eval_h=gsl_vector_alloc(N);
		gsl_matrix* evec_g=gsl_matrix_alloc(N,N);
		gsl_matrix* evec_h=gsl_matrix_alloc(N,N);
		if (bverbose) *gout<<"Memory allocation "<<std::endl;
		//eigenvalues and eigenvectors for both matrices
		gsl_eigen_symmv (gm_Ag_d, eval_g,evec_g,gesw);
		if (bverbose) *gout<<"Ag eigen vectors"<<std::endl;
		gsl_eigen_symmv (gm_Ah_d, eval_h,evec_h,gesw);
		
		if (bverbose) *gout<<"Ah eigen vectors"<<std::endl;
		gsl_eigen_symmv_sort (eval_g, evec_g, GSL_EIGEN_SORT_VAL_DESC);
		gsl_eigen_symmv_sort (eval_h, evec_h, GSL_EIGEN_SORT_VAL_DESC);
		//eigenvalues recalculation
		
		gsl_matrix_view gmv_Lg=gsl_matrix_view_vector(eval_g,N,1);
		gsl_matrix_view gmv_Lh=gsl_matrix_view_vector(eval_h,1,N);
		gsl_matrix_view gmv_L=gsl_matrix_view_array(C_old->data,N,N);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, &gmv_Lg.matrix,&gmv_Lh.matrix,0,C_old);
		gsl_matrix_scale(C_old,-1);
		gsl_matrix_add_constant(C_old,1);

		if (pdebug.ivalue) gsl_matrix_printout(&gmv_Lg.matrix,"gmv_Lg.matrix",pdebug.strvalue); 
		if (pdebug.ivalue) gsl_matrix_printout(&gmv_Lh.matrix,"gmv_Lh.matrix",pdebug.strvalue); 
		if (pdebug.ivalue) gsl_matrix_printout(C_old,"C_old",pdebug.strvalue); 
		gsl_matrix_memcpy(C,gm_ldh);
		gsl_matrix_transpose(C);
		//multiplcation chain
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, C,gm_Ah_d,0,gm_temp);
		gsl_matrix_transpose(gm_Ag_d);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, gm_Ag_d,gm_temp,0,C);
		gsl_matrix_transpose(gm_Ag_d);
		
		gsl_matrix_div_elements(C,C_old);
		
		//gsl_matrix_transpose(gm_Ah_d);
		//gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, C,gm_Ah_d,0,gm_temp);
		gsl_matrix_transpose(C);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, gm_Ah_d,C,0,gm_temp);
		gsl_matrix_transpose(gm_temp);
		gsl_blas_dgemm(CblasNoTrans,CblasNoTrans, 1, gm_Ag_d,gm_temp,0,C);
		
		gsl_matrix_free(gm_Ag_d);
		gsl_matrix_free(gm_Ah_d);		
		gsl_vector_free(eval_g);
		gsl_vector_free(eval_h);
		gsl_matrix_free(evec_g);
		gsl_matrix_free(evec_h);
		gsl_matrix_free(C);	
	};
	if (dalpha_ldh>0) gsl_matrix_scale(gm_ldh,1/dalpha_ldh);
	if (pdebug.ivalue) gsl_matrix_printout(C,"C=rank matrix",pdebug.strvalue); 
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_matrix_free(C_old);
	gsl_matrix_free(gm_temp);
	
	//scaling for hungarian
	double dscale_factor =gsl_matrix_max_abs(C);
	dscale_factor=(dscale_factor>EPSILON)?dscale_factor:EPSILON;
	dscale_factor=10000/dscale_factor;
	gsl_matrix_scale(C,-dscale_factor);
	gsl_matrix_transpose(C);
	if (pdebug.ivalue) gsl_matrix_printout(C,"scale(C)",pdebug.strvalue); 
	gsl_matrix* gm_P=gsl_matrix_alloc(N,N);
	gsl_matrix_hungarian(C,gm_P,NULL,NULL,false,(bblast_match_end?gm_ldh:NULL),false);
        if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);
	
	match_result mres;
	mres.gm_P=gm_P;
	gsl_matrix_free(C);

	//initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,gm_P,cscore_matrix));
	//other output parameters
	mres.dres=mres.vd_trace.at(1);
	mres.inum_iteration=2;
	//transpose matrix save
	mres.gm_P=gm_P;
	mres.gm_P_exact=NULL;
	return mres;
}

