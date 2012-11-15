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
#include "algorithm_umeyama.h"

match_result algorithm_umeyama::match(graph& g, graph& h,gsl_matrix* gm_P_i,gsl_matrix* gm_ldh,double dalpha_ldh)
{
	if (bverbose)
		*gout<<"Umeyama algorithm"<<std::endl;	
	bool bblast_match_end=(get_param_i("blast_match_proj")==1);
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
	
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	
	if (bverbose) *gout<<"Ah eigen vectors"<<std::endl;
	gsl_eigen_symmv_sort (eval_g, evec_g, GSL_EIGEN_SORT_VAL_DESC);
	gsl_eigen_symmv_sort (eval_h, evec_h, GSL_EIGEN_SORT_VAL_DESC);
	
	if (pdebug.ivalue){ gsl_matrix_printout(eval_g,"eval_g",pdebug.strvalue); 
			    gsl_matrix_printout(eval_h,"eval_h",pdebug.strvalue); 
	 		    gsl_matrix_printout(evec_g,"evec_g",pdebug.strvalue);
			    gsl_matrix_printout(evec_h,"evec_h",pdebug.strvalue);};
	gsl_matrix_abs(evec_g);
	gsl_matrix_abs(evec_h);
	if (pdebug.ivalue) gsl_matrix_printout(evec_g,"abs(evec_g)",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(evec_h,"abs(evec_h)",pdebug.strvalue);  
	//loss matrix construction
	gsl_matrix* C=gsl_matrix_alloc(N,N);
	if (bverbose) *gout<<"Loss function matrix allocation"<<std::endl;
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,evec_g,evec_h,0,C);
	if (pdebug.ivalue) gsl_matrix_printout(C,"C=abs(evec_g)*abs(evec_h')",pdebug.strvalue); 
	//label cost matrix
	update_C_hungarian(C,-1);
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
	if (bverbose) *gout<<"Hungarian solved"<<std::endl;
	match_result mres;
	mres.gm_P=gm_P;
	

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

