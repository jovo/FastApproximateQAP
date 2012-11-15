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
#include "algorithm_qcv.h"

match_result algorithm_qcv::match(graph& g, graph& h,gsl_matrix* gm_P_i,gsl_matrix* _gm_ldh,double dalpha_ldh)
{
	if (bverbose)
		*gout<<"QCV matching"<<std::endl;
	bool bblast_match_end=(get_param_i("blast_match_proj")==1);	
	double dfw_xeps=get_param_d("algo_fw_xeps");
	double dfw_feps=get_param_d("algo_fw_feps");
	double dhung_max=get_param_d("hungarian_max");
	double bgreedy=(get_param_i("hungarian_greedy")==1);
	//best path solutions
	bool bbest_path_proj=(get_param_i("best_path_proj_sol")==1);
	bool bbest_path_blast_proj=(get_param_i("best_path_blast_proj_sol")==1);
	bool bbest_path_greedy=(get_param_i("best_path_greedy_sol")==1);
	bool bbest_path_blast_greedy=(get_param_i("best_path_blast_greedy_sol")==1);
	bool bbest_path=bbest_path_proj or bbest_path_blast_proj or bbest_path_greedy or bbest_path_blast_greedy;
	gsl_matrix* gm_P_bp_temp=NULL;
	gsl_matrix* gm_P_bp=NULL;
	double fbest_path=1e+300;
	
	if (bbest_path)
	{
		gm_P_bp_temp=gsl_matrix_alloc(N,N);
		gm_P_bp=gsl_matrix_alloc(N,N);
	};

	//some duplicate variables
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag_d,"Ag",pdebug.strvalue); 
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah_d,"Ah",pdebug.strvalue);
	//Frank-Wolfe algorithm
	bool bstop_algo=false;
	gsl_vector* gv_C=gsl_vector_alloc(N*N);
	gsl_vector* gv_temp=gsl_vector_alloc(N*N);
	gsl_matrix* C;
	gsl_matrix_view gmv_C;
	gsl_vector_view gvv_P,gvv_dP,gvv_P_prev;
	gsl_matrix* gm_P=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_P_prev=gsl_matrix_alloc(N,N);
	gvv_P_prev=gsl_vector_view_array(gm_P_prev->data,N*N);
	gsl_vector* gv_col_inc=gsl_vector_alloc(N);
	gsl_vector_set_zero(gv_col_inc);
	gsl_permutation* gp_sol=gsl_permutation_alloc(N);
	gsl_permutation_init(gp_sol);
	
	
	gsl_matrix* gm_dP=gsl_matrix_alloc(N,N);
	gvv_dP=gsl_vector_view_array(gm_dP->data,N*N);
	gsl_matrix* gm_dP_2=gsl_matrix_alloc(N,N);
	gsl_vector_view gvv_dP_2=gsl_vector_view_array(gm_dP_2->data,N*N);
	gsl_matrix_set_zero(gm_dP_2);
	
	if (gm_P_i==NULL)
		gsl_matrix_set_all(gm_P,1.0/N);
	else
		gsl_matrix_memcpy(gm_P,gm_P_i);

	gsl_matrix_memcpy(gm_P_prev,gm_P);
	//perm matrix transformation into vector
	gvv_P=gsl_vector_view_array(gm_P->data,N*N);
	
	//and in opposite direction for gradient
	gmv_C=gsl_matrix_view_vector(gv_C,N,N);
	C=&gmv_C.matrix;
	gsl_vector*gv_debug_trace=gsl_vector_alloc(3);//debug trace information
	gsl_vector_set_zero(gv_debug_trace);
	gsl_matrix * gm_temp=gsl_matrix_alloc(N,N);
	gsl_matrix * gm_ddP_prev=gsl_matrix_alloc(N,N);
	gsl_vector_view gvv_ldh;
	if (dalpha_ldh>0)
	{ 
		gvv_ldh=gsl_vector_view_array(gm_ldh->data,N*N);
		if (pdebug.ivalue) gsl_matrix_printout(gm_ldh,"gm_ldh",pdebug.strvalue);
	};
	double dt0,dt1,dt2,dt3;
	double ddP_norm,dP_norm, dfvalue,dfvalue_prev;
	dfvalue_prev=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
	if (bverbose) *gout<<"QCVL: main cycle"<<std::endl;
	while(!bstop_algo)
	{
	dt0=clock();
	//the default gsl representation is made by rows
	//gradient estimation: Agh*P
	if (pdebug.ivalue) gsl_matrix_printout(&gvv_P.vector,"gvv_P",pdebug.strvalue);
	qcv_gradient_opt(gm_Ag_d,gm_Ah_d,&gvv_P.vector,gv_C,gm_temp);
	dt1=clock();
	double r=(dt1-dt0)/CLOCKS_PER_SEC;
	if (pdebug.ivalue) gsl_matrix_printout(gv_C,"gv_C",pdebug.strvalue);
		
	//result save
	gsl_matrix_transpose(C);
	if (pdebug.ivalue) gsl_matrix_printout(C,"C=gradient",pdebug.strvalue); 
	gsl_matrix_scale(C,2);
	update_C_hungarian(C);
	double dscale_factor =gsl_matrix_max_abs(C);
	dscale_factor=(dscale_factor>EPSILON)?dscale_factor:EPSILON;
	dscale_factor=dhung_max/dscale_factor;
	gsl_matrix_scale(C,dscale_factor);
	//gsl_matrix_transpose(C);
	if (pdebug.ivalue) gsl_matrix_printout(C,"scale(C)",pdebug.strvalue); 
	
	//hungarian, before the true C matrix must be transposed
	gsl_matrix_transpose(C);
	dt1=clock();
	gsl_matrix_hungarian(C,gm_P,gv_col_inc,gp_sol,false);
	dt2=clock();
	gsl_matrix_transpose(C);
	gsl_matrix_scale(C,1/dscale_factor);
	update_C_hungarian(C,1,true);//return to the original value
	gsl_matrix_scale(C,0.5);

	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);	
	if (pdebug.ivalue) gsl_matrix_printout(gm_P_prev,"gm_P_prev",pdebug.strvalue);  
	
	//line search
	gsl_matrix_memcpy(gm_dP,gm_P);
	gsl_matrix_sub(gm_dP,gm_P_prev);
	
	if (pdebug.ivalue) gsl_matrix_printout(gm_dP,"gm_dP",pdebug.strvalue); 
	
	double a,b1,b2,bldh;
	gsl_matrix_transpose(gm_dP);
	gsl_matrix_transpose(gm_P_prev);
	gsl_matrix_transpose(C);
	gsl_blas_ddot(&gvv_dP.vector,gv_C,&b1);
	gsl_matrix_transpose(gm_dP);
	qcv_gradient_opt(gm_Ag_d,gm_Ah_d,&gvv_dP.vector,gv_temp,gm_temp);
	gsl_matrix_transpose(gm_dP);

	gsl_blas_ddot(&gvv_P_prev.vector,gv_temp,&b2);
	if (dalpha_ldh>0){
	gsl_matrix_transpose(gm_ldh);
	gsl_blas_ddot(&gvv_dP.vector,&gvv_ldh.vector,&bldh);
	gsl_matrix_transpose(gm_ldh);}
	else bldh=0;
	bldh=-dalpha_ldh*bldh;
	b1+=b2+bldh;
	gsl_blas_ddot(&gvv_dP.vector,gv_temp,&a);
	if (a>0)
	{
		double alpha = -b1/(2*a);
		gsl_matrix_transpose(gm_dP);
		gsl_matrix_transpose(gm_P_prev);
		gsl_matrix_transpose(C);
		if ((alpha<1) and (alpha>0))
			{	
				gsl_matrix_scale(gm_dP,(1-alpha));
				gsl_matrix_sub(gm_P,gm_dP);
			};
		if (!(alpha>0))
			gsl_matrix_memcpy(gm_P,gm_P_prev);
		alpha=(alpha>1)?1:alpha;alpha=(alpha<0)?0:alpha;
		if (bverbose) *gout<<"dalpha=("<<alpha<<"); ";
	}
	else
	{
	     dfvalue=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
	     if (dfvalue>dfvalue_prev) 
			{gsl_matrix_memcpy(gm_P,gm_P_prev);
			 dfvalue=dfvalue_prev;};
	};
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_step_finish",pdebug.strvalue); //sparsity 
	long lsparse=0;
	for (long j=0;j<N*N;j++)
		lsparse+=(gm_P->data[j]<1e-7);
	if (bverbose) *gout<<"#zeros="<<lsparse<<", #nonzeros="<<N*N-lsparse<<std::endl;
			
	//stop criterion
	dP_norm=gsl_matrix_norm(gm_P_prev,1);
	gsl_matrix_sub(gm_P_prev,gm_P);
	dfvalue=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
	ddP_norm=gsl_matrix_norm(gm_P_prev,1);
	bstop_algo=((ddP_norm<dfw_xeps*N) and ((abs(dfvalue-dfvalue_prev)<dfw_feps*abs(dfvalue_prev)) or (ddP_norm==0)) or (abs(dfvalue-dfvalue_prev)<1e-20));
	dt3=clock();
	if (bverbose) *gout<<"x="<<dP_norm<<", f="<<dfvalue<<", dx="<<ddP_norm<<", df="<<(dfvalue_prev-dfvalue)<<", "<<"grad="<<(dfvalue_prev-dfvalue)/ddP_norm<<". Timing="<<(dt1-dt0)/CLOCKS_PER_SEC<<" "<<(dt2-dt1)/CLOCKS_PER_SEC<<" "<<(dt3-dt2)/CLOCKS_PER_SEC<<std::endl;
	
	dfvalue_prev=dfvalue;
	//now we test different projection to estimate the best permutation
	if (bbest_path_proj)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);
		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,NULL,bgreedy);
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};
	
	if (bbest_path_blast_proj)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);
		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,gm_ldh,bgreedy);
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};
	
	if (bbest_path_greedy)
	{
		//permuation projection
		gsl_matrix_transpose_memcpy(gm_P_bp_temp,gm_P);
		gsl_matrix_scale(gm_P_bp_temp,-10000);
		gsl_matrix_hungarian(gm_P_bp_temp,gm_P_prev,NULL,NULL,false,gm_ldh,true);
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P_prev,gm_temp,true);
		if (df_bp_new<fbest_path)
		{
			fbest_path=df_bp_new;
			gsl_matrix_memcpy(gm_P_bp,gm_P_prev);
		};
	};
	
	if (bbest_path_blast_greedy)
	{
	};
	gsl_matrix_memcpy(gm_P_prev,gm_P);
	};//end Frank-Wolfe cycle
	
	//permuation projection
	match_result mres;
	mres.gm_P_exact=gsl_matrix_alloc(N,N);
	gsl_matrix_memcpy(mres.gm_P_exact,gm_P);
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_exact",pdebug.strvalue);
	gsl_matrix_transpose_memcpy(gm_dP,gm_P);
	gsl_matrix_scale(gm_dP,-10000);
	gsl_matrix_hungarian(gm_dP,gm_P,NULL,NULL,false,(bblast_match_end?gm_ldh:NULL),bgreedy);
	if (bbest_path){
		double df_bp_new=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
		if (df_bp_new>fbest_path)
		{
			*gout<<"Best path solution is used"<<std::endl;
			gsl_matrix_memcpy(gm_P,gm_P_bp);
		};	
		gsl_matrix_free(gm_P_bp);
		gsl_matrix_free(gm_P_bp_temp);
	};
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P_projected",pdebug.strvalue); 
	
	//memory deallocation
	gsl_matrix_free(gm_temp);
	gsl_matrix_free(gm_dP);
	gsl_matrix_free(gm_P_prev);
	gsl_matrix_free(gm_ddP_prev);
	gsl_vector_free(gv_C);
	gsl_vector_free(gv_temp);
	gsl_vector_free(gv_debug_trace);
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_vector_free(gv_col_inc);
	gsl_permutation_free(gp_sol);
	
	mres.gm_P=gm_P;
	
	//initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,gm_P,cscore_matrix));
	//other output parameters
	mres.dres=mres.vd_trace.at(1);
	//transpose matrix save
	mres.gm_P=gm_P;
	return mres;
}

//optimized version of gradient calculations: tensor product tricks
void algorithm_qcv::qcv_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,gsl_matrix * gm_temp)
{
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix *gm_1=gsl_matrix_alloc(N,N);
	gsl_matrix_view gmv_grad=gsl_matrix_view_array(gv_grad->data,N,N);
	gsl_matrix_view gmv_P=gsl_matrix_view_array(gv_P->data,N,N);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,&gmv_P.matrix,0,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_temp,0,&gmv_grad.matrix);
	//I gsl_blas_dgemm(CblasTrans,CblasNoTrans,1,gm_Ah_d,gm_Ah_d,0,gm_temp);
	//I gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,&gmv_P.matrix,gm_temp,0,gm_1);
	gsl_matrix_transpose(&gmv_P.matrix);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,gm_temp,0,gm_1);
	gsl_matrix_transpose(gm_1);

	gsl_matrix_add(&gmv_grad.matrix,gm_1);
	
	//first term
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	//gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,&gmv_P.matrix,gm_Ah_d,0,gm_temp);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);
	gsl_matrix_transpose(gm_temp);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
	//second term
	if (bnosymm) gsl_matrix_transpose(gm_Ag_d);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ah_d,&gmv_P.matrix,0,gm_temp);
	gsl_matrix_transpose(gm_temp);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
	
	gsl_matrix_transpose(&gmv_P.matrix);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	//III gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasNoTrans,CblasNoTrans,-2,gm_Ag_d,gm_temp,1,&gmv_grad.matrix);
	gsl_matrix_free(gm_1);
	
	gsl_matrix_transpose(&gmv_grad.matrix);
	gsl_matrix_scale(&gmv_grad.matrix,1-dalpha_ldh);//label cost function scaling
	gsl_vector_scale(gv_grad,1/df_norm);//normalization
	
}

