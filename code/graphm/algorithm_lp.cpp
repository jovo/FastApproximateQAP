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
#ifdef LP_INCLUDE
#include "algorithm_lp.h"

match_result algorithm_lp::match(graph& g, graph& h,gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh)
{
	if (bverbose)
		*gout<<"Linear programming matching"<<std::endl;	
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	//linear program parameters
	int* ia = new int[N*N*(N+N-1)+N*N+N*N+N*N+1];
	int* ja = new int[N*N*(N+N-1)+N*N+N*N+N*N+1];
	double* ar=new double [N*N*(N+N-1)+N*N+N*N+N*N+1];
	gsl_matrix* gm_Agh=gsl_matrix_alloc(N*N,N*N);
	gsl_matrix_set_zero(gm_Agh);
	long icounter=1;
	//direct lp parameters setting
	for (int i=0;i<N;i++)
		for (int j=0;j<N;j++){
			for (int k=0;k<N;k++){
				 gsl_matrix_set(gm_Agh,i*N+j,j*N+k,gsl_matrix_get(gm_Ag_d,i,k));
				 ia[icounter]=i*N+j+1;ja[icounter]=j*N+k+1;ar[icounter]=gsl_matrix_get(gm_Ag_d,i,k);
				 icounter++;
			};	
			for (int k=0;k<N;k++){
			   //sometimes we have to change yet allocated constaints
			   if (k!=j){
			   ia[icounter]=i*N+j+1;ja[icounter]=i+k*N+1;ar[icounter]=gsl_matrix_get(gm_Agh,i*N+j,i+k*N)-gsl_matrix_get(gm_Ah_d,k,j);
			   icounter++;
			   }
			   else   { ar[icounter-N]=ar[icounter-N]-gsl_matrix_get(gm_Ah_d,k,j); };
			   gsl_matrix_set(gm_Agh,i*N+j,i+k*N,gsl_matrix_get(gm_Agh,i*N+j,i+k*N)-gsl_matrix_get(gm_Ah_d,k,j));
			};
			};
	//ar filtration
	for (int i=0;i<icounter-1;i++)
		if (abs(ar[i+1])<1e-20) ar[i+1]=1e-20;
	//T-S constraints
	for (int i=0;i<N*N;i++)
	{
	  ia[icounter]=i+1;
	  ja[icounter]=N*N+i+1;
	  ar[icounter]=-1;
	  icounter++;
	};
	for (int i=0;i<N*N;i++)
	{
	  ia[icounter]=i+1;
	  ja[icounter]=2*N*N+i+1;
	  ar[icounter]=1;
	  icounter++;
	};
	//probability constraints
	for (int i=0;i<N;i++)
	 for (int j=0;j<N;j++)
	  {
		 ia[icounter]=N*N+i+1;
		 ja[icounter]=N*i+j+1;
		 ar[icounter]=1;
		 icounter++;
	  };
	//glpk package
	LPX* lp;
	lp=lpx_create_prob();	
	lpx_set_prob_name(lp,"glpk");
	lpx_add_rows(lp,N*N+N);
	for (int i=0;i<N*N;i++)
		lpx_set_row_bnds(lp,i+1,LPX_FX,0,0);//problem objective function consraints
	for (int i=0;i<N;i++)
		lpx_set_row_bnds(lp,i+1+N*N,LPX_FX,1,1);//doubly stochastic matrices
	lpx_add_cols(lp,3*N*N);
	for (int i=0;i<3*N*N;i++)
		lpx_set_col_bnds(lp,i+1,LPX_LO,0,0);
	//objective function
	lpx_set_obj_dir(lp,LPX_MIN);
	//the label cost matrix introduce the coefficients for the first part of objective function
	if (dalpha_ldh>0)
	{		
	     for (int i=0;i<N*N;i++)
			lpx_set_obj_coef(lp,i+1,dalpha_ldh*gm_ldh->data[i]);
	     for (int i=0;i<2*N*N;i++)
			lpx_set_obj_coef(lp,N*N+i+1,1-dalpha_ldh);
	}
	else{
		for (int i=0;i<N*N;i++)
			lpx_set_obj_coef(lp,i+1,0);
		for (int i=0;i<2*N*N;i++)
			lpx_set_obj_coef(lp,N*N+i+1,1);
	}
	lpx_load_matrix(lp,N*N*(N+N-1)+N*N+N*N+N*N,ia,ja,ar);
	//problem solution
	lpx_simplex(lp);
	//result save
	gsl_matrix* C=gsl_matrix_alloc(N,N);
	for (int i=0;i<N;i++)
	 for (int j=0;j<N;j++)
		gsl_matrix_set(C,i,j,lpx_get_col_prim(lp,i*N+j+1));
	lpx_delete_prob(lp);
	if (pdebug.ivalue) gsl_matrix_printout(C,"C=lp solution",pdebug.strvalue); 
	match_result mres;
	mres.gm_P_exact=gsl_matrix_alloc(N,N);
	gsl_matrix_transpose_memcpy(mres.gm_P_exact,C);
	//solution projection on the set of permutation matrices
	double dscale_factor =gsl_matrix_max_abs(C);
	dscale_factor=(dscale_factor>EPSILON)?dscale_factor:EPSILON;
	dscale_factor=10000/dscale_factor;
	gsl_matrix_scale(C,-dscale_factor);
	//gsl_matrix_transpose(C);
	if (pdebug.ivalue) gsl_matrix_printout(C,"scale(C)",pdebug.strvalue); 
	gsl_matrix* gm_P=gsl_matrix_hungarian(C);
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue); 
	mres.gm_P=gm_P;
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_matrix_free(gm_Agh);
	delete[] ja;
	delete[] ia;
	delete[] ar;

	//initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,gm_P,cscore_matrix));
	//other output parameters
	mres.dres=mres.vd_trace.at(1);
	mres.inum_iteration=2;
	//transpose matrix save
	mres.gm_P=gm_P;
	return mres;
}

#endif