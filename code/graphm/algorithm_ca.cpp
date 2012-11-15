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
#include "algorithm_ca.h"

match_result algorithm_ca::match(graph& g, graph& h,gsl_matrix* gm_P_i,gsl_matrix* _gm_ldh,double dalpha_ldh)
{
	double dhung_max=get_param_d("hungarian_max");
	bool bblast_match_end=(get_param_i("blast_match_proj")==1);	
	bool bblast_match=(get_param_i("blast_match")==1);
	double dfw_xeps=get_param_d("algo_fw_xeps");
	double bgreedy=(get_param_i("hungarian_greedy")==1);
	double dfvalue;
	if (bverbose)
		*gout<<"C-adaptation strategy"<<std::endl;	
	//some duplicate variables
	gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);
	gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);
	gsl_matrix * gm_C_a=gsl_matrix_alloc(N,N);
	gsl_matrix * gm_C_up=gsl_matrix_alloc(N,N);
	gsl_matrix * gm_P=gsl_matrix_alloc(N,N);
	gsl_matrix * gm_temp=gsl_matrix_alloc(N,N);
	//initialization
	gsl_matrix_memcpy(gm_C_a,gm_P_i);
	//cycle over linear combination of linear and non linear terms
	double astep=1e-30;bool bcont=true;
	if (dalpha_ldh>0)
		gsl_matrix_scale(gm_ldh,dalpha_ldh);

	while (bcont){
	bcont=(astep<1);
	//optimal permutation
	gsl_matrix_scale(gm_C_a,-dhung_max);
	gsl_matrix_transpose(gm_C_a);
	gsl_matrix_hungarian(gm_C_a,gm_P,NULL,NULL,false,(bblast_match?gm_ldh:NULL),bgreedy);
	gsl_matrix_transpose(gm_C_a);
	gsl_matrix_scale(gm_C_a,-1.0/dhung_max);
	//C-matrix update
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_temp,gm_Ah_d,0,gm_C_up);
	gsl_matrix_scale(gm_C_up,1-dalpha_ldh);
	if (dalpha_ldh>0) gsl_matrix_add(gm_C_up,gm_ldh);
	gsl_matrix_scale(gm_C_a,1-astep);
	gsl_matrix_scale(gm_C_up,astep);
	gsl_matrix_add(gm_C_a,gm_C_up);
	dfvalue=f_qcv(gm_Ag_d,gm_Ah_d,gm_P,gm_temp,true);
	if (bverbose) *gout<<"astep="<<astep<<", f="<<dfvalue<<std::endl;
	astep+=0.01;//0.01*astep;
	astep=(astep<1)?astep:1;
	};
	gsl_matrix_free(gm_Ag_d);
	gsl_matrix_free(gm_Ah_d);
	gsl_matrix_free(gm_temp);
	gsl_matrix_free(gm_C_up);
	gsl_matrix_free(gm_C_a);
	if (dalpha_ldh>0)
		gsl_matrix_scale(gm_ldh,1.0/dalpha_ldh);
	
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

