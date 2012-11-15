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
#include "algorithm.h"

algorithm::algorithm(std::string fconfig)
 : rpc(fconfig)
{
gm_ldh=NULL;dalpha_ldh=0;bnosymm=false;
}
algorithm::algorithm()
 : rpc()
{
gm_ldh=NULL;dalpha_ldh=0;bnosymm=false;df_norm=0;N=0;cdesc_matrix='A';cscore_matrix='A';
}

//common framework for graph matching algorithm
match_result algorithm::gmatch(graph& g, graph& h,gsl_matrix *gm_P_i,gsl_matrix* _gm_ldh,double _dalpha_ldh)
{
dalpha_ldh=_dalpha_ldh;
if (_gm_ldh!=NULL) set_ldhmatrix(_gm_ldh);
pdebug = get_param("debugprint");
pdebug_f = get_param("debugprint_file");
pdebug.strvalue=pdebug_f.strvalue;
bverbose=(get_param_i("verbose_mode")==1);
sverbfile=get_param_s("verbose_file");
cdesc_matrix=get_param_c("cdesc_matrix");
cscore_matrix=get_param_c("cscore_matrix");
bnosymm=(get_param_i("nosymm_matrix")==1);
if (sverbfile.compare("cout")==0)
	gout=&std::cout;
else { fverbose.open(sverbfile.c_str(),std::ios_base::app);
       	gout=&fverbose;	};
N=g.getN();

gsl_matrix* gm_Ag_d=g.get_descmatrix(cscore_matrix);
gsl_matrix* gm_Ah_d=h.get_descmatrix(cscore_matrix);
gsl_matrix* gm_temp=gsl_matrix_alloc(N,N);
df_norm=pow(gsl_matrix_norm(gm_Ag_d,2),2)+pow(gsl_matrix_norm(gm_Ah_d,2),2);
df_norm=(df_norm>EPSILON)?df_norm:EPSILON;

time_t t1=time(NULL);
match_result mres=match(g,h,gm_P_i,_gm_ldh,_dalpha_ldh);

mres.dfvalue=f_qcv(gm_Ag_d,gm_Ah_d,mres.gm_P,gm_temp,true);
if (mres.gm_P_exact!=NULL)
	mres.dfvalue_exact=f_qcv(gm_Ag_d,gm_Ah_d,mres.gm_P_exact,gm_temp,true);
else
	mres.dfvalue_exact=mres.dfvalue;
gsl_matrix_free(gm_Ag_d);	
gsl_matrix_free(gm_Ah_d);
gsl_matrix_free(gm_temp);

mres.dtime=difftime(time(NULL),t1);
return mres;
}
void algorithm::update_C_hungarian(gsl_matrix* gm_C,int isign, bool bback)
{
	if (!bback)
	{
		if (dalpha_ldh>0){
		gsl_matrix_scale(gm_ldh,isign*dalpha_ldh);
		gsl_matrix_sub(gm_C,gm_ldh);
		gsl_matrix_scale(gm_ldh,isign*(1/dalpha_ldh));
		};
	}
	else
	{
		if (dalpha_ldh>0){
		gsl_matrix_scale(gm_ldh,isign*dalpha_ldh);
		gsl_matrix_add(gm_C,gm_ldh);
		gsl_matrix_scale(gm_ldh,isign*(1/dalpha_ldh));
		};
	}

}

//||G-PHP'||^2_F - graph distance
double algorithm::graph_dist(graph &g,graph &h,gsl_matrix* gm_P,char cscore_matrix)
{
	parameter pdebug = get_param("debugprint");
	parameter pdebug_f = get_param("debugprint_file");
	pdebug.strvalue=pdebug_f.strvalue;
	long long N=g.getN();	
	gsl_matrix* gm_Ag=g.get_descmatrix(cscore_matrix);
	gsl_matrix* gm_At=gsl_matrix_alloc(N,N);
	gsl_matrix* gm_Ah=gsl_matrix_alloc(N,N);
	gsl_matrix_memcpy(gm_Ah,h.get_adjmatrix());
	if (pdebug.ivalue) gsl_matrix_printout(gm_P,"gm_P",pdebug.strvalue);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ah,"gm_Ah",pdebug.strvalue);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_P,gm_Ah,0,gm_At);
	if (pdebug.ivalue) gsl_matrix_printout(gm_At,"gm_At",pdebug.strvalue);
	gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,gm_At,gm_P,0,gm_Ah);
	gsl_matrix_sub(gm_Ag,gm_Ah);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag,"Ag-Ah",pdebug.strvalue);

	double dret=pow(gsl_matrix_norm(gm_Ag,2),2);
	
	gsl_matrix_free(gm_At);
	gsl_matrix_free(gm_Ah);
	gsl_matrix_free(gm_Ag);
	return dret;
}
double algorithm::graph_dist(graph &g, graph &h,char cscore_matrix){
	parameter pdebug = get_param("debugprint");
	parameter pdebug_f = get_param("debugprint_file");
	pdebug.strvalue=pdebug_f.strvalue;
	
	gsl_matrix* gm_Ag=g.get_descmatrix(cscore_matrix);
	gsl_matrix* gm_Ah=h.get_descmatrix(cscore_matrix);
	gsl_matrix_sub(gm_Ag,gm_Ah);
	if (pdebug.ivalue) gsl_matrix_printout(gm_Ag,"Ag-Ah",pdebug.strvalue);
	double dret=pow(gsl_matrix_norm(gm_Ag,2),2);
	gsl_matrix_free(gm_Ah);
	gsl_matrix_free(gm_Ag);
	return dret;
}

void algorithm::set_ldhmatrix(const gsl_matrix* _gm_A)
{
 if (gm_ldh!=NULL)
 	gsl_matrix_free(gm_ldh);
 gm_ldh=NULL;
 if (_gm_A!=NULL){ 
 gm_ldh=gsl_matrix_alloc(_gm_A->size1,_gm_A->size2);
 gsl_matrix_memcpy(gm_ldh,_gm_A);
 gsl_matrix_scale(gm_ldh,1.0/(gsl_matrix_norm(gm_ldh,2)+EPSILON));
 }
}

//convex representation of the objective function ||AP-PA||^2_F
double algorithm::f_qcv(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,bool bqcv)
{
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
	//gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasNoTrans,CblasNoTrans,1,gm_Ag_d,gm_P,0,gm_temp);
	gsl_matrix_transpose(gm_P);
	gsl_matrix_transpose(gm_temp);
        if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,-1,gm_Ah_d,gm_P,1,gm_temp);
	if (bnosymm) gsl_matrix_transpose(gm_Ah_d);
	gsl_matrix_transpose(gm_P);
	//gsl_blas_dgemm_sym(CblasLeft,CblasUpper,CblasTrans,CblasTrans,-1,gm_Ah_d,gm_P,1,gm_temp);
	double dres= pow(gsl_matrix_norm(gm_temp,2),2);
	dres=dres*(1-dalpha_ldh);
	dres/=df_norm;
	if ((bqcv) and (dalpha_ldh>0))//if just QCV optimization then alpha_ldh must be integrated, in QCVQCC and P cases we use Delta for these purposes
	{gsl_vector_view gvv_ldh=gsl_vector_view_array(gm_ldh->data,N*N);
	gsl_vector_view gvv_P=gsl_vector_view_array(gm_P->data,N*N);
	double dtr;
	gsl_blas_ddot(&gvv_ldh.vector,&gvv_P.vector,&dtr);
	dtr*=dalpha_ldh;
	dres-=dtr;};
	//normalization
	return dres;
}

algorithm::~algorithm()
{
	if (gm_ldh!=NULL)
 		gsl_matrix_free(gm_ldh);
}


