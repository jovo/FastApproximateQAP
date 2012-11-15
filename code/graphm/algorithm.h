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
#ifndef ALGORITHM_H
#define ALGORITHM_H
#define EPSILON 1e-100
#include "rpc.h"
#include "graph.h"
#include <math.h>
#include "hungarian.h"

#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>

#include <vector>
#include <iostream>



/**
Class of graph matching results
*/
class match_result
{
public:
	match_result(){gm_P=NULL;gm_P_exact=NULL;salgo="";}
	std::vector<double> vd_trace;
	int inum_iteration;
	double dres;
	gsl_matrix* gm_P;
	gsl_matrix* gm_P_exact;
	double dtime;
	double dfvalue;
	double dfvalue_exact;
	std::string salgo;
	~match_result(){};

};
/**
Parent class for all graph matching algorithms

	@author Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/
class algorithm : public rpc
{
public:
     algorithm(std::string );
     algorithm();
     match_result gmatch(graph& g, graph& h,gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);//common stuff,
     virtual match_result match(graph& g, graph& h, gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1)=0;//particular method implementation
     double graph_dist(graph &g,graph &h,gsl_matrix* gm_P,char cscore_matrix);
     double graph_dist(graph &g, graph &h,char cscore_matrix);
      ~algorithm();

     const gsl_matrix* get_ldhmatrix(){return gm_ldh;};
     void  set_ldhmatrix(const gsl_matrix* _gm_A);

protected:
    gsl_matrix *gm_ldh;
    double dalpha_ldh;
    void update_C_hungarian(gsl_matrix* gm_C,int isign=1, bool bback=false);
    double f_qcv(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_P,gsl_matrix * gm_temp,bool bqcv=false);
    char cdesc_matrix,cscore_matrix;
    parameter pdebug,pdebug_f;
    bool bverbose;
    std::string sverbfile;
    std::ofstream fverbose;
    long long N;
    double df_norm;
    bool bnosymm;


};

#endif
