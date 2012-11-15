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
#ifndef ALGORITHM_PATH_H
#define ALGORITHM_PATH_H
#include "algorithm.h"
#include "algorithm_qcv.h"

/**
A path following graph matching algorithm

	@author Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/
class algorithm_ssp : public algorithm_qcv
{
public:
    virtual match_result match(graph &g,graph &h,gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);
    virtual void qcc_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dmult,gsl_matrix * gm_temp);
    void qcvqcc_gradient(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_vector* gv_P, gsl_vector* gv_grad,double dlambda,gsl_vector * gv_temp);
    double f_qcvqcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix *gm_Lg_d,gsl_matrix *gm_Lh_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,double dlambda,gsl_matrix * gm_temp,gsl_matrix *gm_temp2);
    //concave function value
     virtual double f_qcc(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_matrix* gm_Delta,gsl_matrix* gm_P,gsl_matrix * gm_temp,gsl_matrix *gm_temp2);

};

#endif
