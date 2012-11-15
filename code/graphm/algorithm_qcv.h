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
#ifndef ALGORITHM_QCV_H
#define ALGORITHM_QCV_H

#include "algorithm.h"

/**
Simple graph matching algorithm: match vertices with identity permutation matrix

	@author Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/
class algorithm_qcv : public algorithm
{
public:
    virtual match_result match(graph &g,graph &h,gsl_matrix* gm_P_i=NULL, gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);
    //optimized version of gradient calculations: tensor product tricks
    void qcv_gradient_opt(gsl_matrix *gm_Ag_d,gsl_matrix *gm_Ah_d,gsl_vector* gv_P, gsl_vector* gv_grad,gsl_matrix * gm_temp);

};

#endif
