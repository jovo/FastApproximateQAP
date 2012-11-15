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
#ifndef ALGORITHM_LP_H
#define ALGORITHM_LP_H

#include "algorithm.h"
extern "C"{
#include "glpk.h"

}

/**
LP relaxation:  H. A. Almohamad and S. 0. Duffuaa. A Linear Programming Approach for the Weighted Graph Matching Problem. 

	@author Mikhail Zaslavskiy <mikhail.zaslavskiy@ensmp.fr>
*/
class algorithm_lp : public algorithm
{
public:
    virtual match_result match(graph &g,graph &h,gsl_matrix* gm_P_i=NULL,gsl_matrix* gm_ldh=NULL,double dalpha_ldh=-1);

};

#endif
