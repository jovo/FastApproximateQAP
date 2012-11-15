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
#include "algorithm_iden.h"

match_result algorithm_iden::match(graph& g, graph& h,gsl_matrix* gm_P_i, gsl_matrix* gm_ldh,double dalpha_ldh)
{
	if (bverbose)
		*gout<<"Identity matching"<<std::endl;	
	//some duplicate variables
	match_result mres;
	mres.gm_P=gsl_matrix_alloc(N,N);
	gsl_matrix_set_identity(mres.gm_P);
	mres.gm_P_exact=NULL;
        //initial score
	mres.vd_trace.push_back(graph_dist(g,h,cscore_matrix));
	//final score
	mres.vd_trace.push_back(graph_dist(g,h,mres.gm_P,cscore_matrix));
	mres.dres=mres.vd_trace.at(1);
	mres.inum_iteration=2;
	return mres;
}

