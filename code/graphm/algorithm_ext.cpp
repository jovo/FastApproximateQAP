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
#include "algorithm_ext.h"

match_result algorithm_NEW::match(graph& g, graph& h,gsl_matrix* gm_P_i,gsl_matrix* gm_ldh,double dalpha_ldh)
{
if (bverbose) *gout<<"The best matching algorithm"<<std::endl;  
match_result mres; //class with results
gsl_matrix* gm_Ag_d=g.get_descmatrix(cdesc_matrix);//get the adjacency matrix of graph g
gsl_matrix* gm_Ah_d=h.get_descmatrix(cdesc_matrix);//get the adjacency matrix of graph h
//the similarity matrix C is defined in the algorithm class memeber gm_ldh
//dalpha_ldh is corresponding is corresponding to the linear combination coefficent alpha

gsl_matrix* P=gsl_matrix_alloc(N,N);
//YOUR OPERATIONS  WITH MATRICES, RESULT IS A PERMUTATION MATRIX P

//do not forget to release the memory
gsl_matrix_free(gm_Ag_d);
gsl_matrix_free(gm_Ah_d);

mres.gm_P=P;//save the solution
mres.gm_P_exact=NULL; //you can save here the matrix which was used as an approximation for P


mres.dres=graph_dist(g,h,mres.gm_P,cscore_matrix);// distance between graph adjacency matrices
return mres;
};