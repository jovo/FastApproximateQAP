/***************************************************************************
 *   Copyright (C) 2007 by Mikhail Zaslavskiy   *
 *      *
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
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include "rpc.h"
#include <math.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <vector>

/**
Graph class, it implements graph matching algorithms

	@author Mikhail Zaslavskiy
*/
class graph : public rpc
{
public:
    //class constructor
    graph(std::string fconfig="config.txt");
    graph(graph & _graph);
    const graph & operator =(graph &);

    //load graph routine
    int load_graph(std::string fgraph,char ftype='A',char cformat='D',std::string fvertexlst_name="");
    //adjacency matrix control
    const gsl_matrix*  get_adjmatrix(){return gm_A;};
    int set_adjmatrix(const gsl_matrix* _gm_A);
    gsl_matrix* get_descmatrix(char dmt);
    int printout(std::string fname_out);
    void printdot(std::string fname_out, gsl_matrix* =NULL);
    ~graph();
    int add_dummy_nodes(int id);
    long long getN();
protected:
    gsl_matrix* gm_A;//adjacency matrix
    long long N;
};
double gsl_matrix_norm(const gsl_matrix* gm,double p);
int gsl_matrix_abs(gsl_matrix* gm);
double gsl_matrix_min(gsl_matrix* gm_A,double ic);//minimum greater than ic
double gsl_vector_sum(gsl_vector* gv);
int gsl_matrix_printout(const gsl_matrix * gm,std::string,std::string);
int gsl_matrix_printout(gsl_vector * gv,std::string,std::string);
int gsl_matrix_printout(gsl_permutation * gv,std::string,std::string);
int gsl_matrix_printout(std::string sout,std::string );
double abs(double);
void gsl_matrix_sum(gsl_matrix* A,int idim,gsl_vector* gv_res);
double gsl_matrix_sum(gsl_matrix* A);
double gsl_matrix_max_abs(gsl_matrix* A);
double gsl_matrix_min_abs(gsl_matrix* gm);
//global min in [0,1]^2 of a quadratic function
void qglob_min(double a_11,double a_12,double a_22,double a_1,double a_2,double &dalpha_1,double &dalpha_2);
#endif
