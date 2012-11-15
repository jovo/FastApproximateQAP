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
#ifndef EXPERIMENT_H
#define EXPERIMENT_H

#include "rpc.h"
#include "math.h"
#include "graph.h"
#include "algorithm.h"
#include "algorithm_umeyama.h"
#include "algorithm_iden.h"
#include "algorithm_rank.h"
#ifdef LP_INCLUDE
#include "algorithm_lp.h"
#endif
#include "algorithm_qcv.h"
#include "algorithm_unif.h"
#include "algorithm_fsol.h"
#include "algorithm_rand.h"
#include "algorithm_path.h"
#include "algorithm_ca.h"
#include "algorithm_ext.h"
#include <strstream>

/**
Experiment class. It implements all extern routines for graph matching experiments

	@author Mikhail Zaslavskiy
*/


class experiment : public rpc
{
public:
    experiment(std::string fconfig="config.txt"):rpc(fconfig){};
    void run_experiment();
    ~experiment();
    void printout(std::string fname_out,std::string sformat);
    void printout(std::string fname_out);
    void printout();
    algorithm* get_algorithm(std::string salgo);
protected:
    std::vector<match_result> v_mres;
    void synchronize(graph &g, graph &h,gsl_matrix** gm_ldh=NULL);

};

#endif
