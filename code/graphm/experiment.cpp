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
#include "experiment.h"
#include "graph.h"

void experiment::run_experiment()
{
	parameter pdebug = get_param("debugprint");
	parameter pdebug_f = get_param("debugprint_file");
	pdebug.strvalue=pdebug_f.strvalue;
	std::string sverbfile=get_param_s("verbose_file");
	bool bverbose=(get_param_i("verbose_mode")==1);
	std::ofstream fverbose;
	if (sverbfile.compare("cout")==0)
		gout=&std::cout;
	else
	{
		fverbose.open(sverbfile.c_str());
	       	gout=&fverbose;
	};

	//used algorithms
	std::string sexp_type=get_param_s("exp_type");
	std::string salgo_match=get_param_s("algo");
	std::istrstream istr(salgo_match.c_str());
	std::vector<std::string>v_salgo_match;
	std::string stemp;
	while (istr>>stemp)
		v_salgo_match.push_back(stemp);
	//used initialization algorithms
	std::string salgo_init=get_param_s("algo_init_sol");
	std::istrstream istrinit(salgo_init.c_str());
	std::vector<std::string>v_salgo_init;
	while (istrinit>>stemp)
		v_salgo_init.push_back(stemp);
	if (v_salgo_init.size()!=v_salgo_match.size())
		{
			printf("Error: algo and algo_init_sol do not have the same size\n");
			exit(0);
		};

	double dalpha_ldh=get_param_d("alpha_ldh");
	if ((dalpha_ldh<0) ||(dalpha_ldh>1))
	{
		printf("Error:alpha_ldh should be between 0 and 1 \n");
		exit(0);
	};

	graph g(get_config());
	graph h(get_config());
	char cformat1='D',cformat2='D';
	if (bverbose) *gout<<"Data loading"<<std::endl;
	g.load_graph(get_param_s("graph_1"),'A',cformat1);
	h.load_graph(get_param_s("graph_2"),'A',cformat2);
	if (get_param_i("graph_dot_print")==1){ g.printdot("g.dot");h.printdot("h.dot");};
	// Matrix of local similarities between graph vertices
	std::string strfile=get_param_s("C_matrix");
	int N1 =g.get_adjmatrix()->size1;
	int N2 =h.get_adjmatrix()->size1;
 	gsl_matrix* gm_ldh=gsl_matrix_alloc(N1,N2);
	FILE *f=fopen(strfile.c_str(),"r");
	if (f!=NULL){
		gsl_set_error_handler_off();
		int ierror=gsl_matrix_fscanf(f,gm_ldh);
		fclose(f);
		gsl_set_error_handler(NULL);
		if (ierror!=0){ printf("Error: C_matrix is not correctly defined \n"); exit(0);};
		}
	else
	{
		if (bverbose) *gout<<"C matrix is set to constant"<<std::endl;
                dalpha_ldh=0;
		gsl_matrix_set_all(gm_ldh,1.0/(N1*N2));
	};
	//inverse the sign if C is a distance matrix
	if (get_param_i("C_matrix_dist")==1)
		gsl_matrix_scale(gm_ldh,-1);

	if (bverbose) *gout<<"Graph synhronization"<<std::endl;
	synchronize(g,h,&gm_ldh);

	//Cycle over all algoirthms
	for (int a=0;a<v_salgo_match.size();a++)
		{
			//initialization
			algorithm* algo_i=get_algorithm(v_salgo_init[a]);
			algo_i->read_config(get_config_string());
			algo_i->set_ldhmatrix(gm_ldh);
			match_result mres_i=algo_i->gmatch(g,h,NULL,NULL,dalpha_ldh);
			delete algo_i;
			//main algorithm
			algorithm* algo=get_algorithm(v_salgo_match[a]);
			algo->read_config(get_config_string());
			algo->set_ldhmatrix(gm_ldh);
			match_result mres_a=algo->gmatch(g,h,(mres_i.gm_P_exact!=NULL)?mres_i.gm_P_exact:mres_i.gm_P, NULL,dalpha_ldh);
			mres_a.vd_trace.clear();
			mres_a.salgo=v_salgo_match[a];
			v_mres.push_back(mres_a);
			delete algo;
		};
	gsl_matrix_free(gm_ldh);
	printout();
}
experiment::~experiment()
{
}

//graph size synhronization by dummy nodes adding
void experiment::synchronize(graph& g, graph& h, gsl_matrix** pgm_ldh)
{
	int iadd_r=0,iadd_c=0;
	if (g.get_adjmatrix()->size1>h.get_adjmatrix()->size1)
		{
			iadd_c=g.get_adjmatrix()->size1-h.get_adjmatrix()->size1;

		};
	if (g.get_adjmatrix()->size1<h.get_adjmatrix()->size1)
		{
			iadd_r=-g.get_adjmatrix()->size1+h.get_adjmatrix()->size1;
		};
	int idummy_nodes= get_param_i("dummy_nodes");
	switch(idummy_nodes){
	case 0://just complete the smallest graph
		break;
	case 1://double the matrix size
		iadd_r=h.get_adjmatrix()->size1;
		iadd_c=g.get_adjmatrix()->size1;
		break;
	};
	h.add_dummy_nodes(iadd_c);
	g.add_dummy_nodes(iadd_r);

	if (pgm_ldh!=NULL){
	gsl_matrix* gm_C_new=gsl_matrix_alloc((*pgm_ldh)->size1+iadd_r,(*pgm_ldh)->size2+iadd_c);
	gsl_matrix_view gmv_C_new=gsl_matrix_submatrix(gm_C_new,0,0,(*pgm_ldh)->size1,(*pgm_ldh)->size2);
	double dmin_C_blast=gsl_matrix_min((*pgm_ldh));
	double dmax_C_blast=gsl_matrix_max((*pgm_ldh));
	double dccoef= get_param_d("dummy_nodes_c_coef");
	dmin_C_blast=(1-dccoef)*dmin_C_blast+dccoef*dmax_C_blast;
	gsl_matrix_set_all(gm_C_new,dmin_C_blast);
	gsl_matrix_memcpy(&gmv_C_new.matrix,(*pgm_ldh));
	gsl_matrix_free((*pgm_ldh));
	*pgm_ldh=gm_C_new;
	};
}

/*!
    \fn experiment::printout(std::string fname_out)
 */
void experiment::printout()
{
	std::string sfile=get_param_s("exp_out_file");
	printout(sfile);
}
void experiment::printout(std::string fname_out)
{
	std::string sformat=get_param_s("exp_out_format");
	printout(fname_out,sformat);
}
void experiment::printout(std::string fname_out,std::string sformat)
{

	std::ofstream fout(fname_out.c_str(),std::ios::app);
	fout<<"***********************************"<<std::endl;
	fout<<"Experiment parameters:"<<std::endl;
	if (sformat.find("Parameters")!=-1)
	{
	 	rpc::printout(fname_out);
	};
	fout<<"Experiment results:"<<std::endl;
	if (sformat.find("Compact")!=-1)
	{
		fout<<"\t Alpha\t\t";
		for (int a=0;a<v_mres.size();a++)
			fout<<v_mres.at(a).salgo<<"\t\t";
		fout.setf(std::ios::scientific);
		fout<<std::endl<<"Gdist \t"<<get_param_d("alpha_ldh");
		for (int a=0;a<v_mres.size();a++)
			{
				fout<<"\t"<<v_mres.at(a).dres;
			};
		fout<<std::endl<<"F_perm \t"<<get_param_d("alpha_ldh");
		for (int a=0;a<v_mres.size();a++)
			{
				fout<<"\t"<<v_mres.at(a).dfvalue;
			};
		fout<<std::endl<<"F_exact\t"<<get_param_d("alpha_ldh");
		for (int a=0;a<v_mres.size();a++)
			{
				fout<<"\t"<<v_mres.at(a).dfvalue_exact;
			};
		fout<<std::endl<<"Time: \t\t";
		for (int a=0;a<v_mres.size();a++)
			{
				fout<<"\t"<<v_mres.at(a).dtime;
			};
		fout<<std::endl;

	};
	if (sformat.find("Permutation")!=-1)
	{
		fout<<"Permutations:"<<std::endl;
		for (int a=0;a<v_mres.size();a++)
			fout<<v_mres.at(a).salgo<<" ";
		fout.precision(2);
		fout<<std::endl;
		for (int i=0;i<v_mres.at(0).gm_P->size1;i++)
		{
			for (int a=0;a<v_mres.size();a++)
			{

					int ipout=-1;
					for (int j=0;j<v_mres.at(a).gm_P->size2;j++)
						if(v_mres.at(a).gm_P->data[i*v_mres.at(a).gm_P->size1+j]==1)
						ipout=j;
					fout<<ipout+1<<" ";
			};
			fout<<std::endl;
		};

	};
}

/*!
    \fn experiment::get_algorithm(std::string salgo)
    Dynamic association of graph matching algorithm
 */

algorithm* experiment::get_algorithm(std::string salgo)
{
   if (salgo.compare("U")==0){ return new algorithm_umeyama;};
   if (salgo.compare("UMEY")==0){ return new algorithm_umeyama;};
   
   if (salgo.compare("I")==0){ return new algorithm_iden;};
   if (salgo.compare("IDEN")==0){ return new algorithm_iden;};
   
   #ifdef LP_INCLUDE
   if (salgo.compare("LP")==0){ return new algorithm_lp;};
   if (salgo.compare("L")==0){ return new algorithm_lp;};
   #endif
   if (salgo.compare("RANK")==0){ return new algorithm_rank;};
   if (salgo.compare("ISORANK")==0){ return new algorithm_rank;};
   if (salgo.compare("isorank")==0){ return new algorithm_rank;};
   if (salgo.compare("IsoRank")==0){ return new algorithm_rank;};
   if (salgo.compare("R")==0){ return new algorithm_rank;};
   
   if (salgo.compare("QCV")==0){ return new algorithm_qcv;};
   if (salgo.compare("Q")==0){ return new algorithm_qcv;};
   
   if (salgo.compare("P")==0){ return new algorithm_path;};
   if (salgo.compare("PATH")==0){ return new algorithm_path;};
   
   if (salgo.compare("rand")==0){ return new algorithm_rand;};
   if (salgo.compare("RAND")==0){ return new algorithm_rand;};
   if (salgo.compare("fsol")==0){ return new algorithm_fsol;};
   
   if (salgo.compare("unif")==0){ return new algorithm_unif;};
   if (salgo.compare("UNIF")==0){ return new algorithm_unif;};
   if (salgo.compare("ca")==0){ return new algorithm_ca;};
   if (salgo.compare("CA")==0){ return new algorithm_ca;};
   if (salgo.compare("ga")==0){ return new algorithm_ca;};
   if (salgo.compare("GA")==0){ return new algorithm_ca;};
   
   if (salgo.compare("EXT")==0){ return new algorithm_NEW;};
   printf("Error: graph matching algorithm is not selected");
   exit(0);
}
