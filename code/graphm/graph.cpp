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
#include "graph.h"
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_eigen.h>
#include "hungarian.h"
//using namespace std;

graph::graph(std::string fconfig)
 : rpc(fconfig)
{
	gm_A=NULL;
}
graph::graph(graph & gr):rpc()
{
	gm_A=NULL;
	const gsl_matrix* gm_t=gr.get_adjmatrix();
	set_adjmatrix(gm_t);
}
const graph& graph::operator=(graph& gh){
if (&gh!=this)
{
	set_adjmatrix(gh.get_adjmatrix());
}
}
//graph loading from txt file
int graph::load_graph(std::string fgraph_name,char ftype,char cformat,std::string fvertexlst_name)
{
std::ifstream fgraph(fgraph_name.c_str());
double de; int N;

switch(ftype)
{
case 'A':case 'a':
	//matrix size
	N=0;char buf[1];
	while (fgraph>>de)//.getline(buf,1))
	{	
		N++;
	};
	N=sqrt(N);
	int ierror=1;
	if (N>0){
	fgraph.seekg(0,std::ios::beg);	
	gsl_matrix* gm_A_l=gsl_matrix_alloc(N,N);
	FILE *f=fopen(fgraph_name.c_str(),"r");
	gsl_set_error_handler_off();
	ierror=gsl_matrix_fscanf(f,gm_A_l);
	fclose(f);
	gsl_set_error_handler(NULL);
	set_adjmatrix(gm_A_l);
	gsl_matrix_free(gm_A_l);
	};
	if (ierror!=0){ printf("Error: graph adjacency matrix is not correctly defined \n"); exit(0);}
	break;
 };
return 1;
}

gsl_matrix* graph::get_descmatrix(char dmt)
{
gsl_matrix* gm_ret=gsl_matrix_alloc(N,N);
int ip;
switch(dmt)
{
        case 'A'://adjacency matrix
                gsl_matrix_memcpy(gm_ret,gm_A);

        break;
        case 'L'://laplacian matrix
                gsl_matrix_memcpy(gm_ret,gm_A);
                gsl_vector* gv_ones = gsl_vector_alloc(N);
                gsl_vector* gv_res = gsl_vector_alloc(N);
                gsl_vector_set_all(gv_ones,1);
                gsl_blas_dgemv(CblasNoTrans,1,gm_ret,gv_ones,0,gv_res);
                for (int i=0;i<N;i++)
                        gsl_matrix_set(gm_ret,i,i,gsl_vector_get(gv_res,i));
                gsl_vector_free(gv_res);
                gsl_vector_free(gv_ones);
        break;
};
return gm_ret;
}

//dummy nodes adding 
int graph::add_dummy_nodes(int id)
{
	int Nn=N+id;
	double ddummy_nodes_fill=get_param_d("dummy_nodes_fill");
	double dmin=gsl_matrix_min(gm_A);
	double dmax=gsl_matrix_max(gm_A);
	dmin=(1-ddummy_nodes_fill)*dmin+ddummy_nodes_fill*dmax;
 	gsl_matrix* gm_A_new=gsl_matrix_alloc(N+id,N+id);
	gsl_matrix_view gmv_A_new=gsl_matrix_submatrix(gm_A_new,0,0,N,N);
	gsl_matrix_set_all(gm_A_new,dmin);
	gsl_matrix_memcpy(&gmv_A_new.matrix,gm_A);
	set_adjmatrix(gm_A_new);
	gsl_matrix_free(gm_A_new);
	return 0;
};

int graph::set_adjmatrix(const gsl_matrix* _gm_A)
{
 if (gm_A!=NULL)
 	gsl_matrix_free(gm_A);
 N=0;gm_A=NULL;
 if (_gm_A!=NULL){ 
 gm_A=gsl_matrix_alloc(_gm_A->size1,_gm_A->size2);
 gsl_matrix_memcpy(gm_A,_gm_A);
 N=gm_A->size1;}
 return 0;
} 


int graph::printout(std::string fname_out)
{
rpc::printout(fname_out);
std::ofstream fout(fname_out.c_str(),std::ios::app);
fout<<"Adjacency Matrix:"<<std::endl;
for (int i=0;i<gm_A->size1;i++)
{
	for (int j=0;j<gm_A->size2;j++)
		fout<<" "<<gsl_matrix_get(gm_A,i,j);
	fout<<std::endl;
}
}

void graph::printdot(std::string fname_out,gsl_matrix* gm_P)
{
	gsl_matrix* gm_temp;
	int *gper;
	if (gm_P!=NULL)
		{
			gm_temp=gsl_matrix_alloc(N,N);
			gsl_matrix* gm_temp2=gsl_matrix_alloc(N,N);
			gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1,gm_P,gm_A,0,gm_temp2);
			gsl_blas_dgemm(CblasNoTrans,CblasTrans,1,gm_temp2,gm_P,0,gm_temp);
			gsl_matrix_free(gm_temp2);
			gper=new int[N];
			for (int i=0;i<N;i++)
				for (int j=0;j<N;j++)
					if (gm_P->data[i*N+j]==1) gper[i]=j;
		}
	else gm_temp=gm_A;
	char c='"';
   	std::ofstream fout(fname_out.c_str());
	fout<<"graph "<<c<< "MyGraph"<<c<<"{"<<std::endl;
	double dcol_step=0.999/pow(1.001*gm_temp->size1,0.33);
	double dr=0,dg=0,db=0;
	for (int i=0;i<gm_temp->size1;i++)
	{
		dr+=dcol_step;
		if (dr>1) {dr=0;dg+=dcol_step;};
		if (dg>1) {dg=0;db+=dcol_step;};
		if (i<gm_temp->size1/3)
		{
			dr=(double)i*3/gm_temp->size1;
			dg=db=0;
		}
		else
		if (i<2*gm_temp->size1/3)
		{
			dg=(double) (i-gm_temp->size1/3)*3/gm_temp->size1;
			dr=db=0;
		}
		else
		{
			db=(double) (i-2*gm_temp->size1/3)*3/gm_temp->size1;
			dr=dg=0;
		};
		if (gm_P==NULL)
			fout<<c<<i<<c<<"[color="<<c<<dr/2+0.5<<", "<<dg/2+0.5<<", "<<db/2+0.5<<c<<", style=filled];"<<std::endl;
		else
			fout<<c<<gper[i]<<"->"<<i<<c<<"[color="<<c<<dr/2+0.5<<", "<<dg/2+0.5<<", "<<db/2+0.5<<c<<", style=filled];"<<std::endl;
		
	}
	for (int i=0;i<gm_temp->size1;i++)
	{
		for (int j=i+1;j<gm_temp->size2;j++)
			if (gsl_matrix_get(gm_temp,i,j))
			if (gm_P==NULL)
				fout<<c<<i<<c<<"--"<<c<<j<<c<<";"<<std::endl;
			else
				fout<<c<<gper[i]<<"->"<<i<<c<<"--"<<c<<gper[j]<<"->"<<j<<c<<";"<<std::endl;
			
	};
	fout<<"}"<<std::endl;
	if (gm_P!=NULL) gsl_matrix_free(gm_temp);
}

graph::~graph()
{
if (gm_A!=NULL)
 	gsl_matrix_free(gm_A);
}


//***************************************************************
//some special routines for gsl library
//***************************************************************
//matrix p-norm
double gsl_matrix_norm(const gsl_matrix* gm,double p)
{
	double res=0;
	for (int i=0;i<gm->size1;i++)
		for (int j=0;j<gm->size2;j++)
			res+=pow(abs(gsl_matrix_get(gm,i,j)),p);
	res=pow(res,(1/p));
	return res;
}
//matrix absolute value
int gsl_matrix_abs(gsl_matrix* gm)
{
	for (int i=0;i<gm->size1;i++)
		for (int j=0;j<gm->size2;j++)
			gsl_matrix_set(gm,i,j,abs(gsl_matrix_get(gm,i,j)));
	return 0;
}
//formated matrix print 
int gsl_matrix_printout(const gsl_matrix * gm,std::string sname,std::string fname_out)
{
std::ofstream fout(fname_out.c_str(),std::ios::app);
fout<<sname<<":"<<std::endl;
for (int i=0;i<gm->size1;i++)
{
	for (int j=0;j<gm->size2;j++)
		fout<<" "<<gsl_matrix_get(gm,i,j);
	fout<<std::endl;
};
return 0;
}
//formated vector print
int gsl_matrix_printout(gsl_vector * gv,std::string sname,std::string fname_out)
{
std::ofstream fout(fname_out.c_str(),std::ios::app);
fout<<sname<<":"<<std::endl;
for (int j=0;j<gv->size;j++)
	fout<<" "<<gsl_vector_get(gv,j);
fout<<std::endl;
return 0;
}
//formated vector print
int gsl_matrix_printout(gsl_permutation * gv,std::string sname,std::string fname_out)
{
std::ofstream fout(fname_out.c_str(),std::ios::app);
fout<<sname<<":"<<std::endl;
for (int j=0;j<gv->size;j++)
	fout<<" "<<gsl_permutation_get(gv,j);
fout<<std::endl;
return 0;
}

//message print
int gsl_matrix_printout(std::string sout,std::string fname_out)
{
std::ofstream fout(fname_out.c_str(),std::ios::app);
fout<<sout<<sout<<std::endl;
return 0;
}
//abs function for doubles
double abs(double x)
{return (x>0)?x:-x;}
double min(double x1,double x2)
{return (x1>x2)?x2:x1;}
double max(double x1,double x2)
{return (x1>x2)?x1:x2;}

//sum of all vector elements
double gsl_vector_sum(gsl_vector* gv)
{
	double dres=0;
	for (int j=0;j<gv->size;j++)
		dres+=gsl_vector_get(gv,j);
	return dres;
}
//sum of matrix elements (over rows or columns)
void gsl_matrix_sum(gsl_matrix* gm_A,int idim, gsl_vector* gv_res)
{
	int N;
	if (idim ==1) N=gm_A->size1; else N=gm_A->size2;
	gsl_vector* gv_ones=gsl_vector_alloc(N);
	gsl_vector_set_all(gv_ones,1);
	if (idim==1)
		gsl_blas_dgemv(CblasTrans,1,gm_A,gv_ones,0,gv_res);
	else
		gsl_blas_dgemv(CblasNoTrans,1,gm_A,gv_ones,0,gv_res);
	gsl_vector_free(gv_ones);
}
double gsl_matrix_sum(gsl_matrix* gm_A)
{
	double dres=0;
	for (int i=0;i<gm_A->size1;i++)
		for (int j=0;j<gm_A->size2;j++)
			dres+=gm_A->data[i*gm_A->size2+j];
	return dres;
}

double gsl_matrix_max_abs(gsl_matrix* A)
	{
		double dmin=abs(gsl_matrix_min(A));
		double dmax=abs(gsl_matrix_max(A));
		return (dmin>dmax)?dmin:dmax;
	};
double gsl_matrix_min(gsl_matrix* gm,double ic)//minimum greater than ic
{
	double dres=ic-1;
	for (int i=0;i<gm->size1;i++)
		for (int j=0;j<gm->size2;j++)
			if (gm->data[i*gm->size2+j]>ic)
				if (!(dres==ic-1))	dres=(dres<gm->data[i*gm->size2+j])?dres:gm->data[i*gm->size2+j];
				else dres=gm->data[i*gm->size2+j];
			
	return dres;
}
double gsl_matrix_min_abs(gsl_matrix* gm)
{
	double dres=gm->data[0];
	for (int i=0;i<gm->size1;i++)
		for (int j=0;j<gm->size2;j++)
			dres=(dres<abs(gm->data[i*gm->size2+j]))?dres:abs(gm->data[i*gm->size2+j]);
			
	return dres;
}

long gsl_numnonzero(gsl_matrix * gm_A,double deps)
{
long Nz=0;
for (int i=0;i<gm_A->size1;i++)
	for (int j=0;j<gm_A->size2;j++)
		if (abs(gm_A->data[i*gm_A->size1+j])>deps)
			Nz++;
return Nz;
}


/*!
    \fn graph::getN()
 */
long long graph::getN()
{
    return N;
}
