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


#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_heapsort.h>
#include "hungarian.h"

typedef struct {
  int num_rows;
  int num_cols;
  int** cost;
  int** assignment;
  int* col_inc; 
} hungarian_problem_t;
int compare_doubles(const void *av,const void *bv)
{
	double* a=(double*) av;
	double* b=(double*) bv;
	if (*a>*b) return 1;
	else if (*a<*b) return -1;
	else return 0;
}


  
#define HUNGARIAN_NOT_ASSIGNED 0 
#define HUNGARIAN_ASSIGNED 1

#define HUNGARIAN_MODE_MINIMIZE_COST   0
#define HUNGARIAN_MODE_MAXIMIZE_UTIL 1
#define INF (0x7FFFFFFF)
#define verbose (0)



int hungarian_imax(int a, int b) {
  return (a<b)?b:a;
}

int hungarian_init(hungarian_problem_t* p, int** cost_matrix, int rows, int cols, int mode) {

  int i,j, org_cols, org_rows;
  int max_cost;
  max_cost = 0;
  
  org_cols = cols;
  org_rows = rows;

  /* is the number of cols  not equal to number of rows ? 
     if yes, expand with 0-cols / 0-cols */
  rows = hungarian_imax(cols, rows);
  cols = rows;
  
  p->num_rows = rows;
  p->num_cols = cols;

  p->cost = (int**)calloc(rows,sizeof(int*));
  p->assignment = (int**)calloc(rows,sizeof(int*));
  p->col_inc = (int*)calloc(cols,sizeof(int));

  for(i=0; i<p->num_rows; i++) {
    p->cost[i] = (int*)calloc(cols,sizeof(int));
    p->assignment[i] = (int*)calloc(cols,sizeof(int));
    for(j=0; j<p->num_cols; j++) {
      p->cost[i][j] =  (i < org_rows && j < org_cols) ? cost_matrix[i][j] : 0;
      p->assignment[i][j] = 0;

      if (max_cost < p->cost[i][j])
	max_cost = p->cost[i][j];
    }
  }


  if (mode == HUNGARIAN_MODE_MAXIMIZE_UTIL) {
    for(i=0; i<p->num_rows; i++) {
      for(j=0; j<p->num_cols; j++) {
	p->cost[i][j] =  max_cost - p->cost[i][j];
      }
    }
  }
  else if (mode == HUNGARIAN_MODE_MINIMIZE_COST) {
    /* nothing to do */
  }
  else 
    //mexPrintf("%s: unknown mode. Mode was set to HUNGARIAN_MODE_MINIMIZE_COST !\n");
  
  return rows;
}




void hungarian_free(hungarian_problem_t* p) {
  int i;
  for(i=0; i<p->num_rows; i++) {
    free(p->cost[i]);
    free(p->assignment[i]);
  }
  free(p->cost);
  free(p->col_inc);
  free(p->assignment);
  p->cost = NULL;
  p->assignment = NULL;
}



void hungarian_solve(hungarian_problem_t* p)
{
  int i, j, m, n, k, l, s, t, q, unmatched, cost;
  int* col_mate;
  int* row_mate;
  int* parent_row;
  int* unchosen_row;
  int* row_dec;
  int* col_inc;
  int* slack;
  int* slack_row;

  cost=0;
  m =p->num_rows;
  n =p->num_cols;

  col_mate = (int*)calloc(p->num_rows,sizeof(int));
  unchosen_row = (int*)calloc(p->num_rows,sizeof(int));
  row_dec  = (int*)calloc(p->num_rows,sizeof(int));
  slack_row  = (int*)calloc(p->num_rows,sizeof(int));

  row_mate = (int*)calloc(p->num_cols,sizeof(int));
  parent_row = (int*)calloc(p->num_cols,sizeof(int));
  col_inc = (int*)calloc(p->num_cols,sizeof(int));
  slack = (int*)calloc(p->num_cols,sizeof(int));

  for (i=0;i<p->num_rows;i++) {
    col_mate[i]=0;
    unchosen_row[i]=0;
    row_dec[i]=0;
    slack_row[i]=0;
  }
  for (j=0;j<p->num_cols;j++) {
    row_mate[j]=0;
    parent_row[j] = 0;
    col_inc[j]=0;
    slack[j]=0;
  }

  

  /* Begin subtract column minima in order to start with lots of zeroes 12 */
  //if (verbose)    //mexPrintf("Using heuristic\n");
  for (l=0;l<n;l++)
    {
      s=p->cost[0][l];
      for (k=1;k<m;k++) 
	if (p->cost[k][l]<s)
	  s=p->cost[k][l];
      cost+=s;
      if (s!=0)
	for (k=0;k<m;k++)
	  p->cost[k][l]-=s;
    }
  /* End subtract column minima in order to start with lots of zeroes 12 */

  /* Begin initial state 16 */
  t=0;
  for (l=0;l<n;l++)
    {
      row_mate[l]= -1;
      parent_row[l]= -1;
      col_inc[l]=0;
      slack[l]=INF;
    }
 for (k=0;k<m;k++)
    {
      s=p->cost[k][0];
      for (l=1;l<n;l++)
	if (p->cost[k][l]<s)
	  s=p->cost[k][l];
      row_dec[k]=s;
      for (l=0;l<n;l++)
	if (s==p->cost[k][l] && row_mate[l]<0)
	  {
	    col_mate[k]=l;
	    row_mate[l]=k;
	    //if (verbose) mexPrintf( "matching col %d==row %d\n",l,k);
	    goto row_done;
	  }
      col_mate[k]= -1;
      //if (verbose) mexPrintf( "node %d: unmatched row %d\n",t,k);

      unchosen_row[t++]=k;
    row_done:
      ;
    }

for (i=0;i<p->num_rows;++i)
    for (j=0;j<p->num_cols;++j)
      p->assignment[i][j]=HUNGARIAN_NOT_ASSIGNED;
  /* End initial state 16 */
 
  /* Begin Hungarian algorithm 18 */
  if (t==0)
    goto done;
  unmatched=t;
  while (1)
    {
     // if (verbose) mexPrintf( "Matched %d rows.\n",m-t);
      q=0;
      while (1)
	{
	  while (q<t)
	    {
	      /* Begin explore node q of the forest 19 */
	      {
		k=unchosen_row[q];
		s=row_dec[k];
		for (l=0;l<n;l++)
		  if (slack[l])
		    {
		      int del;
		      del=p->cost[k][l]-s+col_inc[l];
		      if (del<slack[l])
			{
			  if (del==0)
			    {
			      if (row_mate[l]<0)
				goto breakthru;
			      slack[l]=0;
			      parent_row[l]=k;
			      //if (verbose) mexPrintf( "node %d: row %d==col %d--row %d\n",
				//       t,row_mate[l],l,k);
			      unchosen_row[t++]=row_mate[l];
			    }
			  else
			    {
			      slack[l]=del;
			      slack_row[l]=k;
			    }
			}
		    }
	      }
	      /* End explore node q of the forest 19 */
	      q++;
	    }
 
	  /* Begin introduce a new zero into the matrix */
	  s=INF;
	  for (l=0;l<n;l++)
	    if (slack[l] && slack[l]<s)
	      s=slack[l];
	  for (q=0;q<t;q++)
	    row_dec[unchosen_row[q]]+=s;
	  for (l=0;l<n;l++)
	    if (slack[l])
	      {
		slack[l]-=s;
		if (slack[l]==0)
		  {
		    /* Begin look at a new zero 22 */
		    k=slack_row[l];
		   // if (verbose)
		      //mexPrintf( 
			     //"Decreasing uncovered elements by %d produces zero at //[%d,%d]\n",
//			     s,k,l);
		    if (row_mate[l]<0)
		      {
			for (j=l+1;j<n;j++)
			  if (slack[j]==0)
			    col_inc[j]+=s;
			goto breakthru;
		      }
		    else
		      {
			parent_row[l]=k;
			//if (verbose)
			  //mexPrintf("node %d: row %d==col %d--row %d\n",t,row_mate[l],l,k);
			unchosen_row[t++]=row_mate[l];
		      }
		    /* End look at a new zero 22 */
		  }
	      }
	    else
	      col_inc[l]+=s;
	  /* End introduce a new zero into the matrix 21 */
	}
    breakthru:
      /* Begin update the matching 20 */
      //if (verbose)
	//mexPrintf( "Breakthrough at node %d of %d!\n",q,t);
      while (1)
	{
	  j=col_mate[k];
	  col_mate[k]=l;
	  row_mate[l]=k;
	  //if (verbose)
	    //mexPrintf( "rematching col %d==row %d\n",l,k);
	  if (j<0)
	    break;
	  k=parent_row[j];
	  l=j;
	}
      /* End update the matching 20 */
      if (--unmatched==0)
	goto done;
      /* Begin get ready for another stage 17 */
      t=0;
      for (l=0;l<n;l++)
	{
	  parent_row[l]= -1;
	  slack[l]=INF;
	}
      for (k=0;k<m;k++)
	if (col_mate[k]<0)
	  {
	    //if (verbose)
	      //mexPrintf( "node %d: unmatched row %d\n",t,k);
	    unchosen_row[t++]=k;
	  }
      /* End get ready for another stage 17 */
    }
 done:

  /* Begin doublecheck the solution 23 */
  for (k=0;k<m;k++)
    for (l=0;l<n;l++)
      if (p->cost[k][l]<row_dec[k]-col_inc[l])
	exit(0);
  for (k=0;k<m;k++)
    {
      l=col_mate[k];
      if (l<0 || p->cost[k][l]!=row_dec[k]-col_inc[l])
	exit(0);
    }
  k=0;
  for (l=0;l<n;l++)
    if (col_inc[l])
      k++;
  if (k>m)
    exit(0);
  /* End doublecheck the solution 23 */
  /* End Hungarian algorithm 18 */
  for (l=0;l<n;++l)
      p->col_inc[l]=col_inc[l];
  for (i=0;i<m;++i)
    {
      p->assignment[i][col_mate[i]]=HUNGARIAN_ASSIGNED;
      /*TRACE("%d - %d\n", i, col_mate[i]);*/
    }
  for (k=0;k<m;++k)
    {
      for (l=0;l<n;++l)
	{
	  /*TRACE("%d ",p->cost[k][l]-row_dec[k]+col_inc[l]);*/
	  p->cost[k][l]=p->cost[k][l]-row_dec[k]+col_inc[l];
	}
      /*TRACE("\n");*/
    }
  for (i=0;i<m;i++)
    cost+=row_dec[i];
  for (i=0;i<n;i++)
    cost-=col_inc[i];

  free(slack);
  free(col_inc);
  free(parent_row);
  free(row_mate);
  free(slack_row);
  free(row_dec);
  free(unchosen_row);
  free(col_mate);
}



gsl_matrix* gsl_matrix_hungarian(gsl_matrix* gm_C)
{
	gsl_matrix* gm_P=gsl_matrix_alloc(gm_C->size1,gm_C->size2);
	gsl_matrix_hungarian(gm_C,gm_P);
	return gm_P;
}
void gsl_matrix_hungarian(gsl_matrix* gm_C,gsl_matrix* gm_P,gsl_vector* gv_col_inc, gsl_permutation* gp_sol, int _bprev_init, gsl_matrix *gm_C_denied, bool bgreedy)
{
  
  long dim, startdim, enddim, n1,n2;
  double *C;
  int i,j;
  int **m;
  double *z;
  hungarian_problem_t p, *q;
  int matrix_size;
  double C_min=gsl_matrix_min(gm_C)-1;
  n1 = gm_C->size1;    /* first dimension of the cost matrix */
  n2 = gm_C->size2;    /* second dimension of the cost matrix */
  C = gm_C->data; 
   //greedy solution
   if (bgreedy)
   {
	int ind,ind1,ind2;
	size_t *C_ind=new size_t[n1*n2];
	gsl_heapsort_index(C_ind,C,n1*n2,sizeof(double),compare_doubles);
        bool* bperm_fix_1=new bool[n1]; bool* bperm_fix_2=new bool[n2]; int inummatch=0;
	for (i=0;i<n1;i++) {bperm_fix_1[i]=false;bperm_fix_2[i]=false;};
	gsl_matrix_set_zero(gm_P);
	for (long l=0;l<n1*n2;l++)
	{
		ind=C_ind[l];
		ind1=floor(ind/n1);
		ind2=ind%n2;
		
		if (!bperm_fix_1[ind1] and !bperm_fix_2[ind2])
		{
			bperm_fix_1[ind1]=true; bperm_fix_2[ind2]=true;
			gm_P->data[ind]=1;inummatch++;
		};
		if (inummatch==n1) break;
	};
	delete[] bperm_fix_1;delete[] bperm_fix_2;
	//because C is a transpose matrix
	gsl_matrix_transpose(gm_P);
	return;	
   };
  double C_max=((gsl_matrix_max(gm_C)-C_min>1)?(gsl_matrix_max(gm_C)-C_min):1)*(n1>n2?n1:n2);
  m = (int**)calloc(n1,sizeof(int*)); 
  for (i=0;i<n1;i++)
        {
        	m[i] = (int*)calloc(n2,sizeof(int));  
        	for (j=0;j<n2;j++)
            		m[i][j] = (int) (C[i+n1*j] - C_min);
		if (gm_C_denied!=NULL)
		for (j=0;j<n2;j++){
			if (j==30)
				int dbg=1;
			bool bden=(gm_C_denied->data[n2*i+j]<1e-10);
            		if (bden) m[i][j] =C_max;
			else 
				int dbg=1;
			};
 	};
    //normalization: rows and columns
    double dmin;
    for (i=0;i<n1;i++)
        {
        	dmin=m[i][0];
        	for (j=1;j<n2;j++)
            		dmin= (m[i][j]<dmin)? m[i][j]:dmin;
        	for (j=0;j<n2;j++)
            		m[i][j]-=dmin;
 	};
    for (j=0;j<n2;j++)
        {
        	dmin=m[0][j];
        	for (i=1;i<n1;i++)
            		dmin= (m[i][j]<dmin)? m[i][j]:dmin;
        	for (i=0;i<n1;i++)
            		m[i][j]-=dmin;
 	};
   if ((_bprev_init) &&(gv_col_inc !=NULL))
	{
	//dual solution v substraction
		for (j=0;j<n2;j++)
        		for (i=0;i<n1;i++)
				m[i][j]-=gv_col_inc->data[j];
	//permutation of m columns
		int *mt = new int[n2];
		for (i=0;i<n1;i++)
		{
			for (j=0;j<n2;j++) mt[j]=m[i][j];
			for (j=0;j<n2;j++) m[i][j]=mt[gsl_permutation_get(gp_sol,j)];
		};
		delete[] mt;
		
	};

   
  /* initialize the hungarian_problem using the cost matrix*/
   matrix_size = hungarian_init(&p, m , n1,n2, HUNGARIAN_MODE_MINIMIZE_COST) ;
  /* solve the assignement problem */
  hungarian_solve(&p);
  q = &p;
  //gsl_matrix* gm_P=gsl_matrix_alloc(n1,n2);
  gsl_permutation* gp_sol_inv=gsl_permutation_alloc(n2);
  if (gp_sol!=NULL)
  	gsl_permutation_inverse(gp_sol_inv,gp_sol);
  else
	gsl_permutation_init(gp_sol_inv);
  for (i=0;i<n1;i++)
         for (j=0;j<n2;j++)
              gsl_matrix_set(gm_P,i,j,q->assignment[i][gp_sol_inv->data[j]]);
  //initialization by the previous solution
  if ((_bprev_init) &&(gv_col_inc !=NULL))
        for (j=0;j<n2;j++)
		gv_col_inc->data[j]=q->col_inc[gp_sol_inv->data[j]];
  if ((_bprev_init) && (gp_sol!=NULL))
  {
  for (i=0;i<n1;i++)
         for (j=0;j<n2;j++)
  		if (gsl_matrix_get(gm_P,i,j)==HUNGARIAN_ASSIGNED)
			gp_sol->data[i]=j;
  };
  /* free used memory */
  gsl_permutation_free(gp_sol_inv);
  hungarian_free(&p);
  for (i=0;i<n1;i++)
        free(m[i]);
  free(m);
  //return gm_P;
}



