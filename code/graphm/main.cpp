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

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <stdlib.h>
#include "graph.h"
#include "experiment.h"
#include <cstdlib>
#include <iostream>
#include <stdlib.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_permutation.h>


using namespace std;
int main(int argc, char *argv[])
{
  std::string sfileconf="";
  if (argc>1)
	{ sfileconf=argv[1];
	  FILE* f=fopen(sfileconf.c_str(),"r");
	  if (f==NULL){ printf("Error: Configuration file does not exist\n"); exit(0);};
	  fclose(f);
	};
  experiment exp(sfileconf);
  std::cout<<sfileconf<<std::endl;
  if (argc>2)
	{
		std::string sconfig=(argv[2]);
		cout <<sconfig<<std::endl;
		exp.read_config(sconfig); 
	};
  exp.run_experiment();
  return EXIT_SUCCESS;
}
