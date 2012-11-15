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
#ifndef RPC_H
#define RPC_H
#include <stdlib.h>
#include <string>
#include <vector>
#include <fstream>
#include <strstream>
#include <sstream>
#include <iostream>
#include <algorithm>

/**
Special root class which support parameter control, setting, reading etc.

	@author Mikhail Zaslavskiy
*/
struct parameter
{
std::string pname;
/*parameter type:
1: double
2: char
3: int
4: float
5: string
*/

int ptype;
double dvalue;
char cvalue;
int ivalue;
float fvalue;
std::string strvalue;
};
class rpc{
public:
    rpc(std::string fconfig="config.txt");
    int set_cofig(std::string fconfig) {fname_config=fconfig; return 0;};
    std::string get_config(){return fname_config;};
    int set_param(std::string pname,double pvalue);
    int set_param(std::string pname,char pvalue);
    int set_param(std::string pname,int pvalue);
    int set_param(std::string pname,float pvalue);
    int set_param(std::string pname,std::string pvalue);
    parameter get_param(std::string pname);
    double get_param_d(std::string pname){parameter p=get_param(pname); return p.dvalue;};
    int get_param_i(std::string pname){parameter p=get_param(pname); return p.ivalue;};
    char get_param_c(std::string pname){parameter p=get_param(pname); return p.cvalue;};
    std::string get_param_s(std::string pname){parameter p=get_param(pname); return p.strvalue;};
    int load_config();
    int load_config(std::string fconfig){fname_config=fconfig; return load_config();};
    int read_config(std::string sparams);
    std::string get_config_string();
    ~rpc();
    int printout(std::string fout);
private:
   std::string fname_config;
   std::vector<parameter> vparams;
   std::string sconfig;
protected:
   std::ostream* gout;
};

#endif
