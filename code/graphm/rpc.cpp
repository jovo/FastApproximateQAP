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
#include "rpc.h"

rpc::rpc(std::string fconfig)
{
   if (fconfig.size()==0) fconfig="config.txt";
   fname_config=fconfig;
   load_config();
}
//parameter setting for diferent types
int rpc::set_param(std::string pname,double pvalue)
{
  bool stop_search=0;
  for (int i=0;(i<vparams.size())&(!stop_search);i++)
		if (pname.compare(vparams[i].pname)==0)
		{ vparams[i].dvalue=pvalue; stop_search=1;};
  if (!stop_search)
  {
	parameter new_param;
	new_param.pname=pname;
	new_param.ptype=1;
	new_param.dvalue=pvalue;
	vparams.push_back(new_param);
  };
  return stop_search;
}
int rpc::set_param(std::string pname,int pvalue)
{
  bool stop_search=0;
  for (int i=0;(i<vparams.size())&(!stop_search);i++)
		if (pname.compare(vparams[i].pname)==0)
		{ vparams[i].ivalue=pvalue; stop_search=1;};
  if (!stop_search)
  {
	parameter new_param;
	new_param.pname=pname;
	new_param.ptype=3;
	new_param.ivalue=pvalue;
	vparams.push_back(new_param);
  };
  return stop_search;
}
int rpc::set_param(std::string pname,char pvalue)
{
  bool stop_search=0;
  for (int i=0;(i<vparams.size())&(!stop_search);i++)
		if (pname.compare(vparams[i].pname)==0)
		{ vparams[i].cvalue=pvalue; stop_search=1;};
  if (!stop_search)
  {
	parameter new_param;
	new_param.pname=pname;
	new_param.ptype=2;
	new_param.cvalue=pvalue;
	vparams.push_back(new_param);
  };
  return stop_search;
}

int rpc::set_param(std::string pname,float pvalue)
{
  bool stop_search=0;
  for (int i=0;(i<vparams.size())&(!stop_search);i++)
		if (pname.compare(vparams[i].pname)==0)
		{ vparams[i].fvalue=pvalue; stop_search=1;};
  if (!stop_search)
  {
	parameter new_param;
	new_param.pname=pname;
	new_param.ptype=4;
	new_param.fvalue=pvalue;
	vparams.push_back(new_param);
  };
  return stop_search;
}

int rpc::set_param(std::string pname,std::string pvalue)
{
  bool stop_search=0;
  for (int i=0;(i<vparams.size())&(!stop_search);i++)
		if (pname.compare(vparams[i].pname)==0)
		{ vparams[i].strvalue=pvalue; stop_search=1;};
  if (!stop_search)
  {
	parameter new_param;
	new_param.pname=pname;
	new_param.ptype=5;
	new_param.strvalue=pvalue;
	vparams.push_back(new_param);
  };
  return stop_search;
}

//parameter getting
parameter rpc::get_param(std::string pname)
{
  bool stop_search=0;
  parameter preturn;
  for (int i=0;(i<vparams.size())&(!stop_search);i++)
		{if (pname.compare(vparams[i].pname)==0)
			{
			 preturn=vparams[i];
			 stop_search=1;
			 };}
  if (!stop_search)
  {
	preturn.pname="";
	preturn.ptype='n';
        preturn.dvalue=0;
        preturn.cvalue=' ';
	preturn.ivalue=0;
        preturn.fvalue=0;
	preturn.strvalue="";
  };
  return preturn;
}
//config file processing
int rpc::load_config()
{
	std::ifstream fconfig(fname_config.c_str());
	//assure(fconfig,fname_config.c_str());
	char buf[200];
	while (fconfig.getline(buf,200))
	{
		std::string nstr(buf);
		if ((nstr.size()>0) and (nstr.at(0)!='/')) {
		int ind1=nstr.find('=');
		int ind2=nstr.rfind(' ');
		std::string pname=nstr.substr(0,ind1);
		std::string pvalue=nstr.substr(ind1+1,ind2-ind1-1);
		std::string ptype=nstr.substr(ind2+1);
		parameter new_param;
		new_param.pname=pname;
		switch (ptype.at(0)){
		case 'd':case 'f':
			new_param.ptype=1;
			new_param.dvalue=atof(pvalue.c_str());
			break;
		case 'c':
			new_param.ptype=2;
			new_param.cvalue=pvalue.at(0);
			break;
		case 'i':
			new_param.ptype=3;
			new_param.ivalue=atoi(pvalue.c_str());
			break;
		case 's':
			new_param.ptype=5;
			new_param.strvalue=pvalue;
			break;
		};
		vparams.push_back(new_param);
		};
	};
	if (sconfig.length()>0)
		read_config(sconfig);
	return 0;
}

//config file processing
int rpc::read_config(std::string sparams)
{
        sconfig=sparams;
        std::replace( sparams.begin(), sparams.end(),';', '\n' );
	std::istrstream fconfig(sparams.c_str());
	//assure(fconfig,fname_config.c_str());
	char buf[200];
	while (fconfig.getline(buf,200))
	{
		std::string nstr(buf);
		if ((nstr.size()>0) and (nstr.at(0)!='/')) {
		int ind1=nstr.find('=');
		int ind2=nstr.length()-1;
		std::string pname=nstr.substr(0,ind1);
		std::string pvalue=nstr.substr(ind1+1,ind2-ind1-2);
		std::string ptype=nstr.substr(ind2);
		if (ptype.length()>0){
		switch (ptype.at(0)){
		case 'd':case 'f':
			set_param(pname,atof(pvalue.c_str()));
			break;
		case 'c':
			set_param(pname,pvalue.at(0));
			break;
		case 'i':
			set_param(pname,atoi(pvalue.c_str()));
			break;
		case 's':
			set_param(pname,pvalue);
			break;
		};};
		};
	};
	return 0;
}

int rpc::printout(std::string fname_out)
{
	std::ofstream fout(fname_out.c_str(),std::ios::app);
	std::vector<parameter>::iterator vp_it=vparams.begin();
	while(vp_it!=vparams.end())
	{
		fout<<vp_it->pname<<"=";
		switch (vp_it->ptype){
		case 1:
			fout<<vp_it->dvalue;
			break;
		case 2:
			fout<<vp_it->cvalue;
			break;
		case 3:
			fout<<vp_it->ivalue;
			break;
		case 4:
			fout<<vp_it->fvalue;
			break;
		case 5:
			fout<<vp_it->strvalue;
			break;
		}
		fout<<std::endl;
		vp_it++;
	};
}
std::string rpc::get_config_string()
{
	std::stringstream fout;
	std::vector<parameter>::iterator vp_it=vparams.begin();
	while(vp_it!=vparams.end())
	{
		fout<<vp_it->pname<<"=";
		switch (vp_it->ptype){
		case 1:
			fout<<vp_it->dvalue<<" d";
			break;
		case 2:
			fout<<vp_it->cvalue<<" c";
			break;
		case 3:
			fout<<vp_it->ivalue<<" i";
			break;
		case 4:
			fout<<vp_it->fvalue<<" f";
			break;
		case 5:
			fout<<vp_it->strvalue<<" s";
			break;
		}
		fout<<std::endl;
		vp_it++;
	};
	std::string sret(fout.str());
	return sret;
}
//destructor
rpc::~rpc()
{
}
/*std::string replaceAll(std::string s, std::string f, std::string r) {
  unsigned int found = s.find(f);
  while(found != std::string::npos) {
    s.replace(found, f.length(), r);
    found = s.find(f);
  }
  return s;
}
std::string replaceAll(std::string s, char f, char r) {
  unsigned int found = s.find(f);
  while(found != std::string::npos) {
    s[found]= r;
    found = s.find(f);
  }
  return s;
}
*/

