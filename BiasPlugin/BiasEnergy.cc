/*******************************************************************************/
/* This module implements a structure-based (Go-like) energy term to be used   */
/* with the PROFASI protein simulation software package                        */
/*******************************************************************************/
#include "BiasEnergy.hh"
//#include <Aux/Constants.hh>
//#include <Aux/RMSD_Utils.hh>
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <sstream>

namespace prf
{
  BiasEnergy::BiasEnergy() : Energy()
  {
    Name("BiasEnergy");
  }
  
  void BiasEnergy::set_pars(std::string s)
  {
    filename=s;
  }     
  
  void BiasEnergy::init()
  {
    std::string infile1="",infile2="";
    size_t lst1[1000],lst2[1000];
    size_t lst1_2[1000],lst2_2[1000];
    size_t res1,res2;
    double dist;
    double a1,a2;

    /* initialize */
    
    for (size_t i=0; i < MAXCONT; ++i) lst1[i] = lst2[i] = lst1_2[i] = lst2_2[i] = -1;
    for (size_t i=0; i < MAXCONT; ++i) atm1[i] = atm2[i] = atm_dist[i] = -1;
    for (size_t i=0; i < MAXCONT; ++i) doublet[i] = -1;
    lambda_SC = 0.0;
    ksi_SC = 1.0;

    /* root */

    prf_xml::node_ptr_type root = prf_xml::get_xml_tree(filename);
    if (root==nullptr or root->name()!="parameters") {
      if (root==nullptr) {
	prf::cout<<Name()<<"> No valid XML tree in "<<filename<<"\n";
      } else {
	prf::cout<<Name()<<"> Root node in "<<filename <<" is not <parameter>\n";
      }
    }

    /* contact maps */

    prf_xml::node_ptr_type par=root->child("contact_maps");
    if (par->name()!="contact_maps"){
      prf::cerr<<"No contact map filenames found in "<<filename<<"\n";
    } 
    if (par->child("cmap_sc_1")!=nullptr) {
      infile1=par->child("cmap_sc_1")->value();
    } 
    if (par->child("cmap_sc_2")!=nullptr) {
      infile2=par->child("cmap_sc_2")->value();
    } 

    /* contact strengths */

    par=root->child("sc_param");
    if (par->name()!="sc_param"){
      prf::cerr<<"No parameters found in "<<filename<<"\n";
    } 
    if (par->child("lambda_SC")!=nullptr) {
      lambda_SC=strtod(par->child("lambda_SC")->value().c_str(),nullptr);
    } 
    if (par->child("ksi_SC")!=nullptr) {
      ksi_SC=strtod(par->child("ksi_SC")->value().c_str(),nullptr);
    }
    
    prf::cout <<"<"<<filename<<"> Reading contact map 1 from file " << infile1 << "\n";
    prf::cout <<"<"<<filename<<"> Reading contact map 2 from file " << infile2 << "\n";
    prf::cout <<"<"<<filename<<"> Setting lambda_SC " << lambda_SC << "\n";
    prf::cout <<"<"<<filename<<"> Setting ksi_SC " << ksi_SC << "\n";

    std::ifstream inFile(infile1);
    std::ifstream inFile2(infile2);

    if (!inFile)  prf::cout << Name() << "> No file " << infile1 << "\n";
    if (!inFile2) prf::cout << Name() << "> No file " << infile2 << "\n";

    ncont = 0;
    while(inFile >> res1 >> res2 >> dist) {
      lst1[ncont] = res1;
      lst2[ncont] = res2;
      
      atm1[ncont] = ( p->amino_acid(res1)->OLC() == G ?
		      p->amino_acid(res1)->Calpha().UniqueId() :
		      p->amino_acid(res1)->sidechain_atom(0).UniqueId() );
      atm2[ncont] = ( p->amino_acid(res2)->OLC() == G ?
		      p->amino_acid(res2)->Calpha().UniqueId() :
		      p->amino_acid(res2)->sidechain_atom(0).UniqueId() );
      atm_dist[ncont] = dist;

      lamSC1[ncont] = lambda_SC;	

      ++ncont;
    }
    
    if (ncont > 0) prf::cout << Name() << "> Number of native contacts read 1: " << ncont << "\n";
    if (ncont >= MAXCONT) {prf::cerr << Name() << "Too many native contacts 1 \n"; exit(-1);}
    
    if (ncont > 0) {
      std::ofstream outFile("./contactmap.sidechain.out");
      for (size_t i=0; i < ncont; ++i) {
	outFile << lst1[i] << " " << atm1[i] << " " << p->amino_acid(lst1[i])->TLC() << " "
		<< lst2[i] << " " << atm2[i] << " " << p->amino_acid(lst2[i])->TLC() << " "
		<< atm_dist[i] << "  " << lamSC1[i] << "\n";
      }
      outFile.close();
    }
    
    ncont2 = 0;
    while(inFile2 >> res1 >> res2 >> dist) {
      lst1_2[ncont2] = res1;
      lst2_2[ncont2] = res2;
      
      atm1_2[ncont2] = ( p->amino_acid(res1)->OLC() == G ?
			 p->amino_acid(res1)->Calpha().UniqueId() :
			 p->amino_acid(res1)->sidechain_atom(0).UniqueId() );
      atm2_2[ncont2] = ( p->amino_acid(res2)->OLC() == G ?
			p->amino_acid(res2)->Calpha().UniqueId() :
			p->amino_acid(res2)->sidechain_atom(0).UniqueId() );
      atm_dist2[ncont2] = dist;

      lamSC2[ncont2] = lambda_SC;

      ++ncont2;
    }

    if (ncont2 > 0) prf::cout << Name() << "> Number of native contacts read 2: " << ncont2 << "\n";
    if (ncont2 >= MAXCONT) {prf::cerr << Name() << "Too many native contacts 2\n"; exit(-1);}
    
    if (ncont2 > 0) {
      std::ofstream outFile2("./contactmap2.sidechain.out");
      for (size_t i=0; i < ncont2; ++i) {
	outFile2 << lst1_2[i] << " " << atm1_2[i] << " " << p->amino_acid(lst1_2[i])->TLC() << " "
		 << lst2_2[i] << " " << atm2_2[i] << " " << p->amino_acid(lst2_2[i])->TLC() << " "
		 << atm_dist2[i] << "  " << lamSC2[i] << "\n";
      }
      outFile2.close();
    }

    /* common contacts in contact maps */
    
    for (size_t i=0; i < ncont; ++i) {
      int j;

      a1 = atm1[i];
      a2 = atm2[i];

      for (size_t j=0; j < ncont2; ++j) {
	if ( (a1 == atm1_2[j] && a2 == atm2_2[j]) || (a2 == atm1_2[j] && a1 == atm2_2[j]) ) doublet[i] = j;
      }

      if ((j = doublet[i]) >= 0) prf::cout << Name() << "> Doublet: 1 "
		  			   << lst1[i]   << " " << lst2[i]   << " " << atm_dist[i]  << " 2 "
					   << lst1_2[j] << " " << lst2_2[j] << " " << atm_dist2[j] << "\n";
    }
    
    inFile.close();
    inFile2.close();

    return ;
  }
  
  double BiasEnergy::evaluate() // calculate total bias energy
  {
    size_t a1,a2;
    int j;
    double r,dr;
    double v1[1000],v2[1000];
    
    for (size_t i = 0; i < ncont ; ++i) {
      a1 = atm1[i];
      a2 = atm2[i];
      r = (p->atom(a1).Pos() - p->atom(a2).Pos()).mag();
      dr = r - atm_dist[i];
      v1[i] = - lamSC1[i] * exp(- dr * dr / 2 / ksi_SC);
    }

    for (size_t i = 0; i < ncont2 ; ++i) {
      a1 = atm1_2[i];
      a2 = atm2_2[i];
      r = (p->atom(a1).Pos() - p->atom(a2).Pos()).mag();
      dr = r - atm_dist2[i];
      v2[i] = - lamSC2[i] * exp(- dr * dr / 2 / ksi_SC);
    }

    vval = 0;
    for (size_t i = 0; i < ncont ; ++i) {
      if ((j = doublet[i]) >=0 ) {
	v1[i] = (v1[i] < v2[j] ? v1[i] : v2[j]);
	v2[j] = 0;
      }
      vval += v1[i];
    }

    for (size_t i = 0; i < ncont2 ; ++i) vval += v2[i];
    
    return vval;
  }
}


