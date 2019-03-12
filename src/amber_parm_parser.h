#ifndef AMBER_PARM_PARSER_H
#define AMBER_PARM_PARSER_H

#include <iostream>
#include<vector>
#include<iostream>
#include<fstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<string>
#include<cstring>
#include<time.h>
#include <iostream>
#include <iomanip>
#include <algorithm>
#include <stdio.h>
#include <string>
#include <getopt.h>
#include <stdio.h>
#include <sstream>


using namespace std;


class AMBER_parm_parser
{
public:

    AMBER_parm_parser();
    struct atom_parm{
        string atom_type_from_parm;
        int function;
        double atomic_number;
        double atomic_mass;
        double atomic_sigma;
        double atomic_epsilon;
    };
    
    struct bond_parm{
        string atom_a;
        string atom_b;
        int function;
        double b0;
        double kb;
    };
    struct angle_parm{
        string atom_a;
        string atom_b;
        string atom_c;
        int function;
        double th0;
        double cth;
    };
    struct proper_dihedral_parm{
        string atom_a;
        string atom_b;
        string atom_c;
        string atom_d;
        int function;
        double phase;
        double kd;
        double pn;
    
    };
    struct improper_dihedral_parm{
        string atom_a;
        string atom_b;
        string atom_c;
        string atom_d;
        int function;
        double phase;
        double kd;
        double pn;
    
    };
    

    vector < atom_parm > atomic_parameters_FF99sb;
    vector < atom_parm > atomic_parameters_FF03ws;

    void FF99SB_nonbonded(string input);
    void FF99SB_bonded(string input);
    void FF03ws_nonbonded(string input);
    void FF03ws_bonded(string input);
    vector < bond_parm > bondparms_vector_FF99SB;
    vector < angle_parm > angleparms_vector_FF99SB;
    vector < proper_dihedral_parm > properparms_vector_FF99SB_0;
    vector < proper_dihedral_parm > properparms_vector_FF99SB_1;
    vector < proper_dihedral_parm > properparms_vector_FF99SB_2;
    vector < proper_dihedral_parm > properparms_vector_FF99SB_3;
    vector < proper_dihedral_parm > properparms_vector_FF03ws_0;
    vector < proper_dihedral_parm > properparms_vector_FF03ws_1;
    vector < proper_dihedral_parm > properparms_vector_FF03ws_2;
    vector < proper_dihedral_parm > properparms_vector_FF03ws_3;
    vector < improper_dihedral_parm > improperparms_vector_FF99SB;
    vector < bond_parm > bondparms_vector_FF03ws;
    vector < angle_parm > angleparms_vector_FF03ws;
    vector < proper_dihedral_parm > properparms_vector_FF03ws;
    vector < improper_dihedral_parm > improperparms_vector_FF03ws;


  // vector < double > atomic_mass
    vector< proper_dihedral_parm > ILE_special_99sbparm;
    vector< proper_dihedral_parm > LEU_special_99sbparm;
    vector< proper_dihedral_parm > ASP_special_99sbparm_1;
    vector< proper_dihedral_parm > ASP_special_99sbparm_2;
    vector< proper_dihedral_parm > ASN_special_99sbparm_1;
    vector< proper_dihedral_parm > ASN_special_99sbparm_2;
;
  //  vector < double > atomic_sigma;
  //  vector < double > atomic_epsilon;

};

#endif // AMBER_PARM_PARSER_H
