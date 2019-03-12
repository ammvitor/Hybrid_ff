#ifndef PROTEIN_H
#define PROTEIN_H


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
#include <aminoacid.h>
#include <amber_parm_parser.h>
#include <dssp_parser.h>

using namespace std;

class Protein
{
public:


    struct bond_parm{
        string atom_a;
        string atom_b;
        string atom_a_index;
        string atom_b_index;
        int function;
        double b0;
        double kb;
    };
    struct angle_parm{
        string atom_a;
        string atom_b;
        string atom_c;
        string atom_a_index;
        string atom_b_index;
        string atom_c_index;
        int function;
        double th0;
        double cth;
    };
    struct proper_dihedral_parm{
        string atom_a;
        string atom_b;
        string atom_c;
        string atom_d;
        string atom_a_index;
        string atom_b_index;
        string atom_c_index;
        string atom_d_index;
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
        string atom_a_index;
        string atom_b_index;
        string atom_c_index;
        string atom_d_index;
        int function;
        double phase;
        double kd;
        double pn;
    };


    Protein(string protein_name,string dssp_file);
    void sequence(string protein_name);
    vector < vector < string > > atom_per_residue; // atomic name per residue (i is the residue and j is the name  inside of the residue)
    vector < vector < string > > index_per_residue; // atomic index per residue (i is the residue and j is the index inside of the residue)
    vector < string > residue_sequence; // vector with the 3-letter name of each residue
    vector <aminoacid*> aminoacid_parameters;
    vector< double> indexial_masses;
    void assign_aminoacids();
    void  indexer();
    void  translate_atomname_to_atomtype();
    vector < string > atomtype_serialized;
    bond_parm GetParm_bonds(string atom_a, string atom_b, AMBER_parm_parser* parmset);
    bond_parm GetParm_bonds_03ws(string atom_a, string atom_b, AMBER_parm_parser* parmset);
    angle_parm GetParm_angles(string atom_a, string atom_b, string atom_c, AMBER_parm_parser* parmset);
    angle_parm GetParm_angles_03ws(string atom_a, string atom_b, string atom_c, AMBER_parm_parser* parmset);
    proper_dihedral_parm GetParm_dihedrals(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset, int pn_index, string residuetype);
    proper_dihedral_parm GetParm_dihedrals_03ws(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset, int pn_index, string residuetype);
    improper_dihedral_parm GetParm_improper_dihedrals(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset);
    improper_dihedral_parm GetParm_improper_dihedrals_03ws(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset);
    double assign_mass(string atomtype, AMBER_parm_parser* parm);
    vector < bond_parm> bonds_parameter_list;
    vector < angle_parm> angles_parameter_list;
    vector < proper_dihedral_parm> proper_parameter_list;
    vector < improper_dihedral_parm > improper_parameter_list;
    void pair_assigner(vector < proper_dihedral_parm> proper_parameter_list);
    vector< vector <string> > pairs_ready;


    void parametrizer(string dssp_file);
    string which_HIS();

};

#endif // PROTEIN_H
