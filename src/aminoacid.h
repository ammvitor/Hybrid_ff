#ifndef AMINOACID_H
#define AMINOACID_H

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
#include <dssp_parser.h>


using namespace std;


class aminoacid
{
public:
    aminoacid(string type,int C_or_N,int starting_atomic_index, string forfield_type);
    vector <int> atomic_indexes;
    int C_or_N;
    string aminoacid_type;
    vector < string> atom_names;
    vector < string> atom_types;
    vector <double> charges;
    vector < vector < string> > bonds;
    vector < vector < string> > bonds_conected;
    vector < vector < string > > impropers;
    vector < vector < string > > angles;
    vector < vector < string > > propers;
    vector < vector < string> > bonds_indexial;
    vector < vector < string > > impropers_indexial;
    vector < vector < string > > angles_indexial;
    vector < vector < string > > propers_indexial;
    vector < vector < string> > bonds_atomtype;
    vector < vector < string > > impropers_atomtype;
    vector < vector < string > > angles_atomtype;
    vector < vector < string > > propers_atomtype;
    void define_angles(vector <vector < string > > bonds);
    void define_angles_indexial(vector <vector < string > > bonds);

    void define_proper_torsions(vector <vector < string > > bonds);
    void define_proper_torsions_indexial(vector <vector < string > > bonds);

    void PARM(string residue_name);
    void SPECIAL_PARM(string residue_name);
    void PARM_03ws(string residue_name);

    vector < string> index_per_atom;
    int starting_atomic_index;
    void translator_name_to_index_bonds();
    void indexial_impropers( vector < vector < string > > impropers);
    void bond_angle_cleaner();

};

#endif // AMINOACID_H
