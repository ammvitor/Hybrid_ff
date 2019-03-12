#ifndef WRITER_H
#define WRITER_H

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
#include <protein.h>

using namespace std;

class writer
{
public:
    writer(string outtop,string protein_gro,string dssp_file);
    void savetopol(string outtop,string protein_gro,string dssp_file);

};

#endif // WRITER_H
