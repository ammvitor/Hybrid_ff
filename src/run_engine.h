#ifndef RUN_ENGINE_H
#define RUN_ENGINE_H
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
#include <dssp_parser.h>
#include <amber_parm_parser.h>
#include <protein.h>
#include <writer.h>

using namespace std;


class run_engine
{
public:
    run_engine(string protein_gro_file, string outputfile, string dssp_file);
};

#endif // RUN_ENGINE_H
