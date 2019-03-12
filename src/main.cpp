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

#include <protein.h>
#include <amber_parm_parser.h>
#include <dssp_parser.h>
#include <run_engine.h>


using namespace std;


int main(int argc, char *argv[])
{
    int c;
    string input_name, output_name, dssp_file;

    while ((c = getopt(argc, argv, "i:o:d:h")) != -1)



            switch (c){
            case 'i':
                input_name = string(optarg);
                break;
            case 'd':
                dssp_file = string(optarg);
                break;
            case 'o':
                output_name = string(optarg);
                break;
            case 'h':
                printf( "Usage -i gro_file -o topology_output -d dssp_file\n" );
                exit(0);
           }



      run_engine* run = new run_engine(input_name,output_name,dssp_file);
      return 0;
}

