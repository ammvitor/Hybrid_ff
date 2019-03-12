#ifndef DSSP_PARSER_H
#define DSSP_PARSER_H

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


class DSSP_parser
{
public:
    DSSP_parser(string string);
    void reader(string string);
    vector < string > is_structured;

};

#endif // DSSP_PARSER_H
