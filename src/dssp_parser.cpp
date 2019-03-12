#include "dssp_parser.h"

DSSP_parser::DSSP_parser(string dssp_input)
{
    this->reader(dssp_input);
}


void DSSP_parser::reader(string dssp_input)
{


    ifstream myfile (dssp_input);
    string buff1,line;

    if(!myfile){
        cout << "Cant open" << endl;
    }
    for(int i =0; i < 28; i++){
        getline (myfile,line);
    }
    while(getline (myfile,line)){
    buff1 = line.substr(16,1);

        if(buff1 == " "){
           // cout << "nao estruturado" << endl;
            this->is_structured.push_back("N");
        }
        else{
           // cout << "estruturado" << endl;
            this->is_structured.push_back("Y");

        }
        }
}
