#include "amber_parm_parser.h"

AMBER_parm_parser::AMBER_parm_parser()
{
    this->FF99SB_nonbonded("/home/jvscunha/workspace/gromacs/share/gromacs/top/amber99sb-ildn.ff/ffnonbonded.itp");
    this->FF99SB_bonded("/home/jvscunha/workspace/gromacs/share/gromacs/top/amber99sb-ildn.ff/ffbonded.itp");
    this->FF03ws_bonded("/home/jvscunha/workspace/gromacs/share/gromacs/top/amber03ws.ff/ffbonded.itp");

}

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


// parses the ff99sb nonbonded file , 
void AMBER_parm_parser::FF99SB_nonbonded(string input){

    atom_parm temporary;
    //cout << input << endl;
    ifstream myfile (input);
    string buff1,buff2,buff3,buff4,buff5,buff6,buff7,line;
    if(!myfile){
        cout << "cant open" << endl;
    }

    else{

        while (getline (myfile,line)) {

            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6 >> buff7;
        if( line[0] != ';' &&  line[0] != '[' ){
           // cout << buff1 << endl;

              temporary.atom_type_from_parm =buff1;
              temporary.atomic_number=stod(buff2);
              temporary.atomic_mass = stod(buff3);
              temporary.atomic_sigma=stod(buff6);
              temporary.atomic_epsilon=stod(buff7);
              this->atomic_parameters_FF99sb.push_back(temporary);
        }
        }



    }
}


void AMBER_parm_parser::FF99SB_bonded(string input){



    //cout << input << endl;
    ifstream myfile (input);
    string buff1,buff2,buff3,buff4,buff5,buff6,buff7,buff8,line;
    if(!myfile){
        cout << "Cant open" << endl;
    }

    else{
        int session = 0;
        bond_parm temporary_bond;
        angle_parm temporary_angle;
        proper_dihedral_parm temporary_dihedral;
        improper_dihedral_parm temporary_improper;
        while (session == 0) {


            getline (myfile,line);
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5;



            if( line[0] != ';' &&  line[0] != '[' ){
           // cout << buff1 << " " << buff2 << endl;
            temporary_bond.atom_a=buff1;
            temporary_bond.atom_b=buff2;
            temporary_bond.function=stoi(buff3);
            temporary_bond.b0=stod(buff4);
            temporary_bond.kb=stod(buff5);
            this->bondparms_vector_FF99SB.push_back(temporary_bond);
          //  cout << "bonds " << endl;

            if(line.empty()){
                session =1;
            }


        }
        }
        while (session == 1) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;
            if(line.empty()){
                session =2;
            }


        }
        }
        while (session == 2) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;
                //cout << buff1 << " " << buff2  << " " << buff3  << " " << buff4<< endl;

                temporary_angle.atom_a=buff1;
                temporary_angle.atom_b=buff2;
                temporary_angle.atom_c=buff3;
                temporary_angle.function=stoi(buff4);
                temporary_angle.th0=stod(buff5);
                temporary_angle.cth=stod(buff6);
                this->angleparms_vector_FF99SB.push_back(temporary_angle);
                //cout << "angles " << endl;

            if(line.empty()){
                session =3;
            }


        }
        }
        while (session == 3) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6 >> buff7 >> buff8;
                temporary_improper.atom_a=buff1;
                temporary_improper.atom_b=buff2;
                temporary_improper.atom_c=buff3;
                temporary_improper.atom_d=buff4;
                temporary_improper.function=stoi(buff5);
                temporary_improper.phase=stod(buff6);
                temporary_improper.kd=stod(buff7);
                temporary_improper.pn=stod(buff8);

                this->improperparms_vector_FF99SB.push_back(temporary_improper);
            //cout << "impropers " << endl;
           // cout << buff1 << " " << buff2  << " " << buff3  << " " << buff4<< endl;

            if(line.empty()){
                session =4;
            }


        }
        }
        while (session == 4) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6 >> buff7 >> buff8;


                temporary_dihedral.atom_a=buff1;
                temporary_dihedral.atom_b=buff2;
                temporary_dihedral.atom_c=buff3;
                temporary_dihedral.atom_d=buff4;
                temporary_dihedral.function=stoi(buff5);
                temporary_dihedral.phase=stod(buff6);
                temporary_dihedral.kd=stod(buff7);
                temporary_dihedral.pn=stod(buff8);
                if(buff8 == "0"){
                    this->properparms_vector_FF99SB_0.push_back(temporary_dihedral);
                    // cout << buff8 << endl;

                }
                else if(buff8 == "1"){
                    this->properparms_vector_FF99SB_1.push_back(temporary_dihedral);
                   // cout << buff8 << endl;

                }
                else if(buff8 == "2"){
                    this->properparms_vector_FF99SB_2.push_back(temporary_dihedral);
                   // cout << buff8 << endl;

                }
                else if(buff8 == "3"){
                    this->properparms_vector_FF99SB_3.push_back(temporary_dihedral);
                   // cout << buff8 << endl;
                }
               // cout << buff1 << " " << buff2  << " " << buff3  << " " <<  buff4<< endl;
               //cout <<  this->properparms_vector_FF99SB_0.size() <<  " "  <<  this->properparms_vector_FF99SB_1.size() <<  " "  <<  this->properparms_vector_FF99SB_2.size() <<  " "  <<  this->properparms_vector_FF99SB_3.size() <<  " " << endl;
            // cout << "propers " << endl;

            if(line.empty()){
                session =5;
            }


        }
        }




        temporary_dihedral.atom_a="N";
        temporary_dihedral.atom_b="CA";
        temporary_dihedral.atom_c="CB";
        temporary_dihedral.atom_d="CG2";
        temporary_dihedral.function=9;
        temporary_dihedral.phase=0.0;
        temporary_dihedral.kd=0.8158800;
        temporary_dihedral.pn=1;
        this->ILE_special_99sbparm.push_back(temporary_dihedral);
        temporary_dihedral.kd=-3.53966;
        temporary_dihedral.pn=2;
        this->ILE_special_99sbparm.push_back(temporary_dihedral);

        temporary_dihedral.atom_a="C";
        temporary_dihedral.atom_d="CG";
        temporary_dihedral.kd=2.3890640;
        temporary_dihedral.pn=1;
        this->LEU_special_99sbparm.push_back(temporary_dihedral);

        temporary_dihedral.kd=-1.4978720;
        temporary_dihedral.pn=2;
        this->LEU_special_99sbparm.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.5648400;
        temporary_dihedral.pn=3;
        this->LEU_special_99sbparm.push_back(temporary_dihedral);



        temporary_dihedral.atom_a="N";
        temporary_dihedral.kd= -11.024840;
        temporary_dihedral.pn=1;
        this->ASP_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-4.978960;
        temporary_dihedral.pn=2;
        this->ASP_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.029288;
        temporary_dihedral.pn=3;
        this->ASP_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-1.744728;
        temporary_dihedral.pn=4;
        this->ASP_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.970688;
        temporary_dihedral.pn=5;
        this->ASP_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.891192;
        temporary_dihedral.pn=6;
        this->ASP_special_99sbparm_1.push_back(temporary_dihedral);


        temporary_dihedral.atom_a="CA";
        temporary_dihedral.atom_b="CB";
        temporary_dihedral.atom_c="CG";
        temporary_dihedral.atom_d="OD";

        temporary_dihedral.kd= 0.0;
        temporary_dihedral.pn=1;
        this->ASP_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=-1.853512;
        temporary_dihedral.pn=2;
        this->ASP_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.0;
        temporary_dihedral.pn=3;
        this->ASP_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.577392;
        temporary_dihedral.pn=4;
        this->ASP_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.0;
        temporary_dihedral.pn=5;
        this->ASP_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.054392;
        temporary_dihedral.pn=6;
        this->ASP_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.atom_a="C";
        temporary_dihedral.atom_b="CA";
        temporary_dihedral.atom_c="CB";
        temporary_dihedral.atom_d="CG";

        temporary_dihedral.kd= 2.389064;
        temporary_dihedral.pn=1;
        this->ASN_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-2.493664;
        temporary_dihedral.pn=2;
        this->ASN_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.493712;
        temporary_dihedral.pn=3;
        this->ASN_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-1.744728;
        temporary_dihedral.pn=4;
        this->ASN_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.435136;
        temporary_dihedral.pn=5;
        this->ASN_special_99sbparm_1.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.422584;
        temporary_dihedral.pn=6;
        this->ASN_special_99sbparm_1.push_back(temporary_dihedral);


        temporary_dihedral.atom_a="CA";
        temporary_dihedral.atom_b="CB";
        temporary_dihedral.atom_c="CG";
        temporary_dihedral.atom_d="ND2";

        temporary_dihedral.kd= -4.376464;
        temporary_dihedral.pn=1;
        this->ASN_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.757304;
        temporary_dihedral.pn=2;
        this->ASN_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.146440;
        temporary_dihedral.pn=3;
        this->ASN_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=0.418400;
        temporary_dihedral.pn=4;
        this->ASN_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd= 0.543920;
        temporary_dihedral.pn=5;
        this->ASN_special_99sbparm_2.push_back(temporary_dihedral);

        temporary_dihedral.kd=-0.443504;
        temporary_dihedral.pn=6;
        this->ASN_special_99sbparm_2.push_back(temporary_dihedral);




        while (session == 5) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){


            if(line.empty()){
                session =6;
            }


        }
        }

    }
}

void AMBER_parm_parser::FF03ws_nonbonded(string input){
    vector < string > atom_type_from_parm;
    vector < double > atomic_number;
    vector < double > atomic_mass;
   vector < double > atomic_sigma;
   vector < double > atomic_epsilon;
    atom_parm temporary;
    vector < atom_parm > atomic_parameters;
   //cout << input << endl;
    ifstream myfile (input);
    string buff1,buff2,buff3,buff4,buff5,buff6,buff7,line;
    if(!myfile){
        cout << "Cant Open" << endl;
    }

    else{

        while (getline (myfile,line)) {

            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6 >> buff7;
        if( line[0] != ';' &&  line[0] != '[' ){
           // cout << buff3 << endl;

              temporary.atom_type_from_parm =buff1;
              temporary.atomic_number=stod(buff2);
              temporary.atomic_mass = stod(buff3);
              temporary.atomic_sigma=stod(buff6);
              temporary.atomic_epsilon=stod(buff7);
        }
        }



    }
}


void AMBER_parm_parser::FF03ws_bonded(string input){



    //cout << input << endl;
    ifstream myfile (input);
    string buff1,buff2,buff3,buff4,buff5,buff6,buff7,buff8,line;
    if(!myfile){
        cout << "cant open" << endl;
    }

    else{
        int session = 0;
        bond_parm temporary_bond;
        angle_parm temporary_angle;
        proper_dihedral_parm temporary_dihedral;
        improper_dihedral_parm temporary_improper;
        while (session == 0) {


            getline (myfile,line);
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5;



            if( line[0] != ';' &&  line[0] != '[' ){
          // cout << buff1 << " " << buff2 << endl;
            temporary_bond.atom_a=buff1;
            temporary_bond.atom_b=buff2;
            temporary_bond.function=stoi(buff3);
            temporary_bond.b0=stod(buff4);
            temporary_bond.kb=stod(buff5);
            this->bondparms_vector_FF03ws.push_back(temporary_bond);
            //cout << "bonds " << endl;
           // cout << line <<  endl;

            if(line.empty()){
                session =1;
            }


        }
        }
        while (session == 1) {
            getline (myfile,line);
            //cout << line << endl;
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;
            if(line.empty()){
                session =2;
            }


        }
        }
        while (session == 2) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;
               // cout << buff1 << " " << buff2  << " " << buff3  << " " << buff4<< endl;

                temporary_angle.atom_a=buff1;
                temporary_angle.atom_b=buff2;
                temporary_angle.atom_c=buff3;
                temporary_angle.function=stoi(buff4);
                temporary_angle.th0=stod(buff5);
                temporary_angle.cth=stod(buff6);
                this->angleparms_vector_FF03ws.push_back(temporary_angle);
               // cout << "angles " << endl;

            if(line.empty()){
                session =3;
            }


        }
        }
        temporary_improper.atom_a="N";
        temporary_improper.atom_b="CT";
        temporary_improper.atom_c="C";
        temporary_improper.atom_d="N";
        temporary_improper.function=9;
        temporary_improper.phase=285.5;
        temporary_improper.kd=0.75;
        temporary_improper.pn=1.0;
        while (session == 3) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6 >> buff7 >> buff8;
                temporary_improper.atom_a=buff1;
                temporary_improper.atom_b=buff2;
                temporary_improper.atom_c=buff3;
                temporary_improper.atom_d=buff4;
                temporary_improper.function=stoi(buff5);
                temporary_improper.phase=stod(buff6);
                temporary_improper.kd=stod(buff7);
                temporary_improper.pn=stod(buff8);

                this->improperparms_vector_FF03ws.push_back(temporary_improper);
           // cout << "impropers " << endl;
            //cout << buff1 << " " << buff2  << " " << buff3  << " " << buff4<< endl;

            if(line.empty()){
                session =4;
            }


        }
        }

        while (session == 4) {
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            if( line[0] != ';' &&  line[0] != '[' ){
                iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6 >> buff7 >> buff8;
                temporary_dihedral.atom_a=buff1;
                temporary_dihedral.atom_b=buff2;
                temporary_dihedral.atom_c=buff3;
                temporary_dihedral.atom_d=buff4;
                temporary_dihedral.function=stoi(buff5);
                temporary_dihedral.phase=stod(buff6);
                temporary_dihedral.kd=stod(buff7);
                temporary_dihedral.pn=stod(buff8);
                if(buff8 == "0"){
                    this->properparms_vector_FF03ws_0.push_back(temporary_dihedral);
                    // cout << buff8 << endl;

                }
                else if(buff8 == "1"){
                    this->properparms_vector_FF03ws_1.push_back(temporary_dihedral);
                   // cout << buff8 << endl;

                }
                else if(buff8 == "2"){
                    this->properparms_vector_FF03ws_2.push_back(temporary_dihedral);
                   // cout << buff8 << endl;

                }
                else if(buff8 == "3"){
                    this->properparms_vector_FF03ws_3.push_back(temporary_dihedral);
                   // cout << buff8 << endl;
                }
               // cout << buff1 << " " << buff2  << " " << buff3  << " " << buff4<< endl;
            // cout << "propers " << endl;

            if(line.empty()){
                session =5;
            }


        }
        }


        this->improperparms_vector_FF03ws.push_back(temporary_improper);


    }
}


