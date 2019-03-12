#include "protein.h"


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
    int function;
    double phase;
    double kd;
    double pn;

};


Protein::Protein(string protein_name,string dssp_file)
{
    this->sequence(protein_name);
    this->assign_aminoacids();
    this->indexer();
    cout << "Parametrizing the protein" << endl;
    this->parametrizer(dssp_file);

}

void Protein::sequence(string protein_name){
// Reads the GRO file and creates a vector of strings representing the aminoacids
    cout << protein_name << endl;
    // variable declaration
    ifstream myfile (protein_name);
    string buff1,buff2,buff3,buff4,buff5,buff6,line,residue_index,residue_index_minus_one;
    vector <string> atoms_in_the_residue;
    vector <string> index_per_atom_in_the_residue;
    int number_of_atoms;
    if(!myfile){
        cout << "cant open" << endl;
    }

    //parsing the gro file
    if(myfile){
        int breaker = 0;
        getline (myfile,line);
        getline (myfile,line);
        std::istringstream iss(line.c_str());
        iss >> buff1 ;
        //reads the second line of the GRO file that has the number of atoms.
        number_of_atoms = stoi(buff1.c_str());
        for(int i =0; i < number_of_atoms; i++ ){
            getline (myfile,line);
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;

            residue_index = buff1.substr(0, buff1.size()-3);

            // this reads each residue int he sequence by comparing to the previous one
            if(i ==0){
               this->residue_sequence.push_back(buff1.substr(buff1.size()-3));
               residue_index = buff1.substr(0, buff1.size()-3);
               //assign the previous one as the first one
               residue_index_minus_one = residue_index;
           }
           else if(residue_index != residue_index_minus_one){
                this->atom_per_residue.push_back(atoms_in_the_residue);
                atoms_in_the_residue.clear();
                this->index_per_residue.push_back(index_per_atom_in_the_residue);
                index_per_atom_in_the_residue.clear();
                cout << buff1.substr(buff1.size()-3) <<endl ;



                /*if(buff1.substr(buff1.size()-3) == "HIS"){
                    for(int j =0; j < 17; j++ ){
                     cout << buff2 << endl;
                    if(buff2 == "HD1"){
                        cout << "HD11 LOL" << endl;
                        for(int k =0; k < 9; k++ ){
                            getline (myfile,line);
                            std::istringstream iss(line.c_str());
                            iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;
                            if(buff2 == "HE2"){
                                cout << "HD12 LOL" << endl;

                                this->residue_sequence.push_back("HIP");
                                breaker=1;

                                break;
                            }
                            else if(buff2 == "N"){
                                cout << "NNNNN" << endl;
                                this->residue_sequence.push_back("HID");
                                breaker=1;

                                break;
                            }
                        }

                    }
                    else if(buff2 == "HE2"){
                        cout << "HE2 LOL" << endl;

                        this->residue_sequence.push_back("HIE");
                        breaker=1;
                               break;
                    }
                    if(breaker == 1){
                        break;
                    }

                    getline (myfile,line);
                    std::istringstream iss(line.c_str());
                    iss >> buff1 >> buff2 >> buff3 >> buff4 >> buff5 >> buff6;
                    //this->residue_sequence.push_back(buff1.substr(buff1.size()-3));
                   }
                }*/
                //else{
                    this->residue_sequence.push_back(buff1.substr(buff1.size()-3));
                    cout << "else" << endl;
                    cout << buff1.substr(buff1.size()-3) << endl;

                //}

           }
           index_per_atom_in_the_residue.push_back(buff3);

           atoms_in_the_residue.push_back(buff2);

           residue_index_minus_one = residue_index;

           if(number_of_atoms-1 == i){
               //reading the last atom.
               this->atom_per_residue.push_back(atoms_in_the_residue);
               atoms_in_the_residue.clear();
               this->index_per_residue.push_back(index_per_atom_in_the_residue);
               index_per_atom_in_the_residue.clear();
           }

        }




    }



}
string Protein::which_HIS(){

}

void Protein::assign_aminoacids(){
    int atomic_index = 1;
    /*
    using the assigned sequence in Protein::sequence - it reads the amonioacids.rtp and assingn the
    molecular structural parameters - bonds and repsective dihedrals
    */
    for(int i =0; i < this->residue_sequence.size(); i++ ){
        if(i == 0){
       // creates the n terminal aminoacid object using the assigned name
            aminoacid* aminoacid_test = new aminoacid(this->residue_sequence[i],1,atomic_index,"99sb");
            this->aminoacid_parameters.push_back(aminoacid_test);
            atomic_index = atomic_index+aminoacid_test->atom_names.size();
            cout << "N-TERMINAL - Atoms in this residue " << aminoacid_test->atom_names.size() << endl;
            }
            //creates for the c terminal residue in the sequence
            else if(i == this->residue_sequence.size()-1){
                aminoacid* aminoacid_test = new aminoacid(this->residue_sequence[i],3,atomic_index,"99sb");
                this->aminoacid_parameters.push_back(aminoacid_test);
                cout << "C-TERMINAL - Atoms in this residue " << aminoacid_test->atom_names.size() << endl;
                atomic_index = atomic_index+aminoacid_test->atom_names.size();

            }
            else{
                //creates and appends the  middle aminozcids objects
                aminoacid* aminoacid_test = new aminoacid(this->residue_sequence[i],2,atomic_index,"99sb");
                this->aminoacid_parameters.push_back(aminoacid_test);
                cout << "Atoms in this residue " << aminoacid_test->atom_names.size() << endl;
                atomic_index = atomic_index+aminoacid_test->atom_names.size();

            }
        }

    }


void  Protein::indexer(){
    //reindex bonds and previous bonds


    /*
    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        for(int k = 0; k < this->aminoacid_parameters[i]->atom_names.size() ; k++){
            cout <<  this->aminoacid_parameters[i]->atomic_indexes[k] << " " << this->aminoacid_parameters[i]->atom_names[k]  << " " << this->atom_per_residue[i][k] << endl;
        }
    }

    */


    /*
     *
     *
     *
     *
     */

    for(int i =1; i < this->aminoacid_parameters.size(); i++){
        for(int j = 0; j < this->aminoacid_parameters[i]->bonds.size() ; j++){
            for(int n =0; n < 2 ; n++){
                for(int k = 0; k < this->aminoacid_parameters[i-1]->atom_names.size() ; k++){



               if(this->aminoacid_parameters[i]->bonds[j][n] == "-C" && this->aminoacid_parameters[i-1]->atom_names[k] == "C" ){
                   this->aminoacid_parameters[i]->bonds_conected[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->bonds_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                  // cout << this->aminoacid_parameters[i-1]->atom_names[k] << " " << this->aminoacid_parameters[i]->bonds[j][n] << endl;
                   break;

               }
               else if(this->aminoacid_parameters[i]->bonds[j][n] == "-N" && this->aminoacid_parameters[i-1]->atom_names[k] == "N" ){
                   this->aminoacid_parameters[i]->bonds_conected[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->bonds_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                  // cout << this->aminoacid_parameters[i-1]->atom_names[k] << " " << this->aminoacid_parameters[i]->bonds[j][n] << endl;
                   break;



               }
               else if(this->aminoacid_parameters[i]->bonds[j][n] == "-CA" && this->aminoacid_parameters[i-1]->atom_names[k] == "CA" ){
                   this->aminoacid_parameters[i]->bonds_conected[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->bonds_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                  // cout << this->aminoacid_parameters[i-1]->atom_names[k] << " " << this->aminoacid_parameters[i]->bonds[j][n] << endl;
                   break;



               }
               else if(this->aminoacid_parameters[i]->bonds[j][n] == "-HA" && this->aminoacid_parameters[i-1]->atom_names[k] == "HA" ){
                   this->aminoacid_parameters[i]->bonds_conected[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->bonds_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
//                   cout << this->aminoacid_parameters[i-1]->atom_names[k] << " " << this->aminoacid_parameters[i]->bonds[j][n] << endl;
                   break;



               }
               else if(this->aminoacid_parameters[i]->bonds[j][n] == "-O" && this->aminoacid_parameters[i-1]->atom_names[k] == "O" ){
                   this->aminoacid_parameters[i]->bonds_conected[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->bonds_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   //cout << this->aminoacid_parameters[i-1]->atom_names[k] << " " << this->aminoacid_parameters[i]->bonds[j][n] << endl;
                   break;



               }
               else if(this->aminoacid_parameters[i]->bonds[j][n] == "-H" && this->aminoacid_parameters[i-1]->atom_names[k] == "H" ){
                   this->aminoacid_parameters[i]->bonds_conected[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->bonds_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   //cout << this->aminoacid_parameters[i-1]->atom_names[k] << " " << this->aminoacid_parameters[i]->bonds[j][n] << endl;
                   break;



               }

               }


           }

        }

    }
    /*
    for(int i =0; i < this->aminoacid_parameters.size() ; i++ ){
        for(int j =0; j < this->aminoacid_parameters[i]->angles.size() ; j++ ){
            cout << this->aminoacid_parameters[i]->angles[j][0] << " "
                 << this->aminoacid_parameters[i]->angles[j][1] << " "
                 << this->aminoacid_parameters[i]->angles[j][2] << endl;
        }

    }*/
    cout << "Indexing angles!" << endl;



    for(int i =1; i < this->aminoacid_parameters.size(); i++){
        for(int j = 0; j < this->aminoacid_parameters[i]->angles.size() ; j++){
            for(int n =0; n < 3 ; n++){
                for(int k = 0; k < this->aminoacid_parameters[i-1]->atom_names.size() ; k++){



               if(this->aminoacid_parameters[i]->angles[j][n] == "-C" && this->aminoacid_parameters[i-1]->atom_names[k] == "C" ){

                   this->aminoacid_parameters[i]->angles[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->angles_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);


               }
               else if(this->aminoacid_parameters[i]->angles[j][n] == "-N" && this->aminoacid_parameters[i-1]->atom_names[k] == "N" ){

                   this->aminoacid_parameters[i]->angles[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->angles_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);


               }
               else if(this->aminoacid_parameters[i]->angles[j][n] == "-CA" && this->aminoacid_parameters[i-1]->atom_names[k] == "CA" ){

                   this->aminoacid_parameters[i]->angles[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->angles_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);


               }
               else if((this->aminoacid_parameters[i]->angles[j][n] == "-HA" && this->aminoacid_parameters[i-1]->atom_names[k] == "HA" )){

                   this->aminoacid_parameters[i]->angles[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->angles_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);


               }
               else if(this->aminoacid_parameters[i]->angles[j][n] == "-O" && this->aminoacid_parameters[i-1]->atom_names[k] == "O" ){

                   this->aminoacid_parameters[i]->angles[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->angles_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);


               }
               else if(this->aminoacid_parameters[i]->angles[j][n] == "-H" && this->aminoacid_parameters[i-1]->atom_names[k] == "H" ){

                   this->aminoacid_parameters[i]->angles[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->angles_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);


               }
               }


           }

        }
    }
    /*
    for(int i =0; i < this->aminoacid_parameters.size() ; i++ ){
        for(int j =0; j < this->aminoacid_parameters[i]->angles.size() ; j++ ){
            cout << this->aminoacid_parameters[i]->angles_indexial[j][0] << " "
                 << this->aminoacid_parameters[i]->angles_indexial[j][1] << " "
                 << this->aminoacid_parameters[i]->angles_indexial[j][2] << endl;
        }

    }
    */
    cout << "Indexing Propers and improper dihedrals!"  << endl;

    for(int i =1; i < this->aminoacid_parameters.size(); i++){
        for(int j = 0; j < this->aminoacid_parameters[i]->propers.size() ; j++){
            for(int n =0; n < 4 ; n++){


                for(int k = 0; k < this->aminoacid_parameters[i-1]->atom_names.size() ; k++){


               if(this->aminoacid_parameters[i]->propers[j][n] == "-C" && this->aminoacid_parameters[i-1]->atom_names[k] == "C" ){
                   this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   break;


               }
               else if(this->aminoacid_parameters[i]->propers[j][n] == "-N" && this->aminoacid_parameters[i-1]->atom_names[k] == "N" ){
                   this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   break;


               }
               else if(this->aminoacid_parameters[i]->propers[j][n] == "-CA" && this->aminoacid_parameters[i-1]->atom_names[k] == "CA" ){
                   this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   break;


               }
               else if((this->aminoacid_parameters[i]->propers[j][n] == "-HA" && this->aminoacid_parameters[i-1]->atom_names[k] == "HA") || (this->aminoacid_parameters[i]->propers[j][n] == "-HA" && this->aminoacid_parameters[i-1]->atom_names[k] == "HA1") ){
                   this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   break;


               }
               else if(this->aminoacid_parameters[i]->propers[j][n] == "-O" && this->aminoacid_parameters[i-1]->atom_names[k] == "O" ){
                   this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   break;


               }
               else if(this->aminoacid_parameters[i]->propers[j][n] == "-H" && this->aminoacid_parameters[i-1]->atom_names[k] == "H" ){
                   this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                   this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                   break;


               }

               else if(this->aminoacid_parameters[i]->propers[j][n] == "-R"){



                   if(this->aminoacid_parameters[i-1]->aminoacid_type == "GLY" ||
                      this->aminoacid_parameters[i-1]->aminoacid_type == "CGLY" ||
                      this->aminoacid_parameters[i-1]->aminoacid_type == "NGLY"){

                       for(int p =0; p < this->aminoacid_parameters[i-1]->atom_names.size();p++){
                            if( this->aminoacid_parameters[i-1]->atom_names[p] == "HA2"){
                                this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[p] ;
                                this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[p]);

                            }
                       }
                   }
                   else{
                       for(int p =0; p < this->aminoacid_parameters[i-1]->atom_names.size();p++){
                        if( this->aminoacid_parameters[i-1]->atom_names[p] == "CB"){
                            this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i-1]->atom_names[p] ;
                            this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[p]);

                        }
                    }

                   }




               }
               }


           }

        }
    }
    for(int i =0; i < this->aminoacid_parameters.size()-1; i++){
            for(int j = 0; j < this->aminoacid_parameters[i]->propers.size() ; j++){
                //cout << this->aminoacid_parameters[i]->aminoacid_type << " " << this->aminoacid_parameters[i]->propers.size() << endl;

                for(int n =0; n < 4 ; n++){
                    for(int k = 0; k < this->aminoacid_parameters[i+1]->atom_names.size() ; k++){



                     if(this->aminoacid_parameters[i]->propers[j][n] == "+N" && this->aminoacid_parameters[i+1]->atom_names[k] == "N" ){
                       this->aminoacid_parameters[i]->propers[j][n] = this->aminoacid_parameters[i+1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->propers_indexial[j][n] = to_string(this->aminoacid_parameters[i+1]->atomic_indexes[k]);
                       break;


                     }
                   }


               }

            }
        }



    for(int i =1; i < this->aminoacid_parameters.size(); i++){
            for(int j = 0; j < this->aminoacid_parameters[i]->impropers.size() ; j++){

                for(int n =0; n < 4 ; n++){
                    for(int k = 0; k < this->aminoacid_parameters[i-1]->atom_names.size() ; k++){


                   if(this->aminoacid_parameters[i]->impropers[j][n] == "-C" && this->aminoacid_parameters[i-1]->atom_names[k] == "C" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                       break;


                   }
                   else if(this->aminoacid_parameters[i]->impropers[j][n] == "-N" && this->aminoacid_parameters[i-1]->atom_names[k] == "N" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                       break;


                   }
                   else if(this->aminoacid_parameters[i]->impropers[j][n] == "-CA" && this->aminoacid_parameters[i-1]->atom_names[k] == "CA" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                       break;


                   }
                   else if((this->aminoacid_parameters[i]->impropers[j][n] == "-HA" && this->aminoacid_parameters[i-1]->atom_names[k] == "HA") || (this->aminoacid_parameters[i]->impropers[j][n] == "-HA" && this->aminoacid_parameters[i-1]->atom_names[k] == "HA1") ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                       break;


                   }
                   else if(this->aminoacid_parameters[i]->impropers[j][n] == "-O" && this->aminoacid_parameters[i-1]->atom_names[k] == "O" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                       break;


                   }
                   else if(this->aminoacid_parameters[i]->impropers[j][n] == "-H" && this->aminoacid_parameters[i-1]->atom_names[k] == "H" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[k]);
                       break;


                   }
                   else if(this->aminoacid_parameters[i]->impropers[j][n] == "-R"){
                       if(this->aminoacid_parameters[i-1]->aminoacid_type != "GLY" || this->aminoacid_parameters[i-1]->aminoacid_type != "CGLY" || this->aminoacid_parameters[i-1]->aminoacid_type != "NGLY"){
                       for(int p =0; p < this->aminoacid_parameters[i-1]->atom_names.size();p++){
                       if( this->aminoacid_parameters[i-1]->atom_names[p] == "CB"){
                            this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[p] ;
                            this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[p]);

                       }
                       }
                       }
                       else{
                           for(int p =0; p < this->aminoacid_parameters[i-1]->atom_names.size();p++){

                            if( this->aminoacid_parameters[i-1]->atom_names[p] == "HA2"){
                                this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i-1]->atom_names[p] ;
                                this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i-1]->atomic_indexes[p]);

                            }
                           }
                       }

                   }

                   else if(this->aminoacid_parameters[i]->impropers[j][n] == "+N" && this->aminoacid_parameters[i+1]->atom_names[k] == "N" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i+1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i+1]->atomic_indexes[k]);
                       break;


                   }
                   }


               }

            }
        }

    for(int i =0; i < this->aminoacid_parameters.size()-1; i++){
            for(int j = 0; j < this->aminoacid_parameters[i]->impropers.size() ; j++){
                //cout << this->aminoacid_parameters[i]->aminoacid_type << " " << this->aminoacid_parameters[i]->impropers.size() << endl;

                for(int n =0; n < 4 ; n++){
                    for(int k = 0; k < this->aminoacid_parameters[i+1]->atom_names.size() ; k++){



                     if(this->aminoacid_parameters[i]->impropers[j][n] == "+N" && this->aminoacid_parameters[i+1]->atom_names[k] == "N" ){
                       this->aminoacid_parameters[i]->impropers[j][n] = this->aminoacid_parameters[i+1]->atom_names[k] ;
                       this->aminoacid_parameters[i]->impropers_indexial[j][n] = to_string(this->aminoacid_parameters[i+1]->atomic_indexes[k]);
                       break;


                     }
                   }


               }

            }
        }



    /*
    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        for(int j = 0; j < this->aminoacid_parameters[i]->propers.size() ; j++){
            cout << "already connected " << this->aminoacid_parameters[i]->aminoacid_type <<  endl;
            cout << this->aminoacid_parameters[i]->propers[j][0] << " "
                 << this->aminoacid_parameters[i]->propers[j][1] << " "
                 << this->aminoacid_parameters[i]->propers[j][2] << " "
                 << this->aminoacid_parameters[i]->propers[j][3] << " " << endl;
            cout << this->aminoacid_parameters[i]->propers_indexial[j][0] << " "
                << this->aminoacid_parameters[i]->propers_indexial[j][1] << " "
                << this->aminoacid_parameters[i]->propers_indexial[j][2] << " "
                 << this->aminoacid_parameters[i]->propers_indexial[j][3] << " " << endl;

        }
    }*/

    this->translate_atomname_to_atomtype();

}

void  Protein::translate_atomname_to_atomtype(){

 for(int i =0; i < this->aminoacid_parameters.size(); i++){
     for(int j = 0; j < this->aminoacid_parameters[i]->atom_types.size() ; j++){
        this->atomtype_serialized.push_back(this->aminoacid_parameters[i]->atom_types[j]);
     }
 }

}

Protein::bond_parm Protein::GetParm_bonds(string atom_a, string atom_b, AMBER_parm_parser* parmset){


    bond_parm bond_ready;
    for(int i =0; i < parmset->bondparms_vector_FF99SB.size(); i++){

        if(parmset->bondparms_vector_FF99SB[i].atom_a == atom_a && parmset->bondparms_vector_FF99SB[i].atom_b == atom_b){
           //cout <<  parmset->bondparms_vector_FF99SB[i].atom_a << " " <<  atom_a  << " " << parmset->bondparms_vector_FF99SB[i].atom_b  << " " << atom_b << endl;


            bond_ready.atom_a=parmset->bondparms_vector_FF99SB[i].atom_a;
            bond_ready.atom_b=parmset->bondparms_vector_FF99SB[i].atom_b;
            bond_ready.function=parmset->bondparms_vector_FF99SB[i].function;
            bond_ready.b0=parmset->bondparms_vector_FF99SB[i].b0;
            bond_ready.kb=parmset->bondparms_vector_FF99SB[i].kb;
            return bond_ready;

        }
        else if(parmset->bondparms_vector_FF99SB[i].atom_a == atom_b && parmset->bondparms_vector_FF99SB[i].atom_b == atom_a){

          // cout <<  parmset->bondparms_vector_FF99SB[i].atom_a << " " <<  atom_b  << " " << parmset->bondparms_vector_FF99SB[i].atom_b  << " " << atom_a << endl;
           bond_ready.atom_a=parmset->bondparms_vector_FF99SB[i].atom_b;
           bond_ready.atom_b=parmset->bondparms_vector_FF99SB[i].atom_a;
           bond_ready.function=parmset->bondparms_vector_FF99SB[i].function;
           bond_ready.b0=parmset->bondparms_vector_FF99SB[i].b0;
           bond_ready.kb=parmset->bondparms_vector_FF99SB[i].kb;
           return bond_ready;
        }

    }
}




Protein::bond_parm Protein::GetParm_bonds_03ws(string atom_a, string atom_b, AMBER_parm_parser* parmset){


    bond_parm bond_ready;
    for(int i =0; i < parmset->bondparms_vector_FF03ws.size(); i++){

        if(parmset->bondparms_vector_FF03ws[i].atom_a == atom_a && parmset->bondparms_vector_FF03ws[i].atom_b == atom_b){
           //cout <<  parmset->bondparms_vector_FF03ws[i].atom_a << " " <<  atom_a  << " " << parmset->bondparms_vector_FF03ws[i].atom_b  << " " << atom_b << endl;


            bond_ready.atom_a=parmset->bondparms_vector_FF03ws[i].atom_a;
            bond_ready.atom_b=parmset->bondparms_vector_FF03ws[i].atom_b;
            bond_ready.function=parmset->bondparms_vector_FF03ws[i].function;
            bond_ready.b0=parmset->bondparms_vector_FF03ws[i].b0;
            bond_ready.kb=parmset->bondparms_vector_FF03ws[i].kb;
            return bond_ready;

        }
        else if(parmset->bondparms_vector_FF03ws[i].atom_a == atom_b && parmset->bondparms_vector_FF03ws[i].atom_b == atom_a){

          // cout <<  parmset->bondparms_vector_FF03ws[i].atom_a << " " <<  atom_b  << " " << parmset->bondparms_vector_FF03ws[i].atom_b  << " " << atom_a << endl;
           bond_ready.atom_a=parmset->bondparms_vector_FF03ws[i].atom_b;
           bond_ready.atom_b=parmset->bondparms_vector_FF03ws[i].atom_a;
           bond_ready.function=parmset->bondparms_vector_FF03ws[i].function;
           bond_ready.b0=parmset->bondparms_vector_FF03ws[i].b0;
           bond_ready.kb=parmset->bondparms_vector_FF03ws[i].kb;
           return bond_ready;
        }

    }
}



Protein::angle_parm Protein::GetParm_angles(string atom_a, string atom_b, string atom_c, AMBER_parm_parser* parmset){


    angle_parm angle_ready;
   // cout << "angle " <<  atom_a << " " << atom_b << " " << atom_c <<  endl;

    for(int i =0; i < parmset->angleparms_vector_FF99SB.size(); i++){

        if(parmset->angleparms_vector_FF99SB[i].atom_a == atom_a &&
           parmset->angleparms_vector_FF99SB[i].atom_b == atom_b &&
           parmset->angleparms_vector_FF99SB[i].atom_c == atom_c){

            angle_ready.atom_a=parmset->angleparms_vector_FF99SB[i].atom_a;
            angle_ready.atom_b=parmset->angleparms_vector_FF99SB[i].atom_b;
            angle_ready.atom_c=parmset->angleparms_vector_FF99SB[i].atom_c;
            angle_ready.function=parmset->angleparms_vector_FF99SB[i].function;
            angle_ready.th0=parmset->angleparms_vector_FF99SB[i].th0;
            angle_ready.cth=parmset->angleparms_vector_FF99SB[i].cth;

            return angle_ready;
        }
        else if(parmset->angleparms_vector_FF99SB[i].atom_c == atom_a &&
                parmset->angleparms_vector_FF99SB[i].atom_b == atom_b &&
                parmset->angleparms_vector_FF99SB[i].atom_a == atom_c){

            angle_ready.atom_a=parmset->angleparms_vector_FF99SB[i].atom_c;
            angle_ready.atom_b=parmset->angleparms_vector_FF99SB[i].atom_b;
            angle_ready.atom_c=parmset->angleparms_vector_FF99SB[i].atom_a;
            angle_ready.function=parmset->angleparms_vector_FF99SB[i].function;
            angle_ready.th0=parmset->angleparms_vector_FF99SB[i].th0;
            angle_ready.cth=parmset->angleparms_vector_FF99SB[i].cth;

            return angle_ready;
        }
    }
    cout << "Missing angles parameters for :"  << atom_a << " " << atom_b << " " << atom_c <<  endl;
    exit(0);


}



Protein::angle_parm Protein::GetParm_angles_03ws(string atom_a, string atom_b, string atom_c, AMBER_parm_parser* parmset){


    angle_parm angle_ready;
   // cout << "angle " <<  atom_a << " " << atom_b << " " << atom_c <<  endl;

    for(int i =0; i < parmset->angleparms_vector_FF03ws.size(); i++){

        if(parmset->angleparms_vector_FF03ws[i].atom_a == atom_a &&
           parmset->angleparms_vector_FF03ws[i].atom_b == atom_b &&
           parmset->angleparms_vector_FF03ws[i].atom_c == atom_c){

            angle_ready.atom_a=parmset->angleparms_vector_FF03ws[i].atom_a;
            angle_ready.atom_b=parmset->angleparms_vector_FF03ws[i].atom_b;
            angle_ready.atom_c=parmset->angleparms_vector_FF03ws[i].atom_c;
            angle_ready.function=parmset->angleparms_vector_FF03ws[i].function;
            angle_ready.th0=parmset->angleparms_vector_FF03ws[i].th0;
            angle_ready.cth=parmset->angleparms_vector_FF03ws[i].cth;

            return angle_ready;
        }
        else if(parmset->angleparms_vector_FF03ws[i].atom_c == atom_a &&
                parmset->angleparms_vector_FF03ws[i].atom_b == atom_b &&
                parmset->angleparms_vector_FF03ws[i].atom_a == atom_c){

            angle_ready.atom_a=parmset->angleparms_vector_FF03ws[i].atom_c;
            angle_ready.atom_b=parmset->angleparms_vector_FF03ws[i].atom_b;
            angle_ready.atom_c=parmset->angleparms_vector_FF03ws[i].atom_a;
            angle_ready.function=parmset->angleparms_vector_FF03ws[i].function;
            angle_ready.th0=parmset->angleparms_vector_FF03ws[i].th0;
            angle_ready.cth=parmset->angleparms_vector_FF03ws[i].cth;

            return angle_ready;
        }
    }
    cout << "Missing angles parameters for :"  << atom_a << " " << atom_b << " " << atom_c <<  endl;
    exit(0);


}

Protein::proper_dihedral_parm Protein::GetParm_dihedrals(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset, int pn_index, string residuetype){


    proper_dihedral_parm proper_ready;
    //cout << residuetype << " " <<   atom_a << " " << atom_b << " " << atom_c << " " << atom_d << endl;

    if( residuetype == "ILE" &&
         atom_a == "N" && atom_b == "CT" &&  atom_c == "CT" && atom_d == "CT"){
         if(pn_index == 1){
            proper_ready.atom_a = atom_a;
            proper_ready.atom_b = atom_b;
            proper_ready.atom_c = atom_c;
            proper_ready.atom_d = atom_d;
            proper_ready.pn=parmset->ILE_special_99sbparm[0].pn;
            proper_ready.kd=parmset->ILE_special_99sbparm[0].kd;
            proper_ready.function=parmset->ILE_special_99sbparm[0].function;
            proper_ready.phase=parmset->ILE_special_99sbparm[0].phase;
            return(proper_ready);
         }
         else if (pn_index == 2){
            proper_ready.atom_a = atom_a;
            proper_ready.atom_b = atom_b;
            proper_ready.atom_c = atom_c;
            proper_ready.atom_d = atom_d;
            proper_ready.pn=parmset->ILE_special_99sbparm[1].pn;
            proper_ready.kd=parmset->ILE_special_99sbparm[1].kd;
            proper_ready.function=parmset->ILE_special_99sbparm[1].function;
            proper_ready.phase=parmset->ILE_special_99sbparm[1].phase;
            return(proper_ready);
         }
    }
    else if( residuetype == "ILE" &&
         atom_d == "N" &&  atom_c == "CT" && atom_b == "CT" && atom_a == "CT"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ILE_special_99sbparm[0].pn;
           proper_ready.kd=parmset->ILE_special_99sbparm[0].kd;
           proper_ready.function=parmset->ILE_special_99sbparm[0].function;
           proper_ready.phase=parmset->ILE_special_99sbparm[0].phase;
           return(proper_ready);
        }
        else if (pn_index = 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ILE_special_99sbparm[1].pn;
           proper_ready.kd=parmset->ILE_special_99sbparm[1].kd;
           proper_ready.function=parmset->ILE_special_99sbparm[1].function;
           proper_ready.phase=parmset->ILE_special_99sbparm[1].phase;
           return(proper_ready);
        }

    }

    else if (residuetype == "LEU"&&
             atom_a == "N" && atom_b == "CT" &&  atom_c == "CT" && atom_d == "CT"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->LEU_special_99sbparm[0].pn;
           proper_ready.kd=parmset->LEU_special_99sbparm[0].kd;
           proper_ready.function=parmset->LEU_special_99sbparm[0].function;
           proper_ready.phase=parmset->LEU_special_99sbparm[0].phase;
           return(proper_ready);
        }
        else if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->LEU_special_99sbparm[1].pn;
           proper_ready.kd=parmset->LEU_special_99sbparm[1].kd;
           proper_ready.function=parmset->LEU_special_99sbparm[1].function;
           proper_ready.phase=parmset->LEU_special_99sbparm[1].phase;
           return(proper_ready);
        }
        else if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->LEU_special_99sbparm[2].pn;
           proper_ready.kd=parmset->LEU_special_99sbparm[2].kd;
           proper_ready.function=parmset->LEU_special_99sbparm[2].function;
           proper_ready.phase=parmset->LEU_special_99sbparm[2].phase;
           return(proper_ready);
        }
    }

    else if (residuetype == "LEU"&&
             atom_d == "N" &&  atom_c == "CT" && atom_b == "CT" && atom_a == "CT"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->LEU_special_99sbparm[0].pn;
           proper_ready.kd=parmset->LEU_special_99sbparm[0].kd;
           proper_ready.function=parmset->LEU_special_99sbparm[0].function;
           proper_ready.phase=parmset->LEU_special_99sbparm[0].phase;
           return(proper_ready);
        }
        else if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->LEU_special_99sbparm[1].pn;
           proper_ready.kd=parmset->LEU_special_99sbparm[1].kd;
           proper_ready.function=parmset->LEU_special_99sbparm[1].function;
           proper_ready.phase=parmset->LEU_special_99sbparm[1].phase;
           return(proper_ready);
        }
        else if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->LEU_special_99sbparm[2].pn;
           proper_ready.kd=parmset->LEU_special_99sbparm[2].kd;
           proper_ready.function=parmset->LEU_special_99sbparm[2].function;
           proper_ready.phase=parmset->LEU_special_99sbparm[2].phase;
           return(proper_ready);
        }
    }


    else if (residuetype == "ASP"&&
             atom_a == "N" && atom_b == "CT" &&  atom_c == "CT" && atom_d == "C"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[0].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[0].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[0].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[1].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[1].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[1].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[2].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[2].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[2].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[3].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[3].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[3].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[4].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[4].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[4].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[5].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[5].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[5].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[5].phase;
           return(proper_ready);
        }
    }
    else if (residuetype == "ASP"&&
             atom_d == "N" &&  atom_c == "CT" && atom_b == "CT" && atom_a == "C"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[0].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[0].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[0].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[1].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[1].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[1].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[2].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[2].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[2].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[3].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[3].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[3].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[4].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[4].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[4].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_1[5].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_1[5].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_1[5].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_1[5].phase;
           return(proper_ready);
        }
    }




    else if (residuetype == "ASP"&&
             atom_a == "CT" && atom_b == "CT" &&  atom_c == "C" && atom_d == "02"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[0].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[0].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[0].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[1].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[1].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[1].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[2].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[2].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[2].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[3].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[3].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[3].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[4].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[4].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[4].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[5].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[5].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[5].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[5].phase;
           return(proper_ready);
        }
    }

    else if (residuetype == "ASP"&&
             atom_d == "CT" &&  atom_c == "CT" && atom_b == "C" && atom_a == "O2"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[0].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[0].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[0].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[1].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[1].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[1].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[2].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[2].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[2].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[3].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[3].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[3].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[4].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[4].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[4].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASP_special_99sbparm_2[5].pn;
           proper_ready.kd=parmset->ASP_special_99sbparm_2[5].kd;
           proper_ready.function=parmset->ASP_special_99sbparm_2[5].function;
           proper_ready.phase=parmset->ASP_special_99sbparm_2[5].phase;
           return(proper_ready);
        }

    }






    else if (residuetype == "ASN"&&
             atom_a == "C" && atom_b == "CT" &&  atom_c == "CT" && atom_d == "C"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[0].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[0].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[0].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[1].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[1].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[1].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[2].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[2].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[2].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[3].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[3].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[3].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[4].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[4].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[4].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[5].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[5].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[5].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[5].phase;
           return(proper_ready);
        }

    }
    else if (residuetype == "ASN"&&
             atom_d == "C" &&  atom_c == "CT" && atom_b == "CT" && atom_a == "C"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[0].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[0].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[0].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[1].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[1].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[1].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[2].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[2].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[2].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[3].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[3].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[3].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[4].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[4].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[4].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_1[5].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_1[5].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_1[5].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_1[5].phase;
           return(proper_ready);
        }
    }




    else if (residuetype == "ASN"&&
             atom_a == "CT" && atom_b == "CT" &&  atom_c == "C" && atom_d == "N"){
        if(pn_index == 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[0].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[0].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[0].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[0].phase;
           return(proper_ready);
        }
        if(pn_index == 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[1].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[1].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[1].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[1].phase;
           return(proper_ready);
        }
        if(pn_index == 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[2].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[2].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[2].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[2].phase;
           return(proper_ready);
        }
        if(pn_index == 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[3].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[3].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[3].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[3].phase;
           return(proper_ready);
        }
        if(pn_index == 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[4].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[4].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[4].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[4].phase;
           return(proper_ready);
        }
        if(pn_index == 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[5].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[5].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[5].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[5].phase;
           return(proper_ready);
        }

    }
    else if (residuetype == "ASN"&&
             atom_d == "CT" &&  atom_c == "CT" && atom_b == "CT" && atom_a == "N"){
        if(pn_index = 1){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[0].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[0].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[0].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[0].phase;
           return(proper_ready);
        }
        if(pn_index = 2){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[1].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[1].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[1].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[1].phase;
           return(proper_ready);
        }
        if(pn_index = 3){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[2].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[2].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[2].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[2].phase;
           return(proper_ready);
        }
        if(pn_index = 4){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[3].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[3].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[3].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[3].phase;
           return(proper_ready);
        }
        if(pn_index = 5){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[4].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[4].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[4].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[4].phase;
           return(proper_ready);
        }
        if(pn_index = 6){
           proper_ready.atom_a = atom_a;
           proper_ready.atom_b = atom_b;
           proper_ready.atom_c = atom_c;
           proper_ready.atom_d = atom_d;
           proper_ready.pn=parmset->ASN_special_99sbparm_2[5].pn;
           proper_ready.kd=parmset->ASN_special_99sbparm_2[5].kd;
           proper_ready.function=parmset->ASN_special_99sbparm_2[5].function;
           proper_ready.phase=parmset->ASN_special_99sbparm_2[5].phase;
           return(proper_ready);
        }

    }

    else{
    if(pn_index == 0){
    for(int i =0; i < parmset->properparms_vector_FF99SB_0.size(); i++){
      // cout  << parmset->properparms_vector_FF99SB_0[i].atom_a << " " << parmset->properparms_vector_FF99SB_0[i].atom_b << " "
       //      << parmset->properparms_vector_FF99SB_0[i].atom_c << " " << parmset->properparms_vector_FF99SB_0[i].atom_d << endl;
       // cout  << "atoms input  "<< atom_a << " " << atom_b << " " << atom_c << " " << atom_d << endl;


        if(parmset->properparms_vector_FF99SB_0[i].atom_a == atom_a &&
           parmset->properparms_vector_FF99SB_0[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_0[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_0[i].atom_d == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_0[i].atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_0[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_0[i].atom_c;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_0[i].atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_0[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_0[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_0[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_0[i].pn;

        }

        else if(parmset->properparms_vector_FF99SB_0[i].atom_d == atom_a &&
           parmset->properparms_vector_FF99SB_0[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_0[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_0[i].atom_a == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_0[i].atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_0[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_0[i].atom_b;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_0[i].atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_0[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_0[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_0[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_0[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_0[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_0[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_0[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_0[i].atom_a == "X"){

            proper_ready.atom_a=atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_0[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_0[i].atom_b;
            proper_ready.atom_d=atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_0[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_0[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_0[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_0[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_0[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_0[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_0[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_0[i].atom_a == "X"){


            proper_ready.atom_a=atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_0[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_0[i].atom_c;
            proper_ready.atom_d=atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_0[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_0[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_0[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_0[i].pn;

        }
    }

   }
    else if(pn_index == 1){
    for(int i =0; i < parmset->properparms_vector_FF99SB_1.size(); i++){
       // cout  << parmset->properparms_vector_FF99SB_1[i].atom_a << " " << parmset->properparms_vector_FF99SB_1[i].atom_b << " "
       //       << parmset->properparms_vector_FF99SB_1[i].atom_c << " " << parmset->properparms_vector_FF99SB_1[i].atom_d << endl;



        if(parmset->properparms_vector_FF99SB_1[i].atom_a == atom_a &&
           parmset->properparms_vector_FF99SB_1[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_1[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_1[i].atom_d == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_1[i].atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_1[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_1[i].atom_c;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_1[i].atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_1[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_1[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_1[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_1[i].pn;

        }

        else if(parmset->properparms_vector_FF99SB_1[i].atom_d == atom_a &&
           parmset->properparms_vector_FF99SB_1[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_1[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_1[i].atom_a == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_1[i].atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_1[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_1[i].atom_b;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_1[i].atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_1[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_1[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_1[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_1[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_1[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_1[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_1[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_1[i].atom_a == "X"){

            proper_ready.atom_a=atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_1[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_1[i].atom_b;
            proper_ready.atom_d=atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_1[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_1[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_1[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_1[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_1[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_1[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_1[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_1[i].atom_a == "X"){


            proper_ready.atom_a=atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_1[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_1[i].atom_c;
            proper_ready.atom_d=atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_1[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_1[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_1[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_1[i].pn;

        }
    }
    return proper_ready;

   }
    else if(pn_index == 2){
    for(int i =0; i < parmset->properparms_vector_FF99SB_2.size(); i++){
       // cout  << parmset->properparms_vector_FF99SB_2[i].atom_a << " " << parmset->properparms_vector_FF99SB_2[i].atom_b << " "
       //       << parmset->properparms_vector_FF99SB_2[i].atom_c << " " << parmset->properparms_vector_FF99SB_2[i].atom_d << endl;



        if(parmset->properparms_vector_FF99SB_2[i].atom_a == atom_a &&
           parmset->properparms_vector_FF99SB_2[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_2[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_2[i].atom_d == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_2[i].atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_2[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_2[i].atom_c;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_2[i].atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_2[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_2[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_2[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_2[i].pn;

        }

        else if(parmset->properparms_vector_FF99SB_2[i].atom_d == atom_a &&
           parmset->properparms_vector_FF99SB_2[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_2[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_2[i].atom_a == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_2[i].atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_2[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_2[i].atom_b;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_2[i].atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_2[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_2[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_2[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_2[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_2[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_2[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_2[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_2[i].atom_a == "X"){

            proper_ready.atom_a=atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_2[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_2[i].atom_b;
            proper_ready.atom_d=atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_2[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_2[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_2[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_2[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_2[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_2[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_2[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_2[i].atom_a == "X"){


            proper_ready.atom_a=atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_2[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_2[i].atom_c;
            proper_ready.atom_d=atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_2[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_2[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_2[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_2[i].pn;

        }
    }
    return proper_ready;
   }
    else if(pn_index == 3){
    for(int i =0; i < parmset->properparms_vector_FF99SB_3.size(); i++){
       // cout  << parmset->properparms_vector_FF99SB_3[i].atom_a << " " << parmset->properparms_vector_FF99SB_3[i].atom_b << " "
       //       << parmset->properparms_vector_FF99SB_3[i].atom_c << " " << parmset->properparms_vector_FF99SB_3[i].atom_d << endl;



        if(parmset->properparms_vector_FF99SB_3[i].atom_a == atom_a &&
           parmset->properparms_vector_FF99SB_3[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_3[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_3[i].atom_d == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_3[i].atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_3[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_3[i].atom_c;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_3[i].atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_3[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_3[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_3[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_3[i].pn;

        }

        else if(parmset->properparms_vector_FF99SB_3[i].atom_d == atom_a &&
           parmset->properparms_vector_FF99SB_3[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_3[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_3[i].atom_a == atom_d){

            proper_ready.atom_a=parmset->properparms_vector_FF99SB_3[i].atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_3[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_3[i].atom_b;
            proper_ready.atom_d=parmset->properparms_vector_FF99SB_3[i].atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_3[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_3[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_3[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_3[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_3[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_3[i].atom_c == atom_b &&
           parmset->properparms_vector_FF99SB_3[i].atom_b == atom_c &&
           parmset->properparms_vector_FF99SB_3[i].atom_a == "X"){

            proper_ready.atom_a=atom_d;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_3[i].atom_c;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_3[i].atom_b;
            proper_ready.atom_d=atom_a;

            proper_ready.function=parmset->properparms_vector_FF99SB_3[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_3[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_3[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_3[i].pn;
        }
        else if(parmset->properparms_vector_FF99SB_3[i].atom_d == "X" &&
           parmset->properparms_vector_FF99SB_3[i].atom_c == atom_c &&
           parmset->properparms_vector_FF99SB_3[i].atom_b == atom_b &&
           parmset->properparms_vector_FF99SB_3[i].atom_a == "X"){


            proper_ready.atom_a=atom_a;
            proper_ready.atom_b=parmset->properparms_vector_FF99SB_3[i].atom_b;
            proper_ready.atom_c=parmset->properparms_vector_FF99SB_3[i].atom_c;
            proper_ready.atom_d=atom_d;

            proper_ready.function=parmset->properparms_vector_FF99SB_3[i].function;
            proper_ready.phase=parmset->properparms_vector_FF99SB_3[i].phase;
            proper_ready.kd=parmset->properparms_vector_FF99SB_3[i].kd;
            proper_ready.pn=parmset->properparms_vector_FF99SB_3[i].pn;

        }

   }
    return proper_ready;

}
    }


}

Protein::proper_dihedral_parm Protein::GetParm_dihedrals_03ws(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset, int pn_index, string residuetype){


    proper_dihedral_parm proper_ready;
    //cout << residuetype << " " <<   atom_a << " " << atom_b << " " << atom_c << " " << atom_d << endl;

   if(pn_index == 0){
     for(int i =0; i < parmset->properparms_vector_FF03ws_0.size(); i++){
       // cout  << parmset->properparms_vector_FF03ws_0[i].atom_a << " " << parmset->properparms_vector_FF03ws_0[i].atom_b << " "
        //      << parmset->properparms_vector_FF03ws_0[i].atom_c << " " << parmset->properparms_vector_FF03ws_0[i].atom_d << endl;
        // cout  << "atoms input  "<< atom_a << " " << atom_b << " " << atom_c << " " << atom_d << endl;


         if(parmset->properparms_vector_FF03ws_0[i].atom_a == atom_a &&
            parmset->properparms_vector_FF03ws_0[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_0[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_0[i].atom_d == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_0[i].atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_0[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_0[i].atom_c;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_0[i].atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_0[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_0[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_0[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_0[i].pn;

         }

         else if(parmset->properparms_vector_FF03ws_0[i].atom_d == atom_a &&
            parmset->properparms_vector_FF03ws_0[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_0[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_0[i].atom_a == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_0[i].atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_0[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_0[i].atom_b;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_0[i].atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_0[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_0[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_0[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_0[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_0[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_0[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_0[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_0[i].atom_a == "X"){

             proper_ready.atom_a=atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_0[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_0[i].atom_b;
             proper_ready.atom_d=atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_0[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_0[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_0[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_0[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_0[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_0[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_0[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_0[i].atom_a == "X"){


             proper_ready.atom_a=atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_0[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_0[i].atom_c;
             proper_ready.atom_d=atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_0[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_0[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_0[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_0[i].pn;

         }
     }

    }
     else if(pn_index == 1){
     for(int i =0; i < parmset->properparms_vector_FF03ws_1.size(); i++){
        // cout  << parmset->properparms_vector_FF03ws_1[i].atom_a << " " << parmset->properparms_vector_FF03ws_1[i].atom_b << " "
        //       << parmset->properparms_vector_FF03ws_1[i].atom_c << " " << parmset->properparms_vector_FF03ws_1[i].atom_d << endl;



         if(parmset->properparms_vector_FF03ws_1[i].atom_a == atom_a &&
            parmset->properparms_vector_FF03ws_1[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_1[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_1[i].atom_d == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_1[i].atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_1[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_1[i].atom_c;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_1[i].atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_1[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_1[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_1[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_1[i].pn;

         }

         else if(parmset->properparms_vector_FF03ws_1[i].atom_d == atom_a &&
            parmset->properparms_vector_FF03ws_1[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_1[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_1[i].atom_a == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_1[i].atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_1[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_1[i].atom_b;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_1[i].atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_1[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_1[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_1[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_1[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_1[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_1[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_1[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_1[i].atom_a == "X"){

             proper_ready.atom_a=atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_1[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_1[i].atom_b;
             proper_ready.atom_d=atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_1[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_1[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_1[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_1[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_1[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_1[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_1[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_1[i].atom_a == "X"){


             proper_ready.atom_a=atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_1[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_1[i].atom_c;
             proper_ready.atom_d=atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_1[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_1[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_1[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_1[i].pn;

         }
     }
     return proper_ready;

    }
     else if(pn_index == 2){
     for(int i =0; i < parmset->properparms_vector_FF03ws_2.size(); i++){
        // cout  << parmset->properparms_vector_FF03ws_2[i].atom_a << " " << parmset->properparms_vector_FF03ws_2[i].atom_b << " "
        //       << parmset->properparms_vector_FF03ws_2[i].atom_c << " " << parmset->properparms_vector_FF03ws_2[i].atom_d << endl;



         if(parmset->properparms_vector_FF03ws_2[i].atom_a == atom_a &&
            parmset->properparms_vector_FF03ws_2[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_2[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_2[i].atom_d == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_2[i].atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_2[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_2[i].atom_c;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_2[i].atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_2[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_2[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_2[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_2[i].pn;

         }

         else if(parmset->properparms_vector_FF03ws_2[i].atom_d == atom_a &&
            parmset->properparms_vector_FF03ws_2[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_2[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_2[i].atom_a == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_2[i].atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_2[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_2[i].atom_b;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_2[i].atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_2[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_2[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_2[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_2[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_2[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_2[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_2[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_2[i].atom_a == "X"){

             proper_ready.atom_a=atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_2[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_2[i].atom_b;
             proper_ready.atom_d=atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_2[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_2[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_2[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_2[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_2[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_2[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_2[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_2[i].atom_a == "X"){


             proper_ready.atom_a=atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_2[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_2[i].atom_c;
             proper_ready.atom_d=atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_2[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_2[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_2[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_2[i].pn;

         }
     }
     return proper_ready;
    }
     else if(pn_index == 3){
     for(int i =0; i < parmset->properparms_vector_FF03ws_3.size(); i++){
        // cout  << parmset->properparms_vector_FF03ws_3[i].atom_a << " " << parmset->properparms_vector_FF03ws_3[i].atom_b << " "
        //       << parmset->properparms_vector_FF03ws_3[i].atom_c << " " << parmset->properparms_vector_FF03ws_3[i].atom_d << endl;



         if(parmset->properparms_vector_FF03ws_3[i].atom_a == atom_a &&
            parmset->properparms_vector_FF03ws_3[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_3[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_3[i].atom_d == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_3[i].atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_3[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_3[i].atom_c;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_3[i].atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_3[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_3[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_3[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_3[i].pn;

         }

         else if(parmset->properparms_vector_FF03ws_3[i].atom_d == atom_a &&
            parmset->properparms_vector_FF03ws_3[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_3[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_3[i].atom_a == atom_d){

             proper_ready.atom_a=parmset->properparms_vector_FF03ws_3[i].atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_3[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_3[i].atom_b;
             proper_ready.atom_d=parmset->properparms_vector_FF03ws_3[i].atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_3[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_3[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_3[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_3[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_3[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_3[i].atom_c == atom_b &&
            parmset->properparms_vector_FF03ws_3[i].atom_b == atom_c &&
            parmset->properparms_vector_FF03ws_3[i].atom_a == "X"){

             proper_ready.atom_a=atom_d;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_3[i].atom_c;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_3[i].atom_b;
             proper_ready.atom_d=atom_a;

             proper_ready.function=parmset->properparms_vector_FF03ws_3[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_3[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_3[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_3[i].pn;
         }
         else if(parmset->properparms_vector_FF03ws_3[i].atom_d == "X" &&
            parmset->properparms_vector_FF03ws_3[i].atom_c == atom_c &&
            parmset->properparms_vector_FF03ws_3[i].atom_b == atom_b &&
            parmset->properparms_vector_FF03ws_3[i].atom_a == "X"){


             proper_ready.atom_a=atom_a;
             proper_ready.atom_b=parmset->properparms_vector_FF03ws_3[i].atom_b;
             proper_ready.atom_c=parmset->properparms_vector_FF03ws_3[i].atom_c;
             proper_ready.atom_d=atom_d;

             proper_ready.function=parmset->properparms_vector_FF03ws_3[i].function;
             proper_ready.phase=parmset->properparms_vector_FF03ws_3[i].phase;
             proper_ready.kd=parmset->properparms_vector_FF03ws_3[i].kd;
             proper_ready.pn=parmset->properparms_vector_FF03ws_3[i].pn;

         }

    }
     return proper_ready;
   }

}



Protein::improper_dihedral_parm Protein::GetParm_improper_dihedrals(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset){


    improper_dihedral_parm proper_ready;
 //    cout <<  atom_a << " " << atom_b << " " << atom_c << " " << atom_d << endl;


    for(int i =0; i < parmset->improperparms_vector_FF99SB.size(); i++){
       // cout  << parmset->properparms_vector_FF99SB_3[i].atom_a << " " << parmset->properparms_vector_FF99SB_3[i].atom_b << " "
       //       << parmset->properparms_vector_FF99SB_3[i].atom_c << " " << parmset->properparms_vector_FF99SB_3[i].atom_d << endl;



        if(parmset->improperparms_vector_FF99SB[i].atom_a == atom_a &&
                 parmset->improperparms_vector_FF99SB[i].atom_b == atom_b &&
                 parmset->improperparms_vector_FF99SB[i].atom_c == atom_c &&
                 parmset->improperparms_vector_FF99SB[i].atom_d == atom_d){

                  proper_ready.atom_a=parmset->improperparms_vector_FF99SB[i].atom_a;
                  proper_ready.atom_b=parmset->improperparms_vector_FF99SB[i].atom_b;
                  proper_ready.atom_c=parmset->improperparms_vector_FF99SB[i].atom_c;
                  proper_ready.atom_d=parmset->improperparms_vector_FF99SB[i].atom_d;

                  proper_ready.function=parmset->improperparms_vector_FF99SB[i].function;
                  proper_ready.phase=parmset->improperparms_vector_FF99SB[i].phase;
                  proper_ready.kd=parmset->improperparms_vector_FF99SB[i].kd;
                  proper_ready.pn=parmset->improperparms_vector_FF99SB[i].pn;

              }

              else if(parmset->improperparms_vector_FF99SB[i].atom_d == atom_a &&
                 parmset->improperparms_vector_FF99SB[i].atom_c == atom_b &&
                 parmset->improperparms_vector_FF99SB[i].atom_b == atom_c &&
                 parmset->improperparms_vector_FF99SB[i].atom_a == atom_d){

                  proper_ready.atom_a=parmset->improperparms_vector_FF99SB[i].atom_d;
                  proper_ready.atom_b=parmset->improperparms_vector_FF99SB[i].atom_c;
                  proper_ready.atom_c=parmset->improperparms_vector_FF99SB[i].atom_b;
                  proper_ready.atom_d=parmset->improperparms_vector_FF99SB[i].atom_a;

                  proper_ready.function=parmset->improperparms_vector_FF99SB[i].function;
                  proper_ready.phase=parmset->improperparms_vector_FF99SB[i].phase;
                  proper_ready.kd=parmset->improperparms_vector_FF99SB[i].kd;
                  proper_ready.pn=parmset->improperparms_vector_FF99SB[i].pn;
              }
        else if(parmset->improperparms_vector_FF99SB[i].atom_a == "X" &&
                parmset->improperparms_vector_FF99SB[i].atom_b == atom_b &&
                parmset->improperparms_vector_FF99SB[i].atom_c == atom_c &&
                parmset->improperparms_vector_FF99SB[i].atom_d == atom_d){

                proper_ready.atom_a=atom_a;
                proper_ready.atom_b=atom_b;
                proper_ready.atom_c=atom_c;
                proper_ready.atom_d=atom_d;

                proper_ready.function=parmset->improperparms_vector_FF99SB[i].function;
                proper_ready.phase=parmset->improperparms_vector_FF99SB[i].phase;
                proper_ready.kd=parmset->improperparms_vector_FF99SB[i].kd;
                proper_ready.pn=parmset->improperparms_vector_FF99SB[i].pn;
      }
              else if(parmset->improperparms_vector_FF99SB[i].atom_a == "X" &&
                 parmset->improperparms_vector_FF99SB[i].atom_b == "X" &&
                 parmset->improperparms_vector_FF99SB[i].atom_c == atom_c &&
                 parmset->improperparms_vector_FF99SB[i].atom_d == atom_d){

                  proper_ready.atom_a=atom_a;
                  proper_ready.atom_b=atom_b;
                  proper_ready.atom_c=atom_c;
                  proper_ready.atom_d=atom_d;

                  proper_ready.function=parmset->improperparms_vector_FF99SB[i].function;
                  proper_ready.phase=parmset->improperparms_vector_FF99SB[i].phase;
                  proper_ready.kd=parmset->improperparms_vector_FF99SB[i].kd;
                  proper_ready.pn=parmset->improperparms_vector_FF99SB[i].pn;
              }





   }
    return proper_ready;
}





Protein::improper_dihedral_parm Protein::GetParm_improper_dihedrals_03ws(string atom_a, string atom_b, string atom_c, string atom_d, AMBER_parm_parser* parmset){


    improper_dihedral_parm proper_ready;
 //    cout <<  atom_a << " " << atom_b << " " << atom_c << " " << atom_d << endl;


    for(int i =0; i < parmset->improperparms_vector_FF03ws.size(); i++){
       // cout  << parmset->properparms_vector_FF03ws_3[i].atom_a << " " << parmset->properparms_vector_FF03ws_3[i].atom_b << " "
       //       << parmset->properparms_vector_FF03ws_3[i].atom_c << " " << parmset->properparms_vector_FF03ws_3[i].atom_d << endl;



        if(parmset->improperparms_vector_FF03ws[i].atom_a == atom_a &&
                 parmset->improperparms_vector_FF03ws[i].atom_b == atom_b &&
                 parmset->improperparms_vector_FF03ws[i].atom_c == atom_c &&
                 parmset->improperparms_vector_FF03ws[i].atom_d == atom_d){

                  proper_ready.atom_a=parmset->improperparms_vector_FF03ws[i].atom_a;
                  proper_ready.atom_b=parmset->improperparms_vector_FF03ws[i].atom_b;
                  proper_ready.atom_c=parmset->improperparms_vector_FF03ws[i].atom_c;
                  proper_ready.atom_d=parmset->improperparms_vector_FF03ws[i].atom_d;

                  proper_ready.function=parmset->improperparms_vector_FF03ws[i].function;
                  proper_ready.phase=parmset->improperparms_vector_FF03ws[i].phase;
                  proper_ready.kd=parmset->improperparms_vector_FF03ws[i].kd;
                  proper_ready.pn=parmset->improperparms_vector_FF03ws[i].pn;

              }

              else if(parmset->improperparms_vector_FF03ws[i].atom_d == atom_a &&
                 parmset->improperparms_vector_FF03ws[i].atom_c == atom_b &&
                 parmset->improperparms_vector_FF03ws[i].atom_b == atom_c &&
                 parmset->improperparms_vector_FF03ws[i].atom_a == atom_d){

                  proper_ready.atom_a=parmset->improperparms_vector_FF03ws[i].atom_d;
                  proper_ready.atom_b=parmset->improperparms_vector_FF03ws[i].atom_c;
                  proper_ready.atom_c=parmset->improperparms_vector_FF03ws[i].atom_b;
                  proper_ready.atom_d=parmset->improperparms_vector_FF03ws[i].atom_a;

                  proper_ready.function=parmset->improperparms_vector_FF03ws[i].function;
                  proper_ready.phase=parmset->improperparms_vector_FF03ws[i].phase;
                  proper_ready.kd=parmset->improperparms_vector_FF03ws[i].kd;
                  proper_ready.pn=parmset->improperparms_vector_FF03ws[i].pn;
              }
        else if(parmset->improperparms_vector_FF03ws[i].atom_a == "X" &&
                parmset->improperparms_vector_FF03ws[i].atom_b == atom_b &&
                parmset->improperparms_vector_FF03ws[i].atom_c == atom_c &&
                parmset->improperparms_vector_FF03ws[i].atom_d == atom_d){

                proper_ready.atom_a=atom_a;
                proper_ready.atom_b=atom_b;
                proper_ready.atom_c=atom_c;
                proper_ready.atom_d=atom_d;

                proper_ready.function=parmset->improperparms_vector_FF03ws[i].function;
                proper_ready.phase=parmset->improperparms_vector_FF03ws[i].phase;
                proper_ready.kd=parmset->improperparms_vector_FF03ws[i].kd;
                proper_ready.pn=parmset->improperparms_vector_FF03ws[i].pn;
      }
              else if(parmset->improperparms_vector_FF03ws[i].atom_a == "X" &&
                 parmset->improperparms_vector_FF03ws[i].atom_b == "X" &&
                 parmset->improperparms_vector_FF03ws[i].atom_c == atom_c &&
                 parmset->improperparms_vector_FF03ws[i].atom_d == atom_d){

                  proper_ready.atom_a=atom_a;
                  proper_ready.atom_b=atom_b;
                  proper_ready.atom_c=atom_c;
                  proper_ready.atom_d=atom_d;

                  proper_ready.function=parmset->improperparms_vector_FF03ws[i].function;
                  proper_ready.phase=parmset->improperparms_vector_FF03ws[i].phase;
                  proper_ready.kd=parmset->improperparms_vector_FF03ws[i].kd;
                  proper_ready.pn=parmset->improperparms_vector_FF03ws[i].pn;
              }





   }
    return proper_ready;
}


double Protein::assign_mass(string atomtype, AMBER_parm_parser* parm){
    for(int i =0; i < parm->atomic_parameters_FF99sb.size(); i++){
        if (parm->atomic_parameters_FF99sb[i].atom_type_from_parm == atomtype){
                return parm->atomic_parameters_FF99sb[i].atomic_mass;
        }
    }

}
void Protein::pair_assigner(vector < proper_dihedral_parm> proper_parameter_list){
    vector< vector <string> > pairs;
    vector< vector <string> > pairs_clean;
    vector < string > temp_pairs;
    for(int i =0; i< proper_parameter_list.size(); i++){
            //temp_pairs.push_back(proper_parameter_list[i].atom_a_index);
           // temp_pairs.push_back(proper_parameter_list[i].atom_b_index);
           // pairs.push_back(temp_pairs);
          // temp_pairs.clear();
          //  temp_pairs.push_back(proper_parameter_list[i].atom_a_index);
          // temp_pairs.push_back(proper_parameter_list[i].atom_c_index);
           // pairs.push_back(temp_pairs);
           // temp_pairs.clear();
            temp_pairs.push_back(proper_parameter_list[i].atom_a_index);
            temp_pairs.push_back(proper_parameter_list[i].atom_d_index);
            pairs.push_back(temp_pairs);
            pairs_clean.push_back(temp_pairs);

            temp_pairs.clear();
    }




    for(int i =0; i< pairs.size(); i++){
        for(int j =i+1; j< pairs_clean.size(); j++){
            if(pairs[i][0] == pairs_clean[j][0] && pairs[i][1] == pairs_clean[j][1]){
                pairs_clean[j][0]="A";
                pairs_clean[j][1]="A";

            }
            else if(pairs[i][1] == pairs_clean[j][0] && pairs[i][0] == pairs_clean[j][1]){
                pairs_clean[j][0]="A";
                pairs_clean[j][1]="A";
            }
        }

    }
    pairs.clear();



    for(int j =0; j< pairs_clean.size(); j++){
        if(pairs_clean[j][0] != "A" && pairs_clean[j][1] != "A"){
            temp_pairs.push_back(pairs_clean[j][0]);
            temp_pairs.push_back(pairs_clean[j][1]);
            pairs.push_back(temp_pairs);
            temp_pairs.clear();

        }
    }


    this->pairs_ready = pairs;


}


void Protein::parametrizer(string dssp_file){
   AMBER_parm_parser* parm = new AMBER_parm_parser();
   string dssp_name = dssp_file;
   DSSP_parser* dssp= new DSSP_parser(dssp_name.c_str());


    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        for(int j =0; j < this->aminoacid_parameters[i]->atom_types.size(); j++){
                  this->indexial_masses.push_back(this->assign_mass(this->aminoacid_parameters[i]->atom_types[j],parm));
        }
    }



    int bond_index =0;
    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        if(dssp->is_structured[i] == "N"){
        for(int j =0; j < this->aminoacid_parameters[i]->bonds.size(); j++){

            //cout << this->aminoacid_parameters[i]->aminoacid_type << " ";
            //cout << this->aminoacid_parameters[i]->bonds.size() << endl;
            //cout << this->aminoacid_parameters[i]->bonds_indexial[j][0] << " " << this->aminoacid_parameters[i]->bonds_indexial[j][1] << endl;
            //cout << this->aminoacid_parameters[i]->bonds[j][0] << " " << this->aminoacid_parameters[i]->bonds[j][1] << endl;

           bonds_parameter_list.push_back(this->GetParm_bonds_03ws(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->bonds_indexial[j][0])-1],
                                          this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->bonds_indexial[j][1])-1], parm));
           bonds_parameter_list[bond_index].atom_a_index = this->aminoacid_parameters[i]->bonds_indexial[j][0];
           bonds_parameter_list[bond_index].atom_b_index = this->aminoacid_parameters[i]->bonds_indexial[j][1];

           bond_index++;


        }
    }
        else if(dssp->is_structured[i] == "Y"){
        for(int j =0; j < this->aminoacid_parameters[i]->bonds.size(); j++){

            //cout << this->aminoacid_parameters[i]->aminoacid_type << " ";
            //cout << this->aminoacid_parameters[i]->bonds.size() << endl;
            //cout << this->aminoacid_parameters[i]->bonds_indexial[j][0] << " " << this->aminoacid_parameters[i]->bonds_indexial[j][1] << endl;
            //cout << this->aminoacid_parameters[i]->bonds[j][0] << " " << this->aminoacid_parameters[i]->bonds[j][1] << endl;

           bonds_parameter_list.push_back(this->GetParm_bonds(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->bonds_indexial[j][0])-1],
                                          this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->bonds_indexial[j][1])-1], parm));
           bonds_parameter_list[bond_index].atom_a_index = this->aminoacid_parameters[i]->bonds_indexial[j][0];
           bonds_parameter_list[bond_index].atom_b_index = this->aminoacid_parameters[i]->bonds_indexial[j][1];

           bond_index++;


        }
    }
    }
    /*
    for(int i =0; i < bonds_parameter_list.size(); i++){
        cout << bonds_parameter_list[i].atom_a << " " << bonds_parameter_list[i].atom_b << " " << bonds_parameter_list[i].b0 << endl;
    }
    */
    int angle_index = 0;
    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        if(dssp->is_structured[i] == "N"){

        for(int j =0; j < this->aminoacid_parameters[i]->angles_indexial.size(); j++){

          // cout << this->aminoacid_parameters[i]->aminoacid_type << " ";
         //  cout << this->aminoacid_parameters[i]->angles_indexial.size() << endl;
          // cout << this->aminoacid_parameters[i]->angles_indexial[j][0] << " " << this->aminoacid_parameters[i]->angles_indexial[j][1] << endl;
            //cout << this->aminoacid_parameters[i]->angles_indexial.size() << " " << this->aminoacid_parameters[i]->angles.size() << endl;
            if(stoi(this->aminoacid_parameters[i]->angles_indexial[j][0])-1 >= 0 ||
               stoi(this->aminoacid_parameters[i]->angles_indexial[j][1])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->angles_indexial[j][2])-1 >= 0 )  {
             angles_parameter_list.push_back(this->GetParm_angles_03ws(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->angles_indexial[j][0])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->angles_indexial[j][1])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->angles_indexial[j][2])-1], parm));

                    angles_parameter_list[angle_index].atom_a_index=this->aminoacid_parameters[i]->angles_indexial[j][0];
                    angles_parameter_list[angle_index].atom_b_index=this->aminoacid_parameters[i]->angles_indexial[j][1];
                    angles_parameter_list[angle_index].atom_c_index=this->aminoacid_parameters[i]->angles_indexial[j][2];


             /*
             cout << "this is the angle list "<< this->aminoacid_parameters[i]->angles_indexial[j][0] << " " <<
                    this->aminoacid_parameters[i]->angles_indexial[j][1] << " " <<
                    this->aminoacid_parameters[i]->angles_indexial[j][2] << " " << endl;
            */
             angle_index++;
            }

        }
        }
        else if(dssp->is_structured[i] == "Y"){

        for(int j =0; j < this->aminoacid_parameters[i]->angles_indexial.size(); j++){

          // cout << this->aminoacid_parameters[i]->aminoacid_type << " ";
         //  cout << this->aminoacid_parameters[i]->angles_indexial.size() << endl;
          // cout << this->aminoacid_parameters[i]->angles_indexial[j][0] << " " << this->aminoacid_parameters[i]->angles_indexial[j][1] << endl;
            //cout << this->aminoacid_parameters[i]->angles_indexial.size() << " " << this->aminoacid_parameters[i]->angles.size() << endl;
            if(stoi(this->aminoacid_parameters[i]->angles_indexial[j][0])-1 >= 0 ||
               stoi(this->aminoacid_parameters[i]->angles_indexial[j][1])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->angles_indexial[j][2])-1 >= 0 )  {
             angles_parameter_list.push_back(this->GetParm_angles(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->angles_indexial[j][0])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->angles_indexial[j][1])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->angles_indexial[j][2])-1], parm));

                    angles_parameter_list[angle_index].atom_a_index=this->aminoacid_parameters[i]->angles_indexial[j][0];
                    angles_parameter_list[angle_index].atom_b_index=this->aminoacid_parameters[i]->angles_indexial[j][1];
                    angles_parameter_list[angle_index].atom_c_index=this->aminoacid_parameters[i]->angles_indexial[j][2];



                    /*
                    cout << "this is the angle list "<< this->aminoacid_parameters[i]->angles_indexial[j][0] << " " <<
                           this->aminoacid_parameters[i]->angles_indexial[j][1] << " " <<
                           this->aminoacid_parameters[i]->angles_indexial[j][2] << " " << endl;
                   */

             angle_index++;
            }

        }
        }

    }

    //dihedral degenerecency_killer 2 - bloodbath return
    vector < vector < string > > proper_quartets_in_series;
   vector < string > temp_proper_quartet;
   vector < int > to_delete;

    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        for(int j =0; j < this->aminoacid_parameters[i]->propers_indexial.size(); j++){
            temp_proper_quartet.push_back(this->aminoacid_parameters[i]->propers_indexial[j][0]);
            temp_proper_quartet.push_back(this->aminoacid_parameters[i]->propers_indexial[j][1]);
            temp_proper_quartet.push_back(this->aminoacid_parameters[i]->propers_indexial[j][2]);
            temp_proper_quartet.push_back(this->aminoacid_parameters[i]->propers_indexial[j][3]);
            proper_quartets_in_series.push_back(temp_proper_quartet);
            temp_proper_quartet.clear();
            to_delete.push_back(0);
        }
    }
    vector < vector < string > > proper_quartets_in_series_2 = proper_quartets_in_series;
    for(int i =0; i < proper_quartets_in_series_2.size(); i++){
        for(int j =i; j< proper_quartets_in_series.size(); j++){

        if(proper_quartets_in_series_2[i][0] == proper_quartets_in_series[j][3] &&
           proper_quartets_in_series_2[i][1] == proper_quartets_in_series[j][2] &&
           proper_quartets_in_series_2[i][2] == proper_quartets_in_series[j][1] &&
           proper_quartets_in_series_2[i][3] == proper_quartets_in_series[j][0]){
            to_delete[i]++;
        }
        else if(proper_quartets_in_series_2[i][0] == proper_quartets_in_series[j][0] &&
           proper_quartets_in_series_2[i][1] == proper_quartets_in_series[j][1] &&
           proper_quartets_in_series_2[i][2] == proper_quartets_in_series[j][2] &&
           proper_quartets_in_series_2[i][3] == proper_quartets_in_series[j][3]){
            to_delete[i]++;
        }
        }
    }



    int total_temp=0;
        for(int i =0; i < this->aminoacid_parameters.size(); i++){

            for(int j =0; j < this->aminoacid_parameters[i]->propers_indexial.size(); j++){
                total_temp = total_temp+1;
               // cout << to_delete.size() << " " << to_delete[total_temp-1] << " " << i  <<endl;
               // cout << "lol" << endl;
                if(to_delete[total_temp-1] > 1){
                    this->aminoacid_parameters[i]->propers_indexial[j][0] = "0";
                    this->aminoacid_parameters[i]->propers_indexial[j][1] = "0";
                    this->aminoacid_parameters[i]->propers_indexial[j][2] = "0";
                    this->aminoacid_parameters[i]->propers_indexial[j][3] = "0";

                }
            }
        }





    int size =0;
    proper_dihedral_parm proper_ready;

    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        size = size + int(this->aminoacid_parameters[i]->propers_indexial.size());
        if(dssp->is_structured[i] == "N"){
        for(int j =0; j < this->aminoacid_parameters[i]->propers_indexial.size(); j++){
          // cout << this->aminoacid_parameters[i]->aminoacid_type << " ";
           // cout << "propers size    " << size << endl;
          // cout << this->aminoacid_parameters[i]->angles_indexial[j][0] << " " << this->aminoacid_parameters[i]->angles_indexial[j][1] << endl;
            //cout << this->aminoacid_parameters[i]->angles_indexial.size() << " " << this->aminoacid_parameters[i]->angles.size() << endl;

            if(stoi(this->aminoacid_parameters[i]->propers_indexial[j][0])-1 >= 0 ||
               stoi(this->aminoacid_parameters[i]->propers_indexial[j][1])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->propers_indexial[j][2])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->propers_indexial[j][3])-1 >= 0)  {
            for(int pn_index =0 ; pn_index < 4 ; pn_index++){
             proper_ready = (this->GetParm_dihedrals_03ws(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][0])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][1])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][2])-1],
                                              this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][3])-1], parm, pn_index, this->aminoacid_parameters[i]->aminoacid_type ));
                if(proper_ready.atom_a != ""){
                    proper_ready.atom_a_index = this->aminoacid_parameters[i]->propers_indexial[j][0];
                    proper_ready.atom_b_index = this->aminoacid_parameters[i]->propers_indexial[j][1];
                    proper_ready.atom_c_index = this->aminoacid_parameters[i]->propers_indexial[j][2];
                    proper_ready.atom_d_index = this->aminoacid_parameters[i]->propers_indexial[j][3];
                    proper_parameter_list.push_back(proper_ready);
                }
                }

            }

        }
        }
        else if(dssp->is_structured[i] == "Y"){
        for(int j =0; j < this->aminoacid_parameters[i]->propers_indexial.size(); j++){
          // cout << this->aminoacid_parameters[i]->aminoacid_type << " ";
           // cout << "propers size    " << size << endl;
          // cout << this->aminoacid_parameters[i]->angles_indexial[j][0] << " " << this->aminoacid_parameters[i]->angles_indexial[j][1] << endl;
            //cout << this->aminoacid_parameters[i]->angles_indexial.size() << " " << this->aminoacid_parameters[i]->angles.size() << endl;

            if(stoi(this->aminoacid_parameters[i]->propers_indexial[j][0])-1 >= 0 ||
               stoi(this->aminoacid_parameters[i]->propers_indexial[j][1])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->propers_indexial[j][2])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->propers_indexial[j][3])-1 >= 0)  {
            for(int pn_index =0 ; pn_index < 7 ; pn_index++){
             proper_ready = (this->GetParm_dihedrals(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][0])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][1])-1],
                                             this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][2])-1],
                                              this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->propers_indexial[j][3])-1], parm, pn_index, this->aminoacid_parameters[i]->aminoacid_type ));
                if(proper_ready.atom_a != ""){
                    proper_ready.atom_a_index = this->aminoacid_parameters[i]->propers_indexial[j][0];
                    proper_ready.atom_b_index = this->aminoacid_parameters[i]->propers_indexial[j][1];
                    proper_ready.atom_c_index = this->aminoacid_parameters[i]->propers_indexial[j][2];
                    proper_ready.atom_d_index = this->aminoacid_parameters[i]->propers_indexial[j][3];
                    proper_parameter_list.push_back(proper_ready);
                }
                }

            }

        }
        }

    }

/*
    for(int i =0; i < proper_parameter_list.size(); i++){

    cout << "this is the dihedral list "<< proper_parameter_list[i].atom_a << " " <<
           proper_parameter_list[i].atom_b << " " <<
           proper_parameter_list[i].atom_c << " " <<
           proper_parameter_list[i].atom_d <<  " " <<
            proper_parameter_list[i].atom_a_index <<  " " <<
            proper_parameter_list[i].atom_b_index <<  " " <<
            proper_parameter_list[i].atom_c_index <<  " " <<
            proper_parameter_list[i].atom_d_index <<  " " <<
            proper_parameter_list[i].kd <<  " " <<
           proper_parameter_list[i].pn << endl;

        }
*/
    improper_dihedral_parm improper_ready;


    for(int i =0; i < this->aminoacid_parameters.size(); i++){
        if(dssp->is_structured[i] == "Y"){

        for(int j =0; j < this->aminoacid_parameters[i]->impropers_indexial.size(); j++){
            /*
            cout << this->aminoacid_parameters[i]->impropers_indexial[j][0] << " "
                 << this->aminoacid_parameters[i]->impropers_indexial[j][1] << " "
                 << this->aminoacid_parameters[i]->impropers_indexial[j][2] << " "
                 << this->aminoacid_parameters[i]->impropers_indexial[j][3]  << endl;
            cout << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][0])-1] << " "
                 << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][1])-1] << " "
                 << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][2])-1] << " "
                 << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][3])-1] << endl;*/

            if(stoi(this->aminoacid_parameters[i]->impropers_indexial[j][0])-1 >= 0 ||
               stoi(this->aminoacid_parameters[i]->impropers_indexial[j][1])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->impropers_indexial[j][2])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->impropers_indexial[j][3])-1 >= 0)  {
               improper_ready = (this->GetParm_improper_dihedrals(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][0])-1],
                                                this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][1])-1],
                                                this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][2])-1],
                                                 this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][3])-1], parm));

               improper_ready.atom_a_index = this->aminoacid_parameters[i]->impropers_indexial[j][0];
               improper_ready.atom_b_index = this->aminoacid_parameters[i]->impropers_indexial[j][1];
               improper_ready.atom_c_index = this->aminoacid_parameters[i]->impropers_indexial[j][2];
               improper_ready.atom_d_index = this->aminoacid_parameters[i]->impropers_indexial[j][3];
               improper_parameter_list.push_back(improper_ready);

            }


        }
        }
        else if(dssp->is_structured[i] == "N"){

        for(int j =0; j < this->aminoacid_parameters[i]->impropers_indexial.size(); j++){
            /*
            cout << this->aminoacid_parameters[i]->impropers_indexial[j][0] << " "
                 << this->aminoacid_parameters[i]->impropers_indexial[j][1] << " "
                 << this->aminoacid_parameters[i]->impropers_indexial[j][2] << " "
                 << this->aminoacid_parameters[i]->impropers_indexial[j][3]  << endl;
            cout << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][0])-1] << " "
                 << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][1])-1] << " "
                 << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][2])-1] << " "
                 << this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][3])-1] << endl;*/

            if(stoi(this->aminoacid_parameters[i]->impropers_indexial[j][0])-1 >= 0 ||
               stoi(this->aminoacid_parameters[i]->impropers_indexial[j][1])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->impropers_indexial[j][2])-1 >= 0  ||
               stoi(this->aminoacid_parameters[i]->impropers_indexial[j][3])-1 >= 0)  {
               improper_ready = (this->GetParm_improper_dihedrals_03ws(this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][0])-1],
                                                this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][1])-1],
                                                this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][2])-1],
                                                 this->atomtype_serialized[stoi(this->aminoacid_parameters[i]->impropers_indexial[j][3])-1], parm));

               improper_ready.atom_a_index = this->aminoacid_parameters[i]->impropers_indexial[j][0];
               improper_ready.atom_b_index = this->aminoacid_parameters[i]->impropers_indexial[j][1];
               improper_ready.atom_c_index = this->aminoacid_parameters[i]->impropers_indexial[j][2];
               improper_ready.atom_d_index = this->aminoacid_parameters[i]->impropers_indexial[j][3];
               improper_parameter_list.push_back(improper_ready);

            }


        }
        }
    }


    /*
    for(int i =0; i < improper_parameter_list.size(); i++){


    cout << "this is the improper dihedral list "<< improper_parameter_list[i].atom_a << " " <<
           improper_parameter_list[i].atom_b << " " <<
           improper_parameter_list[i].atom_c << " " <<
           improper_parameter_list[i].atom_d <<  " " <<
           improper_parameter_list[i].kd <<  " " <<
           improper_parameter_list[i].pn << endl;

        }
        */
    this->pair_assigner(this->proper_parameter_list);


}

//function to connect and generate inter-residue parameters
//void  Protein::connector(){








































