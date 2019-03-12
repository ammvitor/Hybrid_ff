#include "aminoacid.h"

aminoacid::aminoacid(string type, int C_or_N, int starting_atomic_index, string forfield_type)
{
    //define the name of the residue and the list - if ILE,ASP,ASN or LEU, assign in a separated set, for the AMBER99SB-ildn speciual torsions
    this->C_or_N = C_or_N;
    this->starting_atomic_index = starting_atomic_index;
    if(forfield_type == "99sb"){
    if(C_or_N == 1){
    type = "N"+type;
    if(type == "NILE"){
        this->SPECIAL_PARM(type);
    }
    else if(type == "NLEU"){
        this->SPECIAL_PARM(type);

    }
    else if(type == "NASP"){
        this->SPECIAL_PARM(type);

    }
    else if(type == "NASN"){
        this->SPECIAL_PARM(type);

    }
    else{
        this->PARM(type);
    }

    }
    else if(C_or_N == 2){
    if(type == "ILE"){
        this->SPECIAL_PARM(type);
    }
    else if(type == "LEU"){
        this->SPECIAL_PARM(type);

    }
    else if(type == "ASP"){
        this->SPECIAL_PARM(type);

    }
    else if(type == "ASN"){
        this->SPECIAL_PARM(type);

    }
    else{
        this->PARM(type);
    }

    }
    else if(C_or_N == 3){
    type = "C"+type;
    if(type == "CILE"){
        this->SPECIAL_PARM(type);
    }
    else if(type == "CLEU"){
        this->SPECIAL_PARM(type);

    }
    else if(type == "CASP"){
        this->SPECIAL_PARM(type);

    }
    else if(type == "CASN"){
        this->SPECIAL_PARM(type);

    }
    else{
        this->PARM(type);
    }
    }
    }



    if(forfield_type == "03ws"){
    if(C_or_N == 1){
    type = "N"+type;

        this->PARM_03ws(type);

    }
    else if(C_or_N == 2){

        this->PARM_03ws(type);


    }
    else if(C_or_N == 3){
    type = "C"+type;

        this->PARM_03ws(type);

    }
    }
}

vector < vector < string >> add_previous_backbone(vector < vector < string >> bonds){
    //add virtual bonds for the middle residues and the Cterminal, to be able to generate the tree from the residue adding the angles and torsions automatically.
    vector < string > bond_temp;
    bond_temp.push_back("-C");
    bond_temp.push_back("-O");
    bonds.push_back(bond_temp);
    bond_temp.clear();
    bond_temp.push_back("-C");
    bond_temp.push_back("-CA");
    bonds.push_back(bond_temp);
    bond_temp.clear();
    bond_temp.push_back("-CA");
    bond_temp.push_back("-HA");
    bonds.push_back(bond_temp);
    bond_temp.clear();
    bond_temp.push_back("-CA");
    bond_temp.push_back("-R");
    bonds.push_back(bond_temp);
    bond_temp.clear();
    bond_temp.push_back("-CA");
    bond_temp.push_back("-N");
    bonds.push_back(bond_temp);
    bond_temp.clear();
    return bonds;
}

//a boolean compararer or vectors, to delete equal vectors.
bool equal_vector(std::vector<int> a, std::vector<int> b) {
  std::sort(a.begin(), a.end());
  std::sort(b.begin(), b.end());

  return std::equal(a.begin(), a.end(), b.begin());
}

static bool IsEqual(vector<int> a, vector<int> b)
{
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    return (a == b);
}
static bool IsEqual(vector<string> a, vector<string> b)
{
    sort(a.begin(), a.end());
    sort(b.begin(), b.end());
    return (a == b);
}


void aminoacid::translator_name_to_index_bonds(){

    /*
     * To organize the topology, gormacs accepts the indexes of the attoms listed in the [atoms] directive
     * therefore, we need to assign the right indexes to write the topology
     *
     */
    vector < string > bond_indexal_temp;
    string N_index;
    // translated the bonds in the aminoacids.rtp to its respective indexes
    for(int i = 0; i < this->bonds.size(); i++){
        for(int j = 0; j < this->atom_names.size(); j++){
        if(this->atom_names[j] == "N"){
           N_index = to_string(this->atomic_indexes[j]);

        }

        if(this->bonds[i][0] == this->atom_names[j]){

            for(int k = 0; k < this->atom_names.size(); k++){
            if(this->bonds[i][1] == this->atom_names[k]){

                bond_indexal_temp.push_back(to_string(this->atomic_indexes[j]));
                bond_indexal_temp.push_back(to_string(this->atomic_indexes[k]));
                this->bonds_indexial.push_back(bond_indexal_temp);
                bond_indexal_temp.clear();
            }
            }
        }

    }

}
//for the previous bonds, add a -1 for -C, adds a -2 for -O, -3 for -CA, -4 for -HA, -5 for -R(the first atom in the sidechain), -6 for the -N ;
    for(int i = 0; i < this->bonds.size(); i++){

        if(this->bonds[i][0] == "-C" || this->bonds[i][0] == "-O" ||  this->bonds[i][0] == "-CA" || this->bonds[i][0] == "-HA" || this->bonds[i][0] == "-R" || this->bonds[i][0] == "-N" || this->bonds[i][0] == "N" ){

            if(this->bonds[i][1] == "-C" || this->bonds[i][1] == "-O" ||  this->bonds[i][1] == "-CA" || this->bonds[i][1] == "-HA" || this->bonds[i][1] == "-R" || this->bonds[i][1] == "-N" || this->bonds[i][1] == "N" ){

                if(this->bonds[i][0] == "-C"){
                    bond_indexal_temp.push_back("-1");
                }
                else if(this->bonds[i][0] == "-O"){
                    bond_indexal_temp.push_back("-2");

                }
                else if(this->bonds[i][0] == "-CA"){
                    bond_indexal_temp.push_back("-3");

                }
                else if(this->bonds[i][0] == "-HA"){
                    bond_indexal_temp.push_back("-4");

                }
                else if(this->bonds[i][0] == "-R"){
                    bond_indexal_temp.push_back("-5");

                }
                else if(this->bonds[i][0] == "N"){
                    bond_indexal_temp.push_back(N_index);

                }
                else if(this->bonds[i][0] == "-N"){
                    bond_indexal_temp.push_back("-6");

                }

                if(this->bonds[i][1] == "-C"){
                    bond_indexal_temp.push_back("-1");
                }
                else if(this->bonds[i][1] == "-O"){
                    bond_indexal_temp.push_back("-2");

                }
                else if(this->bonds[i][1] == "-CA"){
                    bond_indexal_temp.push_back("-3");

                }
                else if(this->bonds[i][1] == "-HA"){
                    bond_indexal_temp.push_back("-4");

                }
                else if(this->bonds[i][1] == "-R"){
                    bond_indexal_temp.push_back("-5");

                }
                else if(this->bonds[i][1] == "N"){
                    bond_indexal_temp.push_back(N_index);

                }
                else if(this->bonds[i][1] == "-N"){
                    bond_indexal_temp.push_back("-6");

                }

                this->bonds_indexial.push_back(bond_indexal_temp);
                bond_indexal_temp.clear();

            }

        }
    }
        //forthe liste improter in aminoacids.rpt. creates a vector with the index instead of the atomname ;
        // add a -1 for -C, adds a -2 for -O, -3 for -CA, -4 for -HA, -5 for -R(the first atom in the sidechain), -6 for the -N
            for(int i = 0; i < this->impropers.size(); i++){
                for(int j = 0; j < 4; j++){

                for(int k =0 ; k < this->atom_names.size(); k++){

                        if(this->impropers[i][j] == this->atom_names[k]){
                            bond_indexal_temp.push_back(to_string(this->atomic_indexes[k]));
                            //cout << "1 " <<  this->atom_names[k]  << endl;
                        }
                        else if(this->impropers[i][j] == "-O"){
                            bond_indexal_temp.push_back("-2");
                                    //cout << "-2" << endl;
                                    break;
                        }
                        else if(this->impropers[i][j] == "-C"){
                            bond_indexal_temp.push_back("-1");
                                   // cout << "-1" << endl;
                                    break;

                        }
                        else if(this->impropers[i][j] == "-CA"){
                            bond_indexal_temp.push_back("-3");
                                            //cout << "-3" << endl;
                                            break;


                        }
                        else if(this->impropers[i][j] == "-HA"){
                            bond_indexal_temp.push_back("-4");
                            //cout << "-4" << endl;
                                break;


                        }
                        else if(this->impropers[i][j] == "-N"){
                            bond_indexal_temp.push_back("-6");
                           // cout << "-6" << endl;
                            break;


                       }
                       else if(this->impropers[i][j] == "+N"){
                            bond_indexal_temp.push_back("-7");
                            //cout << "-7" << endl;
                            break;

                        }


                    }
                    }

                this->impropers_indexial.push_back(bond_indexal_temp);
                bond_indexal_temp.clear();
                }


            /*
            cout <<  this->impropers_indexial.size() << endl;
            for(int i = 0; i < this->impropers_indexial.size(); i++){

              cout << " impropers translation " << endl;
              cout << this->impropers_indexial[i][0] << " " << this->impropers_indexial[i][1] << " " << this->impropers_indexial[i][2] << " " << this->impropers_indexial[i][3] << " " << endl;
              cout << this->impropers[i][0] << " " << this->impropers[i][1] << " " << this->impropers[i][2] << " " << this->impropers[i][3] << " " << endl;

            }
            */

}

void aminoacid::define_angles(vector <vector < string > > bonds){
    //Angles Finder by buiding the tree from the bonds assigned from aminoacids.rtp, will output a matrix with the names

    vector <vector < string> > atom_id_angles;
    vector < string> ang_Atom_triad;
      //Angles Finder
      for(int i =0; i < bonds.size(); i++){
          for(int j = 0 ; j < bonds.size(); j++){
              if(i!=j){

                       if(bonds[i][0] == bonds[j][0]){
                           ang_Atom_triad.push_back(bonds[i][1]);
                           ang_Atom_triad.push_back(bonds[i][0]);
                           ang_Atom_triad.push_back(bonds[j][1]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }
                       if(bonds[i][0] == bonds[j][1]){
                           ang_Atom_triad.push_back(bonds[j][0]);
                           ang_Atom_triad.push_back(bonds[i][0]);
                           ang_Atom_triad.push_back(bonds[i][1]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }
                       if(bonds[i][1] == bonds[j][0]){
                           ang_Atom_triad.push_back(bonds[j][1]);
                           ang_Atom_triad.push_back(bonds[j][0]);
                           ang_Atom_triad.push_back(bonds[i][0]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }
                       if(bonds[i][1] == bonds[j][1]){
                           ang_Atom_triad.push_back(bonds[i][0]);
                           ang_Atom_triad.push_back(bonds[i][1]);
                           ang_Atom_triad.push_back(bonds[j][0]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }


              }

          }
      }
      vector <vector < string> > atom_id_angles_2 = atom_id_angles;
      vector <vector < string> > atom_id_angles_3;

      //Angles Degenerescency killer
      while(atom_id_angles.size() > 0){

          for(int k = 1; k < atom_id_angles_2.size(); k++){


             if(IsEqual(atom_id_angles[0],atom_id_angles_2[k])){

             atom_id_angles_3.push_back(atom_id_angles[0]);
             atom_id_angles.erase(atom_id_angles.begin());
             atom_id_angles_2.erase(atom_id_angles_2.begin());

             atom_id_angles_2.erase(atom_id_angles_2.begin()+k-1);
             atom_id_angles.erase(atom_id_angles.begin()+k-1);
             break;

         }

       }

     }
      /*
      for(int i =0; i < atom_id_angles_3.size(); i++){
          cout << " inthe difine anglesssssssssssssssssssss " << atom_id_angles_3[i][0] << " " << atom_id_angles_3[i][1] << " " <<  atom_id_angles_3[i][2] << endl;
      }
      */
      if(this->C_or_N==2 || this->C_or_N==3){
          atom_id_angles_3.pop_back();
      }
      /*
      for(int i =0; i < atom_id_angles_3.size(); i++){
          cout << " inthe difine poppeesss " << atom_id_angles_3[i][0] << " " << atom_id_angles_3[i][1] << " " <<  atom_id_angles_3[i][2] << endl;
      }*/

      this->angles = atom_id_angles_3;


}

void aminoacid::define_angles_indexial(vector <vector < string > > bonds){
    vector <vector < string> > atom_id_angles;
    vector < string> ang_Atom_triad;
    //Angles Finder by buiding the tree from the bonds assigned from aminoacids.rtp, will output a matrix with the indexes
      for(int i =0; i < bonds.size(); i++){
          for(int j = 0 ; j < bonds.size(); j++){
              if(i!=j){

                       if(bonds[i][0] == bonds[j][0]){
                           ang_Atom_triad.push_back(bonds[i][1]);
                           ang_Atom_triad.push_back(bonds[i][0]);
                           ang_Atom_triad.push_back(bonds[j][1]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }
                       if(bonds[i][0] == bonds[j][1]){
                           ang_Atom_triad.push_back(bonds[j][0]);
                           ang_Atom_triad.push_back(bonds[i][0]);
                           ang_Atom_triad.push_back(bonds[i][1]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }
                       if(bonds[i][1] == bonds[j][0]){
                           ang_Atom_triad.push_back(bonds[j][1]);
                           ang_Atom_triad.push_back(bonds[j][0]);
                           ang_Atom_triad.push_back(bonds[i][0]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }
                       if(bonds[i][1] == bonds[j][1]){
                           ang_Atom_triad.push_back(bonds[i][0]);
                           ang_Atom_triad.push_back(bonds[i][1]);
                           ang_Atom_triad.push_back(bonds[j][0]);
                           atom_id_angles.push_back(ang_Atom_triad);
                           ang_Atom_triad.clear();

                       }


              }

          }
      }
      vector <vector < string> > atom_id_angles_2 = atom_id_angles;
      vector <vector < string> > atom_id_angles_3;

      //Angles Degenerescency killer
      while(atom_id_angles.size() > 0){

          for(int k = 1; k < atom_id_angles_2.size(); k++){


             if(IsEqual(atom_id_angles[0],atom_id_angles_2[k])){

             atom_id_angles_3.push_back(atom_id_angles[0]);
             atom_id_angles.erase(atom_id_angles.begin());
             atom_id_angles_2.erase(atom_id_angles_2.begin());

             atom_id_angles_2.erase(atom_id_angles_2.begin()+k-1);
             atom_id_angles.erase(atom_id_angles.begin()+k-1);
             break;

         }

       }

     }
      if(this->C_or_N==2 || this->C_or_N==3){
          atom_id_angles_3.pop_back();
      }
      this->angles_indexial = atom_id_angles_3;


}

void aminoacid::define_proper_torsions(vector <vector < string > > bonds){
    vector <vector < string> > atom_id_dihedral;
    vector < string> ang_Atom_quad;



    //dihedral Finder by buiding the tree from the bonds assigned from aminoacids.rtp, will output a matrix with the atom names
    for(int i =0; i < bonds.size(); i++){
        for(int j = 0 ; j < bonds.size(); j++){
            for(int k = 0; k < bonds.size(); k++){


            if(i!=k && k!=j && j!=i){

                     if(bonds[i][0] == bonds[j][0]){

                         if(bonds[i][1] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][1] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }


                     }


                     if(bonds[i][0] == bonds[j][1]){

                         if(bonds[i][1] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][1] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }

                     }


                     if(bonds[i][1] == bonds[j][0]){

                         if(bonds[i][0] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][0] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][1]){

                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }


                     }


                     if(bonds[i][1] == bonds[j][1]){

                         if(bonds[i][0] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][0] == bonds[k][1]){

                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][1]){

                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();
                     }
            }



            }

        }
      }
    }

    vector <vector < string> > atom_id_dihedral_2 = atom_id_dihedral;
    vector <vector < string> > atom_id_dihedral_3;
    vector <int> del_index;

    //Dihedral Degenerescency killer


    while(atom_id_dihedral.size() > 0){

        for(int k = 1; k < atom_id_dihedral.size(); k++){


            if(IsEqual(atom_id_dihedral_2[0],atom_id_dihedral[k])){
             del_index.push_back(k);

            }

       }

         if(del_index.size() > 0){
             atom_id_dihedral_3.push_back(atom_id_dihedral[0]);
             atom_id_dihedral_2.erase(atom_id_dihedral_2.begin());
             atom_id_dihedral.erase(atom_id_dihedral.begin());
         for(int i =0; i < del_index.size(); i++){
            atom_id_dihedral.erase(atom_id_dihedral.begin()+del_index[i]-i-1);
             atom_id_dihedral_2.erase(atom_id_dihedral_2.begin()+del_index[i]-i-1);
        }
         }
         else{
             atom_id_dihedral_3.push_back(atom_id_dihedral[0]);
             atom_id_dihedral_2.erase(atom_id_dihedral_2.begin());
             atom_id_dihedral.erase(atom_id_dihedral.begin());
         }

         del_index.clear();


     }

    this->propers = atom_id_dihedral_3;

}

void aminoacid::define_proper_torsions_indexial(vector <vector < string > > bonds){
    vector <vector < string> > atom_id_dihedral;
    vector < string> ang_Atom_quad;



    //Dihedral Finder
    for(int i =0; i < bonds.size(); i++){
        for(int j = 0 ; j < bonds.size(); j++){
            for(int k = 0; k < bonds.size(); k++){


            if(i!=k && k!=j && j!=i){

                     if(bonds[i][0] == bonds[j][0]){

                         if(bonds[i][1] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][1] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }


                     }


                     if(bonds[i][0] == bonds[j][1]){

                         if(bonds[i][1] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][1] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][0]){
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }

                     }


                     if(bonds[i][1] == bonds[j][0]){

                         if(bonds[i][0] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][0] == bonds[k][1]){
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][1] == bonds[k][1]){

                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }


                     }


                     if(bonds[i][1] == bonds[j][1]){

                         if(bonds[i][0] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[i][0] == bonds[k][1]){

                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[i][1]);
                             ang_Atom_quad.push_back(bonds[j][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][0]){

                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();

                         }
                         if(bonds[j][0] == bonds[k][1]){

                             ang_Atom_quad.push_back(bonds[k][0]);
                             ang_Atom_quad.push_back(bonds[k][1]);
                             ang_Atom_quad.push_back(bonds[j][1]);
                             ang_Atom_quad.push_back(bonds[i][0]);
                             atom_id_dihedral.push_back(ang_Atom_quad);
                             ang_Atom_quad.clear();
                     }
            }



            }

        }
      }
    }

    vector <vector < string> > atom_id_dihedral_2 = atom_id_dihedral;
    vector <vector < string> > atom_id_dihedral_3;
    vector <int> del_index;

    //Dihedral Degenerescency killer


    while(atom_id_dihedral.size() > 0){

        for(int k = 1; k < atom_id_dihedral.size(); k++){


            if(IsEqual(atom_id_dihedral_2[0],atom_id_dihedral[k])){
             del_index.push_back(k);

            }

       }

         if(del_index.size() > 0){
             atom_id_dihedral_3.push_back(atom_id_dihedral[0]);
             atom_id_dihedral_2.erase(atom_id_dihedral_2.begin());
             atom_id_dihedral.erase(atom_id_dihedral.begin());
         for(int i =0; i < del_index.size(); i++){
            atom_id_dihedral.erase(atom_id_dihedral.begin()+del_index[i]-i-1);
             atom_id_dihedral_2.erase(atom_id_dihedral_2.begin()+del_index[i]-i-1);
        }
         }
         else{
             atom_id_dihedral_3.push_back(atom_id_dihedral[0]);
             atom_id_dihedral_2.erase(atom_id_dihedral_2.begin());
             atom_id_dihedral.erase(atom_id_dihedral.begin());
         }

         del_index.clear();


     }

    this->propers_indexial = atom_id_dihedral_3;

}


void aminoacid::PARM(string residue_name){

    string buff1, buff2, buff3, buff4;
    string test_line;
    this->aminoacid_type = residue_name;
    cout << "Residue " <<  residue_name << endl;

    vector < vector < string > > bonds;
    vector < string > bonds_temp;

    string line ;
    ifstream myfile ("/home/jvscunha/workspace/gromacs/share/gromacs/top/amber99sb-ildn.ff/aminoacids.rtp");
    if(myfile){
        string comparative = "[ "+residue_name+" ]";
        while ( line != comparative )
        {
        getline (myfile,line);

    }
        getline (myfile,line);
        getline (myfile,line);

        while ( line != " [ bonds ]" )
        {

        std::istringstream iss(line.c_str());
        iss >> buff1 >> buff2 >> buff3 >> buff4 ;
        this->atom_names.push_back(buff1);
        this->atom_types.push_back(buff2);
        this->atomic_indexes.push_back(this->starting_atomic_index);
        this->starting_atomic_index = this->starting_atomic_index+1;
        this->charges.push_back(stod(buff3));

        getline (myfile,line);

        }

        getline (myfile,line);

        while ( line != " [ impropers ]" )
        {
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2  ;
            //cout << buff1 << " " << buff2 << endl;
            bonds_temp.push_back(buff1);
            bonds_temp.push_back(buff2);
            this->bonds.push_back(bonds_temp);
            bonds_temp.clear();
            getline (myfile,line);


    }


        if(this->C_or_N ==1){
             this->bonds.pop_back();

        }
        if(this->C_or_N==2 || this->C_or_N==3){
         this->bonds = add_previous_backbone(this->bonds);
        }

        this->bonds_conected = this->bonds;

        getline (myfile,line);

       test_line = line;
       vector <string> improper_temp;
        test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());
        while( test_line !=""){
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 ;
            //cout << buff1 << " " << buff2 << " " << buff3 << " " << buff4 << endl;
            improper_temp.push_back(buff1);
            improper_temp.push_back(buff2);
            improper_temp.push_back(buff3);
            improper_temp.push_back(buff4);

            this->impropers.push_back(improper_temp);
            improper_temp.clear();
            getline (myfile,line);
            test_line = line;
            test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());
        }
    }

    this->define_angles(this->bonds);

    this->define_proper_torsions(this->bonds);


    this->translator_name_to_index_bonds();
    this->define_angles_indexial(this->bonds_indexial);
    this->define_proper_torsions_indexial(this->bonds_indexial);


    /*
    for(int i =0; i < this->angles.size() ; i++){
        cout << this->angles[i][0] << " " << this->angles[i][1] << " " << this->angles[i][2] << " " << endl;
    }*/

    if(this->C_or_N==2 || this->C_or_N==3){
        this->bond_angle_cleaner();
    }
}

void aminoacid::SPECIAL_PARM(string residue_name){
    string buff1, buff2, buff3, buff4;
    string test_line;
    this->aminoacid_type = residue_name;
    cout << "Residue " <<  residue_name << endl;

    vector < vector < string > > bonds;
    vector < string > bonds_temp;

    string line ;
    ifstream myfile ("/home/jvscunha/workspace/gromacs/share/gromacs/top/amber99sb-ildn.ff/aminoacids.rtp");
    if(myfile){
        string comparative = "[ "+residue_name+" ]";
        while ( line != comparative )
        {
        getline (myfile,line);

    }
        getline (myfile,line);
        getline (myfile,line);

        while ( line != " [ bonds ]" )
        {

        std::istringstream iss(line.c_str());
        iss >> buff1 >> buff2 >> buff3 >> buff4 ;
        this->atom_names.push_back(buff1);
        this->atom_types.push_back(buff2);
        this->atomic_indexes.push_back(this->starting_atomic_index);
        this->starting_atomic_index = this->starting_atomic_index+1;
        this->charges.push_back(stod(buff3));
        getline (myfile,line);

        }

        getline (myfile,line);
        while ( line != " [ dihedrals ]" )
        {
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2  ;
            bonds_temp.push_back(buff1);
            bonds_temp.push_back(buff2);
            this->bonds.push_back(bonds_temp);
            bonds_temp.clear();
            getline (myfile,line);


    }

        if(this->C_or_N==2 || this->C_or_N==3){
        this->bonds = add_previous_backbone(this->bonds);
        }

      getline (myfile,line);

        while ( line != " [ impropers ]" )
        {
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 ;
            getline (myfile,line);


    }
        getline (myfile,line);



        vector <string> improper_temp;

        test_line = line;
        test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());
        while( test_line !=""){
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 ;
            improper_temp.push_back(buff1);
            improper_temp.push_back(buff2);
            improper_temp.push_back(buff3);
            improper_temp.push_back(buff4);
            this->impropers.push_back(improper_temp);
            getline (myfile,line);
            test_line = line;
            test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());
        }

    }

    this->define_angles(this->bonds);
    this->define_proper_torsions(this->bonds);
    this->translator_name_to_index_bonds();

    this->define_angles_indexial(this->bonds_indexial);

    this->define_proper_torsions_indexial(this->bonds_indexial);


    if(this->C_or_N==2 || this->C_or_N==3){
        this->bond_angle_cleaner();
    }

}

void aminoacid::indexial_impropers( vector < vector < string > > impropers){


}

void aminoacid::bond_angle_cleaner(){
    /*
    cout << "cleaning " << this->bonds.size() << " " <<  this->bonds_indexial.size() << endl;
    for(int i =0; i < this->bonds.size(); i++){
        cout << this->bonds[i][0] << " " << this->bonds[i][1] << endl;
    }*/

    this->bonds.erase(this->bonds.end() - 1);
    this->bonds.erase(this->bonds.end() - 1);
    this->bonds.erase(this->bonds.end() - 1);
    this->bonds.erase(this->bonds.end() - 1);
    this->bonds.erase(this->bonds.end() - 1);


    this->bonds_indexial.erase(this->bonds_indexial.end() - 1);
    this->bonds_indexial.erase(this->bonds_indexial.end() - 1);
    this->bonds_indexial.erase(this->bonds_indexial.end() - 1);
    this->bonds_indexial.erase(this->bonds_indexial.end() - 1);
    this->bonds_indexial.erase(this->bonds_indexial.end() - 1);

    /*
    for(int i =0; i < this->bonds.size(); i++){
        cout << " bonds names " << this->bonds[i][0] << " " << this->bonds[i][1] << endl;
    }
    for(int i =0; i < this->bonds_indexial.size(); i++){
        cout << "bonds indexiais  " <<  this->bonds_indexial[i][0] << " " << this->bonds_indexial[i][1] << endl;
    }


    cout << "cleaning " << this->bonds.size() << " " <<  this->bonds_indexial.size() << endl;
    */

    this->bonds_conected = this->bonds;



    for(int i = 0; i < this->angles.size(); i++){
        for(int k =0 ; k < 3; k++){
            if(this->angles[i][k] == "-R" || this->angles[i][k] == "-N" || this->angles[i][k] == "-HA"){

                this->angles.erase(this->angles.begin()+i);

                i =0 ;
                break;

            }
            if((this->angles[i][0] == "-O" && this->angles[i][2] == "-CA") || (this->angles[i][2] == "-O" && this->angles[i][0] == "-CA")){


                this->angles.erase(this->angles.begin()+i);
                i =0 ;
                break;
            }


        }

    }

}

void aminoacid::PARM_03ws(string residue_name){

    string buff1, buff2, buff3, buff4;
    string test_line;
    this->aminoacid_type = residue_name;
    cout << "Residue " <<  residue_name << endl;

    vector < vector < string > > bonds;
    vector < string > bonds_temp;

    string line ;
    ifstream myfile ("/home/jvscunha/workspace/gromacs/share/gromacs/top/amber99sb-ildn.ff/aminoacids.rtp");
    if(myfile){
        string comparative = "[ "+residue_name+" ]";
        while ( line != comparative )
        {
        getline (myfile,line);

    }
        getline (myfile,line);
        getline (myfile,line);

        while ( line != " [ bonds ]" )
        {

        std::istringstream iss(line.c_str());
        iss >> buff1 >> buff2 >> buff3 >> buff4 ;
        this->atom_names.push_back(buff1);
        this->atom_types.push_back(buff2);
        this->atomic_indexes.push_back(this->starting_atomic_index);
        this->starting_atomic_index = this->starting_atomic_index+1;
        this->charges.push_back(stod(buff3));

        getline (myfile,line);

        }

        getline (myfile,line);

        while ( line != " [ impropers ]" )
        {
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2  ;
           //cout << buff1 << " " << buff2 << endl;
            bonds_temp.push_back(buff1);
            bonds_temp.push_back(buff2);
            this->bonds.push_back(bonds_temp);
            bonds_temp.clear();
            getline (myfile,line);


    }


        if(this->C_or_N ==1){
               this->bonds.pop_back();
        }
        if(this->C_or_N==2 || this->C_or_N==3){
         this->bonds = add_previous_backbone(this->bonds);
        }

        this->bonds_conected = this->bonds;

        getline (myfile,line);

       test_line = line;
       vector <string> improper_temp;
        test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());

       for(int i =0; i< 2; i++){
            std::istringstream iss(line.c_str());
            iss >> buff1 >> buff2 >> buff3 >> buff4 ;
            //cout << buff1 << " " << buff2 << " " << buff3 << " " << buff4 << endl;
            improper_temp.push_back(buff1);
            improper_temp.push_back(buff2);
            improper_temp.push_back(buff3);
            improper_temp.push_back(buff4);
            this->impropers.push_back(improper_temp);
            improper_temp.clear();
            getline (myfile,line);
            test_line = line;
            test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());
        }


       //cout << buff1 << " " << buff2 << " " << buff3 << " " << buff4 << endl;
       improper_temp.push_back(buff1);
       improper_temp.push_back(buff2);
       improper_temp.push_back(buff3);
       improper_temp.push_back(buff4);
       this->impropers.push_back(improper_temp);
       improper_temp.clear();
       getline (myfile,line);
       test_line = line;
       test_line.erase(remove(test_line.begin(), test_line.end(), ' '),test_line.end());

    }

    this->define_angles(this->bonds);

    this->define_proper_torsions(this->bonds);

    this->translator_name_to_index_bonds();
    this->define_angles_indexial(this->bonds_indexial);
    this->define_proper_torsions_indexial(this->bonds_indexial);
    if(this->C_or_N==2 || this->C_or_N==3){
        this->bond_angle_cleaner();
    }
}
