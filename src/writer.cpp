#include "writer.h"

writer::writer(string outtop,string protein_gro,string dssp_file)
{
    this->savetopol(outtop,protein_gro,dssp_file);
}

void writer::savetopol(string outtop,string protein_gro,string dssp_file){
    Protein* receptor = new Protein(protein_gro,dssp_file);

     FILE * output_topology;
     output_topology = fopen (outtop.c_str(),"w");
     fprintf(output_topology, "[ defaults ]\n");
     fprintf(output_topology, "; nbfunc        comb-rule       gen-pairs       fudgeLJ fudgeQQ\n");
     fprintf(output_topology, "1               2               yes             0.5     0.8333\n");
     fprintf(output_topology, "\n");





     fprintf(output_topology, "#include \"amber99sb-ildn.ff/ffnonbonded.itp\"\n");
     fprintf(output_topology, "\n");
     fprintf(output_topology, "[ moleculetype ]\n");
     fprintf(output_topology, "Protein          3\n");
     fprintf(output_topology, "\n");
     fprintf(output_topology, "[ atoms ]\n");
     int atom_index = 1;
     for(int i =0; i < receptor->aminoacid_parameters.size(); i++ ){
         for(int j =0; j < receptor->aminoacid_parameters[i]->atom_names.size() ; j++ ){
            // cout << atom_index << " " << receptor->aminoacid_parameters[i]->atom_types[j] << " " << i << " " <<  receptor->aminoacid_parameters[i]->aminoacid_type << " " << atom_index << " "
            //      << receptor->aminoacid_parameters[i]->charges[j] << " " << receptor->indexial_masses[atom_index-1] << endl;
             fprintf(output_topology, "%d   %s   %d    %s  %s  %d    %f    %f\n",
                     atom_index,
                     receptor->aminoacid_parameters[i]->atom_types[j].c_str(),
                     i+1,
                     receptor->aminoacid_parameters[i]->aminoacid_type.c_str(),
                     receptor->aminoacid_parameters[i]->atom_names[j].c_str(),
                     atom_index,
                     receptor->aminoacid_parameters[i]->charges[j],
                     receptor->indexial_masses[atom_index-1]);

             atom_index++;
         }

     }
     fprintf(output_topology, "\n");
     fprintf(output_topology, "[ pairs ]\n");
     for(int i =0; i < receptor->pairs_ready.size(); i++ ){
            // cout << atom_index << " " << receptor->aminoacid_parameters[i]->atom_types[j] << " " << i << " " <<  receptor->aminoacid_parameters[i]->aminoacid_type << " " << atom_index << " "
            //      << receptor->aminoacid_parameters[i]->charges[j] << " " << receptor->indexial_masses[atom_index-1] << endl;
             fprintf(output_topology, "%s   %s   1\n",
                     receptor->pairs_ready[i][0].c_str(),
                     receptor->pairs_ready[i][1].c_str());
    }

     fprintf(output_topology, "\n");
     fprintf(output_topology, "[ bonds ]\n");

     for(int i =0; i < receptor->bonds_parameter_list.size(); i++ ){
            // cout << atom_index << " " << receptor->aminoacid_parameters[i]->atom_types[j] << " " << i << " " <<  receptor->aminoacid_parameters[i]->aminoacid_type << " " << atom_index << " "
            //      << receptor->aminoacid_parameters[i]->charges[j] << " " << receptor->indexial_masses[atom_index-1] << endl;
             fprintf(output_topology, "%s   %s   %d    %f    %f\n",
                     receptor->bonds_parameter_list[i].atom_a_index.c_str(),
                     receptor->bonds_parameter_list[i].atom_b_index.c_str(),
                     receptor->bonds_parameter_list[i].function,
                     receptor->bonds_parameter_list[i].b0,
                     receptor->bonds_parameter_list[i].kb);



     }
     fprintf(output_topology, "\n");
     fprintf(output_topology, "[ angles ]\n");

     for(int i =0; i < receptor->angles_parameter_list.size(); i++ ){
            // cout << atom_index << " " << receptor->aminoacid_parameters[i]->atom_types[j] << " " << i << " " <<  receptor->aminoacid_parameters[i]->aminoacid_type << " " << atom_index << " "
            //      << receptor->aminoacid_parameters[i]->charges[j] << " " << receptor->indexial_masses[atom_index-1] << endl;
             fprintf(output_topology, "%s   %s  %s  %d    %f    %f\n",
                     receptor->angles_parameter_list[i].atom_a_index.c_str(),
                     receptor->angles_parameter_list[i].atom_b_index.c_str(),
                     receptor->angles_parameter_list[i].atom_c_index.c_str(),
                     receptor->angles_parameter_list[i].function,
                     receptor->angles_parameter_list[i].th0,
                     receptor->angles_parameter_list[i].cth);



     }

     fprintf(output_topology, "\n");
     fprintf(output_topology, "[ dihedrals ]\n");
     for(int i =0; i < receptor->proper_parameter_list.size(); i++ ){
            // cout << atom_index << " " << receptor->aminoacid_parameters[i]->atom_types[j] << " " << i << " " <<  receptor->aminoacid_parameters[i]->aminoacid_type << " " << atom_index << " "
            //      << receptor->aminoacid_parameters[i]->charges[j] << " " << receptor->indexial_masses[atom_index-1] << endl;
             fprintf(output_topology, "%s   %s  %s %s %d    %f    %f %f \n",
                     receptor->proper_parameter_list[i].atom_a_index.c_str(),
                     receptor->proper_parameter_list[i].atom_b_index.c_str(),
                     receptor->proper_parameter_list[i].atom_c_index.c_str(),
                     receptor->proper_parameter_list[i].atom_d_index.c_str(),
                     receptor->proper_parameter_list[i].function,
                     receptor->proper_parameter_list[i].phase,
                     receptor->proper_parameter_list[i].kd,
                     receptor->proper_parameter_list[i].pn);



     }

     for(int i =0; i < receptor->improper_parameter_list.size(); i++ ){
            // cout << atom_index << " " << receptor->aminoacid_parameters[i]->atom_types[j] << " " << i << " " <<  receptor->aminoacid_parameters[i]->aminoacid_type << " " << atom_index << " "
            //      << receptor->aminoacid_parameters[i]->charges[j] << " " << receptor->indexial_masses[atom_index-1] << endl;
             fprintf(output_topology, "%s   %s  %s %s %d    %f    %f %f \n",
                     receptor->improper_parameter_list[i].atom_a_index.c_str(),
                     receptor->improper_parameter_list[i].atom_b_index.c_str(),
                     receptor->improper_parameter_list[i].atom_c_index.c_str(),
                     receptor->improper_parameter_list[i].atom_d_index.c_str(),
                     receptor->improper_parameter_list[i].function,
                     receptor->improper_parameter_list[i].phase,
                     receptor->improper_parameter_list[i].kd,
                     receptor->improper_parameter_list[i].pn);



     }
     fprintf(output_topology, "\n");
     fprintf(output_topology, "#ifdef POSRES\n#include \"posre.itp\"\n#endif\n\n#include \"amber99sb-ildn.ff/tip3p.itp\"\n#ifdef POSRES_WATER\n[ position_restraints ]\n1    1       1000       1000       1000\n#endif\n#include \"amber99sb-ildn.ff/ions.itp\"\n[ system ]\nProtein\n[ molecules ]\nProtein    1\n");





}
