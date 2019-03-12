#include "run_engine.h"

run_engine::run_engine(string protein_gro_file,string outputfile,string dssp_file)
{
    //Protein* receptor = new Protein(protein_gro_file);
    writer* output = new writer(outputfile,protein_gro_file,dssp_file);
    //AMBER_parm_parser* parm = new AMBER_parm_parser();

}

