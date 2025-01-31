
#include<iostream>
#include<string>
#include"MP2D.h"
#include<sstream>
#include<fstream>
#include <cstring>
#include <cassert>
#include <stdio.h>
#include <stdlib.h>


using namespace std;

int main(int argc, char ** argv) {
  // This opens the input file
  const int maxlength = 100;
  char filename[maxlength];
  if (!argv[1]) {
    cerr << endl << "MP2D library v1.2   Aug 2024   https://github.com/Chandemonium/MP2D" << endl;
    cerr << "  Computes MP2D and SCS-MP2D dispersion corrections for an MP2 energy." << endl;
    cerr << "  Optionally also computes the MP2D nuclear gradient." << endl;
    cerr << "  Input requires a Cartesian XYZ file and, for SCS-MP2D, some energy data (see below)." << endl;
    cerr << endl << "Syntax:       MP2D <file.xyz> <options>" << endl;
    cerr << endl << "Options:" << endl;
    cerr << "  --scsmp2d   Request SCS-MP2D calculation.  Must also provide --hf, --mp2os, and --mp2ss values." << endl;
    cerr << "  --hf        Provide HF energy for SCS-MP2D calculation (in hartrees)." << endl;
    cerr << "  --mp2os     Provide MP2 opposite-spin correlation energy for SCS-MP2D calculation (in hartrees)." << endl;
    cerr << "  --mp2ss     Provide MP2 same-spin correlation energy for SCS-MP2D calculation (in hartrees)." << endl;
    cerr << endl;
    cerr << "  --gradient  Compute the nuclear gradient of the MP2D dispersion correction (MP2D only)." << endl;
    cerr << endl;
    cerr << "  --param     Specify custom parameter mp2d_parameters.dat and/or scsmp2d_parameters.dat file path" << endl;
    cerr << "  --uchf      Specify custom UCHF C6 coefficient file" << endl;
    cerr << "  --cks       Specify custom CKS C6 coefficient file" << endl;
    cerr << "  --w         Specify double damping short-range width parameter" << endl;
    cerr << "  --rcut      Specify double damping short-range cutoff distance" << endl;
    cerr << "  --TTa1      Specify Tang-Toennies a1 parameter" << endl;
    cerr << "  --TTa2      Specify Tang-Toennies a2 parameter" << endl;
    cerr << "  --s8        Specify s8 scaling coeffient" << endl;
    cerr << "  --cos       Specify MP2 opposite-spin correlation scaling coeffient" << endl;
    cerr << "  --css       Specify MP2 same-spin correlation scaling coeffient" << endl;
 
    
    cerr << endl;
    cerr << "Examples:" << endl;
    cerr << "  MP2D:       MP2D benzene.xyz" << endl;
    cerr << "  SCS-MP2D:   MP2D benzene.xyz --scsmp2d --hf=-230.78162420426395 --mp2ss=-0.2399095404665765 --mp2os=-0.7228154950870085" << endl;
    cerr << endl;
    
    return 1;
    }
    string firstarg = argv[1];
    if (firstarg == "--version") {
        cout << "@mp2d_VERSION@" << endl;
        return 0;
  }
  strcpy(filename, argv[1]);
  
  ifstream infile;
  infile.open(filename);
  assert(infile.is_open());
  //    infile.close();
  
  MP2D::mp2d().Initialize(infile, argc, argv);
    
}
