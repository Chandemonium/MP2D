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
    const int maxlength = 50;
    char filename[maxlength];
    if (!argv[1]) {
        cerr << "No input file provided" << endl;
        return 1;
    }
    strcpy(filename, argv[1]);
    
    ifstream infile;
    infile.open(filename);
    assert(infile.is_open());
//    infile.close();

    MP2D::mp2d().Initialize(infile, argc, argv);

}



