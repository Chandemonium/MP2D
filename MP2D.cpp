#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include"MP2D.h"
#include<iterator>
#include<unistd.h>
#include<string.h>
#include<iomanip> 
#include<algorithm> 
#include<functional> 

//using namespace std::placeholders; 

MP2D::MP2D() : Natoms_per_monomer(NULL), AtomicSymbols(NULL), 
AtomicNumbers(NULL), xyz(NULL), symbols(NULL), RcovAB(NULL) // r0AB(NULL), RcovAB(NULL)

{}

void MP2D::Initialize(ifstream &infile, int argc, char *argv[]) {

  AngToBohr = 1.8897259886;
  HartreeToKcal = 627.5095; 

  // First read any parameters from the command line
  GetUserParameters(argc,argv);

  // Now read the default empirical parameters.  If the user specified some parameters in the command line,
  // these functionals will replace the default values with the user's values.
 if (SCSMP2D_U == true) {
   GetSCSMP2DParameters(infile);
 }
 else{
   GetMP2DParameters(infile);
 }
 
  Ntot = BeforeGetTotalNumberOfAtoms(infile);
  GetCoordinates(infile);
  *GetAtomicNumber();
  *GetMultipoleExpectationValue();
  GetCovalentRadii();
  GetC6Coefficients();
  CoordinationNumber();
  CoordinationNumberGradient();
  ComputeEnergyandPrintResults();


}

//Destructor
MP2D::~MP2D() {

    delete [] Natoms_per_monomer;
    delete [] xyz;
}    


// This function allows the user to specify parameter values for the calculations, and paths for the coefficients.
void MP2D::GetUserParameters(int argc, char* argv[]) {

  // This loop matches the user specified inputs by their first several characters, and sets a few variables based on the matches that are made.
  for (int i = 0; i<argc; i++) {
    string str = argv[i]; 
    // With this keyword you can set an alternative parameter file instead of using mp2d_parameter.dat
    if (str.substr(0,7) == "--param") {
      Param_U = true;
      string PARAM_U = str.substr(8);
      stringstream PARAM(PARAM_U);
      PARAM >> Param;
      cout << Param << endl;
    }
    // With this keyword you can set a path to where the UCHF coefficients are located
    if (str.substr(0,6) == "--uchf") {
      UCHF_U = true;
      string Uchf_U = str.substr(7);
      stringstream UCHF(Uchf_U);
      UCHF >> UCHF_Path;
      cout << UCHF_Path << endl;
    }
    // With this path you can set a path to where the CKS coefficients are located
    else if (str.substr(0,5) == "--cks") {
      // cout << "yay again";
      CKS_U = true;
      string Cks_U = str.substr(6);
      stringstream CKS(Cks_U);
      CKS >> CKS_Path;
      //cout << CKS_Path << endl;
    }
    // This keyword lets you specify an alternative width for the double damping function
    else if (str.substr(0,3) == "--w") {
      w_U = true;
      string W_U = str.substr(4);
      stringstream WIDTH(W_U);
      Width = 0.0;
      WIDTH >> Width;
      //cout << Width << endl;
    }
    // This keyword lets you specify an alternative cutoff radius for the double damping function
    else if (str.substr(0,6) == "--rcut") {
      rcut_U = true;
      string Rcut_U = str.substr(7);
      stringstream RCUT(Rcut_U);
      Rcut = 0.0;
      RCUT >> Rcut;
      cout << Rcut << endl;
    }
    // This keyword lets you specify an alternative a1 parameter for the Tang-Toennies damping function
    else if (str.substr(0,7) == "--TT_a1") {
      TT_a1_U = true;
      string TTa1_U = str.substr(8);
      stringstream TTA1(TTa1_U);
      TT_a1 = 0.0;
      TTA1 >> TT_a1;
    }
    // This keyword lets you specify an altrernative a2 parameter for the Tang-Toennies damping function
    else if (str.substr(0,7) == "--TT_a2") {
      TT_a2_U = true;
      string TTa2_U = str.substr(8);
      stringstream TTA2(TTa2_U);
      TT_a2 = 0.0;
      TTA2 >> TT_a2;
      //cout << TT_a2 << endl;
    }
    // This keyword lets you specify an alternative s8 parameter which compensates for the lack of higher order terms in the energy expansion
    else if (str.substr(0,4) == "--s8") {
      s8_U = true;
      string S8_U = str.substr(5);
      stringstream S8(S8_U);
      s8 = 0.0;
      S8 >> s8;
      //cout << s8 << endl;
    }
    // This keyword lets you specify an alternative Css same-spin correlation scaling parameter
    else if (str.substr(0,5) == "--css") {
      Css_U = true;
      string tmp_U = str.substr(6);
      stringstream tmp(tmp_U);
      userCss = 0.0;
      tmp >> userCss;
      //cout << "User requests Css = " << userCss  << endl;
    }
    // This keyword lets you specify an alternative Cos opposite-spin correlation scaling parameter
    else if (str.substr(0,5) == "--cos") {
      Cos_U = true;
      string tmp_U = str.substr(6);
      stringstream tmp(tmp_U);
      userCos = 0.0;
      tmp >> userCos;
      //cout << "User requests Css = " << userCss  << endl;
    }
    
    // This keyword lets you specify a Hartree-Fock energy for SCS-MP2D
    else if (str.substr(0,4) == "--hf") {
      hf_U = true;
      string hf_U = str.substr(5);
      stringstream HF(hf_U);
      hf = 0.0;
      HF >> hf;
      //cout << s8 << endl;
    }
    // This keyword lets you specify a mp2 os energy for SCS-MP2D
    else if (str.substr(0,7) == "--mp2os") {
      mp2os_U = true;
      string mp2os_U = str.substr(8);
      stringstream MP2OS(mp2os_U);
      mp2_os = 0.0;
      MP2OS >> mp2_os;
      //cout << s8 << endl;
    }
    // This keyword lets you specify a mp2 ss energy for SCS-MP2D
    else if (str.substr(0,7) == "--mp2ss") {
      mp2ss_U = true;
      string mp2ss_U = str.substr(8);
      stringstream MP2SS(mp2ss_U);
      mp2_ss = 0.0;
      MP2SS >> mp2_ss;
      //cout << s8 << endl;
    }
    // Specify that this is a SCS-MP2D calculation
    else if (str.substr(0,9) == "--SCSMP2D" || str.substr(0,9) == "--scsmp2d") {
      SCSMP2D_U =true;
    }
    // This parameter lets you specify whether or not you would like to print out the gradient when you run the calculation
    else if (str.substr(0,10) == "--gradient") {
      Gradient_U =true;
    }
    
    else {} 
    str.clear();
  }


  // Error handling for keywords

  // 1) Only allow user-specified Cos/Css settings for SCS-MP2D
  if ( (Cos_U || Css_U ) && !SCSMP2D_U) {
    cerr << "ERROR: To use custom Cos / Css spin-component scaling coefficients, you must run an SCS-MP2D calculation." << endl << "To do so, you must specify the --SCSMP2D flags and give values for the --hf, --mp2os, and --mp2ss flags." << endl;
    exit(1);
  }

  // 2) Gradients only for MP2D, not SCS-MP2D
  if (SCSMP2D_U && Gradient_U) {
    cerr << "ERROR: Gradients are only implemented for MP2D, not for SCS-MP2D." << endl;
    exit(1);
  }

  // 3) General case where a keyword doesn't match a a known keyword.
  if (argv[2] && Param_U==false && UCHF_U == false && CKS_U == false && w_U == false && rcut_U == false && TT_a1_U == false && TT_a2_U == false && s8_U == false && Cos_U == false && Css_U == false && hf_U == false && mp2os_U == false && mp2ss_U == false &&  SCSMP2D_U == false  && Gradient_U == false    ) {
    cerr << "Calculation terminated, the provided input does not match any known keyword. Please refer to README.txt" << endl;
    exit(1);
  }
}

// Reads empirical MP2D parameters from file
void MP2D::GetSCSMP2DParameters(ifstream& infile) {

  ifstream inFile;  
  string psiP;
  
  // Set the path to the parameter file. The default is mp2d_parameters.dat
  if (Param_U == true) {
    psiP = Param;
  }
  // Check if MP2D_PARAM_PATH is set. If it is then set psiP. If not it 
  // will throw an error down the line.
  else if (getenv("MP2D_PARAM_PATH")!=NULL) {
    psiP = getenv("MP2D_PARAM_PATH");
  }
  else {}
  if (psiP.empty()) {
    cerr << "Calculation terminated, please add an MP2D_PARAM_PATH environment varaible. Refer to README.txt" << endl;
    exit (1);
  }
  
  string file2 = "/scsmp2d_parameters.dat";
  string path = psiP +file2;
  
  inFile.open(path.c_str());
  string line;
  
  string aa;
  string bb;
  double cc;
  
  if (!inFile) {
    cerr << "Unable to open file scsmp2d_parameters.dat" << endl;
    exit (1);
    }
  
  
  // loop through the parameter file to grab the numerical parameters for the damping functions, and the energy expression. Also, get the paths to the coefficient files.
  while(getline(inFile,line)) {    
    istringstream iss(line);
    iss >> aa >> bb >> cc;
    
    if (aa =="a_one") {a_one = cc;}
    if (aa == "a_two") {a_two = cc;}
    if (aa == "rcut") {rcut = cc;}
    if (aa == "width") {width = cc;}
    if (aa == "s_8") {s_8 = cc;}
    if (aa == "Cos") {Cos = cc;}
    if (aa == "Css") {Css = cc;}
    if (aa == "CKS_Path:") {GrimmePath = bb;}
    if (aa == "UCHF_Path:") {UCHFPath = bb;}   
    
  }
  inFile.close();
    
  // If appropriate, replace the parameter file values with any specified command-line parameter values.
  if (rcut_U == true) {rcut = Rcut;}
  if (w_U == true) {width = Width;}
  if (TT_a1_U == true) {a_one = TT_a1;}
  if (TT_a2_U == true) {a_two = TT_a2;}
  if (s8_U == true) {s_8 = s8;}
  if (Cos_U == true) {Cos = userCos;}
  if (Css_U == true) {Css = userCss;}

  
  // Set HF, same-spin MP2 correlation, and opposite-spin MP2 correlation energies
  HFT = hf;
  MP2_os = mp2_os;
  MP2_ss = mp2_ss;

  if ( HFT == 0.0 || MP2_os == 0.0 || MP2_ss == 0.0 ) {
    cerr << "ERROR: The HF and MP2 spin component energies must be supplied for SCS-MP2D." << endl << "Use the --hf, --mp2os, and --mp2ss flags to specify the values (in hartrees)." << endl;
    exit (1);
  }
  
}

// Reads empirical MP2D parameters from file
void MP2D::GetMP2DParameters(ifstream& infile) {

  ifstream inFile;  
  string psiP;
  
  // Set the path to the parameter file. The default is mp2d_parameters.dat
  if (Param_U == true) {
    psiP = Param;
  }
  // Check if MP2D_PARAM_PATH is set. If it is then set psiP. If not it 
  // will throw an error down the line.
  else if (getenv("MP2D_PARAM_PATH")!=NULL) {
    psiP = getenv("MP2D_PARAM_PATH");
  }
  else {}
  if (psiP.empty()) {
    cerr << "Calculation terminated, please add an MP2D_PARAM_PATH environment varaible. Refer to README.txt" << endl;
    exit (1);
  }
  
  string file2 = "/mp2d_parameters.dat";
  string path = psiP +file2;
  
  inFile.open(path.c_str());
  string line;
  
  string aa;
  string bb;
  double cc;
  
  if (!inFile) {
    cerr << "Unable to open file mp2d_parameters.dat" << endl;
    exit (1);
    }
  
  
  // loop through the parameter file to grab the numerical parameters for the damping functions, and the energy expression. Also, get the paths to the coefficient files.
  while(getline(inFile,line)) {    
    istringstream iss(line);
    iss >> aa >> bb >> cc;
    
    if (aa =="a_one") {a_one = cc;}
    if (aa == "a_two") {a_two = cc;}
    if (aa == "rcut") {rcut = cc;}
    if (aa == "width") {width = cc;}
    if (aa == "s_8") {s_8 = cc;}
    if (aa == "Cos") {Cos = cc;}
    if (aa == "Css") {Css = cc;}
    if (aa == "CKS_Path:") {GrimmePath = bb;}
    if (aa == "UCHF_Path:") {UCHFPath = bb;}   
    
  }
  inFile.close();
    
    // If appropriate, replace the parameter file values with any specified command-line parameter values.
  if (rcut_U == true) {rcut = Rcut;}
  if (w_U == true) {width = Width;}
  if (TT_a1_U == true) {a_one = TT_a1;}
  if (TT_a2_U == true) {a_two = TT_a2;}
  if (s8_U == true) {s_8 = s8;}

}




// Read the total number of atoms from the first line of the XYZ input file
int MP2D::BeforeGetTotalNumberOfAtoms(ifstream &infile) {

  Rewind(infile);
  infile >> NumberOfAtoms;
  
  //printf("XYZ file contains %d atoms\n",NumberOfAtoms);
  return NumberOfAtoms;
}

// Read atomic symbols and Cartesian coordinates
void MP2D::GetCoordinates(ifstream& infile) {
  string line;
  
  // Stores the atomic symbols    
  symbols = new string[Ntot];
  
  // Stores the coordinates of the system
  XYZ.resize(3*Ntot);
  
  Rewind(infile);
  
  
  int atom = 0;
  
  // Skip the line for the number of atoms and the line for comments
  getline(infile,line);
  getline(infile,line);
  
  // starting from the third line, loop through and store the coordinates
  
  for (int i=0; i< NumberOfAtoms; i++) {     
    infile >> symbols[atom];
    
    // Convert coordinates to Bohr as they are pulled in
    infile >> XYZ[3*atom];
    XYZ[3*atom] = XYZ[3*atom]*AngToBohr;
    infile >> XYZ[3*atom+1];
    XYZ[3*atom+1] = XYZ[3*atom+1]*AngToBohr;
    infile >> XYZ[3*atom+2];
    XYZ[3*atom+2] = XYZ[3*atom+2]*AngToBohr;
    
    atom++;
  } 
}

int *MP2D::GetAtomicNumber() {
  
  // atomic numbers of elements that we have UCHF coefficients for
  // the index of a given element corresponds to its atomic number
  string AtomicSymbols[36] = {
    "XX",
    "H","XX",
    "XX","XX","B","C","N","O","F","Ne", // period 2
    "XX","XX","XX","XX","P","S","Cl","Ar", // period 3
    "XX","XX","XX","XX","XX","XX","XX","XX","XX","XX","XX","XX", // period 4
    "XX","XX","XX","XX","Br"}; // period 4
  
  int i, j, at_num = 0;
  
  atomic_number = new int[Ntot];
  // The array atomic_number is filled with the atomic numbers of each 
  // atom that appears in the input.
  int counter = 0;
  for (i=0;i<Ntot;i++) {
    for (j=0;j<36;j++) {
      if (symbols[i]==AtomicSymbols[j]) {
	at_num = j;
	atomic_number[i] = at_num;
	counter++;
      }
    }
    
  }
  if (counter!=Ntot) {
    cerr << "Calculation terminated, the input contains atoms that are not currently supported by MP2D. Please refer to README.txt." << endl;
    exit(1);
  }
  
  // Here the returned argument is a pointer to the array atomic_number[]
  return atomic_number;
}

double *MP2D::GetMultipoleExpectationValue() {
  
  double MultipoleExpectationValues[36] = { 0.0, 8.0589, 0.0, 0.0, 0.0, 11.8799, 7.8715, 5.5588, 4.7566, 3.8025
,3.1036, 0.0, 0.0, 0.0, 0.0, 9.5361, 8.1652, 6.7463, 5.6004, 0.0 ,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 7.1251};

  int i, j;
  r2_r4 = new double[Ntot]; 
  
  for (int i=0;i<Ntot;i++) {
        r2_r4[i] = MultipoleExpectationValues[atomic_number[i]];
  }
  return r2_r4;
} 

void MP2D::GetCovalentRadii() {
  // Covalent radii for H,B,C,N,O,F,Ne,P,S,Cl,Ar, and Br
  // Each covalent radii is placed at the index corresponding to the element's
  // atomic number
  double CovalentRadii[36] =  {
    0.0,
    0.32,0.0,
    0.0,0.0,0.77,0.75,0.71,0.63,0.64,0.67, // period 2
    0.0,0.0,0.0,0.0,1.1,1.02,0.99,0.96, // period 3
    0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0, // period 4
    0.0,0.0,0.0,0.0,1.14}; // period 4
  
  double Rcov[Ntot];
  for (int i=0; i<Ntot; i++ ) {
    Rcov[i] = CovalentRadii[atomic_number[i]];
  }
  
  RcovAB = new double*[Ntot];
  for(int i=0; i < Ntot; i++) {
    RcovAB[i] = new double[Ntot];
  }
  for (int i=0; i<Ntot; i++) {
    for (int j= 0; j<Ntot; j++) {
      if (i==j) {
	RcovAB[i][j] = 0.0;
      }
      else {
	RcovAB[i][j] = (Rcov[i] + Rcov[j])*AngToBohr;
      }
    }
  } 
}

double MP2D::ComputeDistance(double a, double b, double c, double d, double e, double f) {

  double R_AB = sqrt( pow(a - b, 2) + pow(c - d, 2) + pow(e - f, 2));
  return R_AB;
  
}

void MP2D::GetC6Coefficients() {
  
  C6_counter = 0;
  
  ifstream inFile2;
  string psiP;
  
  if (CKS_U == true) {
    psiP = CKS_Path;
  }
  else {
    psiP = GrimmePath;
  }
  if (psiP.empty()) {
    printf("Please set MP2D_PARAM_PATH shell environment variable to point directory containing mp2d_parameters.dat file\n");    
  }
  
  string path = psiP;
  CKSPATH = path;
  
  inFile2.open(path.c_str());
  string line2;
  if (!inFile2) {
    cerr << "unable to open file: " << path << endl;
    exit (1);
  }
  
  while (getline(inFile2, line2)) {
    C6_counter++;
  }
  inFile2.close();
  
  // define matrix here with dimensions index x 5
  
  matrix_Grimme.resize(C6_counter);
  
  for (int i=0; i<C6_counter; i++) {
    matrix_Grimme[i].resize(5);
  }
  
  inFile2.open(path.c_str());
  string Line2;
  if (!inFile2) {
    cerr << "unable to open file: " << path << endl;   
    exit (1);
  }
  
  int counter_2 =0;
  
  while (getline(inFile2, line2)) {
    istringstream ss(line2);
    
    double atomA, atomB, CN_A, CN_B, C6CKS;
    ss >> atomA >> atomB >> CN_A >> CN_B >> C6CKS;
    
    matrix_Grimme[counter_2][0] = atomA;
    matrix_Grimme[counter_2][1] = atomB;
    matrix_Grimme[counter_2][2] = CN_A;
    matrix_Grimme[counter_2][3] = CN_B;
    matrix_Grimme[counter_2][4] = C6CKS;
    counter_2++;
  }
  inFile2.close();
  
  // Where the UCHF coefficients are stored
  UCHF_counter = 0;
  
  ifstream inFile;
  string psiP2;
  if (UCHF_U == true) {
    psiP2 = UCHF_Path;
  }
  else {
    psiP2 = UCHFPath;
  }
  
  if (psiP2.empty()) {
    printf("Please set MP2D_PARAM_PATH shell environment variable to point directory containing mp2d_parameters.dat file\n");    
  }
  
  
  string path2 = psiP2;
  UCHFPATH = path2;
  
  inFile.open(path2.c_str());
  string line;
  if (!inFile) {
    cerr << "Unable to open file: " << path2 << endl;
    exit (1);
  }
  
  while (getline(inFile, line)) {
    UCHF_counter++;
  }
  inFile.close();
  
  matrix_UCHF.resize(UCHF_counter);
  for (int i=0; i<UCHF_counter; i++) {
    matrix_UCHF[i].resize(5);
  }
  
  inFile.open(path2.c_str());
  if (!inFile) {
    cerr << "unable to open file: " << path2 << endl;
    exit (1);
  }
  
  counter_2 =0;
  
  while (getline(inFile, line)) {
    istringstream ss(line);
    
    double atomA, atomB, CN_A, CN_B, UCHFCKS;
    ss >> atomA >> atomB >> CN_A >> CN_B >> UCHFCKS;
    
    matrix_UCHF[counter_2][0] = atomA;
    matrix_UCHF[counter_2][1] = atomB;
    matrix_UCHF[counter_2][2] = CN_A;
    matrix_UCHF[counter_2][3] = CN_B;
    matrix_UCHF[counter_2][4] = UCHFCKS;
    counter_2++;
  }
  inFile.close();
    
}

double MP2D::Factorial(int n) {
  if ( n>1 ) {
    return n*Factorial(n-1);
  }
  else { return 1; }
}

double MP2D::DoubleDamping(double R0ab, double R, int atnumA, int atnumB ) {
  
  double modR = 0.0;
    
  if (R > R0ab*(rcut + width/2.0)) {
    modR = R;
  }
  
  else if (R < R0ab*(rcut - width/2.0)) {
    modR = rcut*R0ab;
  }
  
  else {
    double r_prime = R0ab*rcut;
    double w_prime = R0ab*width;
    double x = (R - (r_prime - w_prime/2.0))/w_prime;
    double interpolateR = (-2.5*pow(x,8) + 10*pow(x,7) - 14*pow(x,6) + 7*pow(x,5))*w_prime;
    modR = rcut*R0ab + interpolateR;
  }

  return modR;
}

double MP2D::DoubleDampingDeriv(double R0ab, double R, int atnumA, int atnumB ) {

  double modR_Deriv = 0.0;
  
  if (R >= R0ab*(rcut + width/2.0)) {
    modR_Deriv = 1.0;
  }
  
  else if (R <= R0ab*(rcut - width/2.0)) {
    modR_Deriv = 0.0;
  }
  
  else {
    double r_prime = R0ab*rcut;
    double w_prime = R0ab*width;
    double x = (R - (r_prime - w_prime/2.0))/w_prime;
    
    double interpolateR_Deriv = (-20*pow(x,7) + 70*pow(x,6) - 84*pow(x,5) + 35*pow(x,4));
    modR_Deriv = interpolateR_Deriv;
  }
  
  return modR_Deriv;
}

double MP2D::TangToennies(int n, double r0ab, double Rab) {

    double damp_sum = 0.0;
 
    double S = (r0ab*(a_one/pow(AngToBohr,2)) + a_two/AngToBohr)*Rab;

    for (int i = 0; i < n+1; i++) {

        damp_sum += pow(S,i)/Factorial(i);

    }

    return 1 - exp(-1*S)*damp_sum;
    
}

double MP2D::TangToenniesDamp(int n, double r0ab, double Rab) {

  double damp_sum = 0.0;
  double damp_sum_deriv = 0.0;
  
  double S = (r0ab*(a_one/pow(AngToBohr,2)) + a_two/AngToBohr)*Rab;
  double S_prime = (r0ab*(a_one/pow(AngToBohr,2)) + a_two/AngToBohr);
  
  for (int i = 0; i < n+1; i++) {
    damp_sum += pow(S,i)/Factorial(i);
    damp_sum_deriv += ((pow(S_prime,i)*pow(Rab,i-1))/Factorial(i))*i;
  }

  double ttdamp = (S_prime*exp(-1*S)*damp_sum - exp(-1*S)*damp_sum_deriv)*1.8897;
  return ttdamp;
}

void MP2D::CoordinationNumber() {

  double function = 0.0;
  double CN;

  Coordination_Number.resize(Ntot);

  for (int i=0; i < NumberOfAtoms; i++) {
    CN = 0.0;
    for (int j=0; j < NumberOfAtoms; j++ ) {
      double R_AB = ComputeDistance(XYZ[3*i], XYZ[3*j], XYZ[3*i +1], XYZ[3*j +1], XYZ[3*i +2], XYZ[3*j + 2]);
      if (i!=j) {
	if (R_AB <= 0.95*RcovAB[i][j]) {
	  function = 1.0;
	}
	else if (R_AB >= 1.75*RcovAB[i][j]) {
	  function = 0.0;
	}
	else {
	  double numerator = R_AB - 0.95*RcovAB[i][j];
	  double xprime = numerator/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
	  function = 1.0 - (-20*pow(xprime,7) + 70*pow(xprime, 6) - 84*pow(xprime,5) + 35*pow(xprime,4));
	}
      }
      
      else {
	function = 0.0;
      }
      
      CN = CN + function;
    }
    Coordination_Number[i] = CN;  
  }
}

void MP2D::CoordinationNumberGradient() {
  
  double function_gradient = 0.0;
  
  CNgrad.resize(Ntot);
  dC6_CKSGrad.resize(Ntot);
  dC6_UCHFGrad.resize(Ntot);
  for (int i=0; i < Ntot; i++ ) {
    CNgrad[i].resize(3*Ntot);
    dC6_CKSGrad.resize(3*Ntot);
    dC6_UCHFGrad.resize(3*Ntot);
  }
  
  // Initialize CNgrad with 0.0 elements
  for (int i=0; i < NumberOfAtoms; i++) {
    for (int j=0; j < NumberOfAtoms; j++ ) {
      CNgrad[i][j] = 0.0;        
    }
  }

  for (int i=0; i < NumberOfAtoms; i++) {
    for (int j=0; j < NumberOfAtoms; j++ ) {
      double R_AB = ComputeDistance(XYZ[3*i], XYZ[3*j], XYZ[3*i +1], XYZ[3*j +1], XYZ[3*i +2], XYZ[3*j + 2]);
      if (i!=j) {
	if (R_AB <= 0.95*RcovAB[i][j]) { 
	  function_gradient = 0.0;
	}
	else if ( R_AB >= 1.75*RcovAB[i][j]) {
	  function_gradient = 0.0;
	}
	else {
	  double x = (R_AB - 0.95*RcovAB[i][j])/(1.75*RcovAB[i][j] - 0.95*RcovAB[i][j]);
	  double x_grad = 1.25/RcovAB[i][j];
	  function_gradient = (140*pow(x,6) - 420*pow(x,5) + 420*pow(x,4) - 140*pow(x,3))*x_grad;
	}
	
	double dRdxi = function_gradient*(XYZ[3*i]-XYZ[3*j])/R_AB;
	double dRdyi = function_gradient*(XYZ[3*i +1]-XYZ[3*j +1])/R_AB;
	double dRdzi = function_gradient*(XYZ[3*i +2]-XYZ[3*j +2])/R_AB;
	double dRdxj = -dRdxi;
	double dRdyj = -dRdyi;
	double dRdzj = -dRdzi;
	
	CNgrad[i][3*i] += dRdxi;
	CNgrad[i][3*i+1] += dRdyi;
	CNgrad[i][3*i+2] += dRdzi;
	CNgrad[i][3*j] += dRdxj;
	CNgrad[i][3*j+1] += dRdyj;
	CNgrad[i][3*j+2] += dRdzj;
      }
      
    }
  }  
}

double MP2D::ComputeC6_UCHF(int atnumA, int atnumB, double CN_A, double CN_B) {

  double C6_UCHF = 0.0;
  double sum_L_ij=0.0;
  double sum_numerator=0.0;
  double sum_C6ref_Lij =0.0;
  
  for( int w=0; w<UCHF_counter-1; w++) {
    double atA, atB, cnA, cnB, c6uchf;
    
    atA  = matrix_UCHF[w][0];
    atB  = matrix_UCHF[w][1];
    
    cnA  = matrix_UCHF[w][2];
    cnB  = matrix_UCHF[w][3];
    c6uchf  = matrix_UCHF[w][4]; 
    
    if ((atA == atnumA && atB == atnumB) || (atA==atnumB && atB==atnumA)) {
      if (atnumA==atnumB && cnA != cnB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB)))) + exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
      } 
      
      else if (atA ==atnumA && atB ==atnumB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))));
      }
      else if (atA==atnumB && atB ==atnumA) {
	L_ij = exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
      }
      else {}
      sum_L_ij += L_ij;
      
      sum_C6ref_Lij += L_ij*c6uchf;
      
      double numerator = L_ij*c6uchf;
      sum_numerator += numerator;
      
      C6_UCHF = sum_numerator/sum_L_ij;
    }
    
  }
  return C6_UCHF;
} 

double MP2D::ComputeC6_CKS(int atnumA, int atnumB, double CN_A, double CN_B) {

  double C6_CKS = 0.0;
  double sum_L_ij=0.0;
  double sum_numerator=0.0;
  double sum_C6ref_Lij =0.0;
  
  for( int w=0; w<C6_counter-1; w++) {
    double atA, atB, cnA, cnB, c6uchf;
    
    atA  = matrix_Grimme[w][0];
    atB  = matrix_Grimme[w][1];
    
    cnA  = matrix_Grimme[w][2];
    cnB  = matrix_Grimme[w][3];
    c6uchf  = matrix_Grimme[w][4]; 

    if ((atA == atnumA && atB == atnumB) || (atA==atnumB && atB==atnumA)) {
      // unless the coordination numbers are equal, the contribution needs to be counted twice.
      if (atnumA==atnumB && cnA != cnB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB)))) + exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
      } 
      else if (atA ==atnumA && atB ==atnumB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))));
      }
      else if (atA==atnumB && atB ==atnumA) {
	L_ij = exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
      }
      else {}
      
      sum_L_ij += L_ij;
      sum_C6ref_Lij += L_ij*c6uchf;
      double numerator = L_ij*c6uchf;
      sum_numerator += numerator;
      C6_CKS = sum_numerator/sum_L_ij;
    }
  }
  return C6_CKS;
} 


valarray<double> MP2D::ComputeC6_UCHF_Gradient(int atnumA, int atnumB, double CN_A, double CN_B, valarray<double> vectA, valarray<double> vectB) {
  double C6_UCHF = 0.0;
  double C6_UCHF_Deriv = 0.0;
  double sum_L_ij=0.0;
  double sum_numerator=0.0;
  double sum_C6ref_Lij =0.0;
  
  valarray<double> L_ij_deriv;
  L_ij_deriv.resize(3*Ntot);
  
  valarray<double> sum_L_ij_deriv;
  sum_L_ij_deriv.resize(3*Ntot);
  
  valarray<double> sum_C6ref_Lij_deriv;
  sum_C6ref_Lij_deriv.resize(3*Ntot);
  valarray<double> C6_UCHFGrad;
  C6_UCHFGrad.resize(3*Ntot);
  
  for( int w=0; w<UCHF_counter-1; w++) {
    double atA, atB, cnA, cnB, c6uchf;
    
    atA  = matrix_UCHF[w][0];
    atB  = matrix_UCHF[w][1];
    
    cnA  = matrix_UCHF[w][2];
    cnB  = matrix_UCHF[w][3];
    c6uchf  = matrix_UCHF[w][4]; 
    
    if ((atA == atnumA && atB == atnumB) || (atA==atnumB && atB==atnumA)) {
      if (atnumA==atnumB && cnA != cnB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB)))) + exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
	L_ij_deriv = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))))*(-8*(((CN_A-cnA)*vectA)+((CN_B-cnB)*vectB))) + exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))))*(-8*(((CN_A-cnB)*vectA)+((CN_B-cnA)*vectB)));
      } 
      
      else if (atA ==atnumA && atB ==atnumB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))));
	L_ij_deriv = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))))*(-8*(((CN_A-cnA)*vectA)+((CN_B-cnB)*vectB)));
      }
      else if (atA==atnumB && atB ==atnumA) {
	L_ij = exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
	L_ij_deriv = exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))))*(-8*(((CN_A-cnB)*vectA)+((CN_B-cnA)*vectB)));
      }
      else {}

      sum_L_ij += L_ij;
      sum_C6ref_Lij += L_ij*c6uchf;
      sum_L_ij_deriv += L_ij_deriv;
      sum_C6ref_Lij_deriv += L_ij_deriv*c6uchf;
      
      //for (int i=0; i<Ntot; i++ ) {
      //
      //}
      
      double numerator = L_ij*c6uchf;
      sum_numerator += numerator;
      
      C6_UCHF = sum_numerator/sum_L_ij;
    }
  }
  
  C6_UCHFGrad = sum_C6ref_Lij_deriv/sum_L_ij - sum_C6ref_Lij*sum_L_ij_deriv/pow(sum_L_ij,2);
  
  return C6_UCHFGrad;
} 

valarray<double> MP2D::ComputeC6_CKS_Gradient(int atnumA, int atnumB, double CN_A, double CN_B, valarray<double> vectA, valarray<double> vectB) {

  double C6_UCHF = 0.0;
  double C6_UCHF_Deriv = 0.0;
  double sum_L_ij=0.0;
  double sum_numerator=0.0;
  double sum_C6ref_Lij =0.0;
  
  valarray<double> L_ij_deriv;
  L_ij_deriv.resize(3*Ntot);
  
  valarray<double> sum_L_ij_deriv;
  sum_L_ij_deriv.resize(3*Ntot);
  
  valarray<double> sum_C6ref_Lij_deriv;
  sum_C6ref_Lij_deriv.resize(3*Ntot);
  valarray<double> C6_UCHFGrad;
  C6_UCHFGrad.resize(3*Ntot);
  
  for( int w=0; w<C6_counter-1; w++) {
    double atA, atB, cnA, cnB, c6uchf;
    
    atA  = matrix_Grimme[w][0];
    atB  = matrix_Grimme[w][1];

    cnA  = matrix_Grimme[w][2];
    cnB  = matrix_Grimme[w][3];
    c6uchf  = matrix_Grimme[w][4]; 
    
    if ((atA == atnumA && atB == atnumB) || (atA==atnumB && atB==atnumA)) {
      if (atnumA==atnumB && cnA != cnB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB)))) + exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
	L_ij_deriv = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))))*(-8*(((CN_A-cnA)*vectA)+((CN_B-cnB)*vectB))) + exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))))*(-8*(((CN_A-cnB)*vectA)+((CN_B-cnA)*vectB)));
      } 
      
      else if (atA ==atnumA && atB ==atnumB) {
	L_ij = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))));
	L_ij_deriv = exp(-4*(((CN_A - cnA)*(CN_A - cnA)) + ((CN_B - cnB)*(CN_B - cnB))))*(-8*(((CN_A-cnA)*vectA)+((CN_B-cnB)*vectB)));
      }
      else if (atA==atnumB && atB ==atnumA) {
	L_ij = exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))));
	L_ij_deriv = exp(-4*(((CN_A - cnB)*(CN_A - cnB)) + ((CN_B - cnA)*(CN_B - cnA))))*(-8*(((CN_A-cnB)*vectA)+((CN_B-cnA)*vectB)));
      }
                       
      else {}
      
      sum_L_ij += L_ij;
      sum_C6ref_Lij += L_ij*c6uchf;
      
      sum_L_ij_deriv += L_ij_deriv;
      sum_C6ref_Lij_deriv += L_ij_deriv*c6uchf;
      
      //for (int i=0; i<Ntot; i++ ) {
      //	
      //}

      double numerator = L_ij*c6uchf;
      sum_numerator += numerator;
      
      C6_UCHF = sum_numerator/sum_L_ij;
    }
  }
    
  C6_UCHFGrad = sum_C6ref_Lij_deriv/sum_L_ij - sum_C6ref_Lij*sum_L_ij_deriv/pow(sum_L_ij,2);
  
  return C6_UCHFGrad;
} 

void MP2D::ComputeEnergyandPrintResults() {

  double CutoffRadii[36][36] = {
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	2.1823,	0.0000,	0.0000,	0.0000,	2.5141,	2.4492,	2.3667,	2.1768,	2.0646,	1.9892,	0.0000,	0.0000,	0.0000,	0.0000,	2.8304,	2.6190,	2.4757,	2.3725,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.6026,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	2.5141,	0.0000,	0.0000,	0.0000,	3.2160,	2.9531,	2.7776,	2.6482,	2.6233,	2.4994,	0.0000,	0.0000,	0.0000,	0.0000,	3.1866,	3.0651,	2.9768,	2.9093,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.1029,
    0.0000,	2.4492,	0.0000,	0.0000,	0.0000,	2.9531,	2.9103,	2.7063,	2.5697,	2.4770,	2.4091,	0.0000,	0.0000,	0.0000,	0.0000,	3.1245,	2.9879,	2.8848,	2.8040,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.0108,
    0.0000,	2.3667,	0.0000,	0.0000,	0.0000,	2.7776,	2.7063,	2.6225,	2.4846,	2.3885,	2.3176,	0.0000,	0.0000,	0.0000,	0.0000,	3.0465,	2.9054,	2.7952,	2.7071,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.9227,
    0.0000,	2.1768,	0.0000,	0.0000,	0.0000,	2.6482,	2.5697,	2.4846,	2.4817,	2.3511,	2.2571,	0.0000,	0.0000,	0.0000,	0.0000,	2.8727,	2.8805,	2.7457,	2.6386,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.8694,
    0.0000,	2.0646,	0.0000,	0.0000,	0.0000,	2.6233,	2.4770,	2.3885,	2.3511,	2.2996,	2.1946,	0.0000,	0.0000,	0.0000,	0.0000,	2.7664,	2.7330,	2.6881,	2.5720,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.8109,
    0.0000,	1.9892,	0.0000,	0.0000,	0.0000,	2.4994,	2.4091,	2.3176,	2.2571,	2.1946,	2.1374,	0.0000,	0.0000,	0.0000,	0.0000,	2.6926,	2.6331,	2.5728,	2.5139,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	2.6929,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	2.8304,	0.0000,	0.0000,	0.0000,	3.1866,	3.1245,	3.0465,	2.8727,	2.7664,	2.6926,	0.0000,	0.0000,	0.0000,	0.0000,	3.5017,	3.3180,	3.1916,	3.0982,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.3218,
    0.0000,	2.6190,	0.0000,	0.0000,	0.0000,	3.0651,	2.9879,	2.9054,	2.8805,	2.7330,	2.6331,	0.0000,	0.0000,	0.0000,	0.0000,	3.3180,	3.3107,	3.1523,	3.0352,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.2815,
    0.0000,	2.4757,	0.0000,	0.0000,	0.0000,	2.9768,	2.8848,	2.7952,	2.7457,	2.6881,	2.5728,	0.0000,	0.0000,	0.0000,	0.0000,	3.1916,	3.1523,	3.1046,	2.9730,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.2346,
    0.0000,	2.3725,	0.0000,	0.0000,	0.0000,	2.9093,	2.8040,	2.7071,	2.6386,	2.5720,	2.5139,	0.0000,	0.0000,	0.0000,	0.0000,	3.0982,	3.0352,	2.9730,	2.9148,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.0994,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,
    0.0000,	2.6026,	0.0000,	0.0000,	0.0000,	3.1029,	3.0108,	2.9227,	2.8694,	2.8109,	2.6929,	0.0000,	0.0000,	0.0000,	0.0000,	3.3218,	3.2815,	3.2346,	3.0994,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	0.0000,	3.3662,
  } ;  

  double E6UCHF = 0.0;
  double E6CKS = 0.0;
  double E8UCHF = 0.0;
  double E8CKS = 0.0;
  
  //Define the Gradient
  vector<double> Gradient;
  Gradient.resize(3*Ntot);
  for (int i=0; i<Ntot; i++ ) {
    Gradient[3*i] = 0.0;
    Gradient[3*i +1] = 0.0;
    Gradient[3*i +2] = 0.0;
  }
  
  valarray<double> C6Grad;
  C6Grad.resize(3*Ntot);
  
  for (int i=0; i < NumberOfAtoms; i++) {
    for (int j=i+1; j<Ntot; j++) {
      int atnumA = atomic_number[i];
      int atnumB = atomic_number[j];
      
      double R_AB = ComputeDistance(XYZ[3*i], XYZ[3*j], XYZ[3*i +1], XYZ[3*j +1], XYZ[3*i +2], XYZ[3*j + 2]);
      double ModR = DoubleDamping(CutoffRadii[atnumA][atnumB]*AngToBohr,R_AB, atnumA, atnumB);
      double ModR_D = DoubleDampingDeriv(CutoffRadii[atnumA][atnumB]*AngToBohr,R_AB, atnumA, atnumB);
      
      // C6 contribution
      double damp = TangToennies(6, CutoffRadii[atnumA][atnumB]*AngToBohr, ModR);
      double damp6D = TangToenniesDamp(6, CutoffRadii[atnumA][atnumB]*AngToBohr, ModR);
      
      double C6_UCHF = ComputeC6_UCHF(atnumA, atnumB, Coordination_Number[i], Coordination_Number[j]);
      double C6_CKS = ComputeC6_CKS(atnumA, atnumB, Coordination_Number[i], Coordination_Number[j]);
      
      double E6_UCHF = -1*(C6_UCHF*damp)/(pow(ModR,6));
      double E6_CKS = -1*(C6_CKS*damp)/(pow(ModR,6));
      
      E6UCHF += E6_UCHF;
      E6CKS += E6_CKS;
      
      // C8 contribution
      double damp8D = TangToenniesDamp(8, CutoffRadii[atnumA][atnumB]*AngToBohr, ModR);
      double damp8 = TangToennies(8, CutoffRadii[atnumA][atnumB]*AngToBohr, ModR);
      
      double Q_A = 0.5*sqrt(atomic_number[i])* r2_r4[i];
      double Q_B = 0.5*sqrt(atomic_number[j])* r2_r4[j];
      
      double E8_UCHF = -1*(s_8*3*C6_UCHF*sqrt(Q_A*Q_B)*damp8)/pow(ModR,8);
      double E8_CKS = -1*(s_8*3*C6_CKS*sqrt(Q_A*Q_B)*damp8)/pow(ModR,8);
      
      E8UCHF += E8_UCHF;
      E8CKS += E8_CKS;
      
      // Build the Gradient
      double deriv6_CKS = (E6_CKS/damp)*damp6D*ModR_D/1.8897 + 6*(C6_CKS)/pow(ModR,7)*damp*ModR_D;
      double deriv6_UCHF = (E6_UCHF/damp)*damp6D*ModR_D/1.8897 + 6*(C6_UCHF)/pow(ModR,7)*damp*ModR_D;
      double deriv8_CKS = (E8_CKS/damp8)*damp8D*ModR_D/1.8897 + 8*(s_8*3*C6_CKS*sqrt(Q_A*Q_B))/pow(ModR,9)*damp8*ModR_D;
      double deriv8_UCHF = (E8_UCHF/damp8)*damp8D*ModR_D/1.8897 + 8*(s_8*3*C6_UCHF*sqrt(Q_A*Q_B))/pow(ModR,9)*damp8*ModR_D;
      
      //actual arguments form MP2D gradient
      double argument_x = (deriv6_CKS + deriv8_CKS)-(deriv6_UCHF+deriv8_UCHF);
      double argument_y = (deriv6_CKS + deriv8_CKS)-(deriv6_UCHF+deriv8_UCHF);
      double argument_z = (deriv6_CKS + deriv8_CKS)-(deriv6_UCHF+deriv8_UCHF);  
      
      // C6 Gradients
      valarray<double> C6_UCHF_Deriv = ComputeC6_UCHF_Gradient(atnumA, atnumB, Coordination_Number[i], Coordination_Number[j], CNgrad[i], CNgrad[j]);
      C6_UCHF_Deriv = (-1.0*C6_UCHF_Deriv*damp)/pow(ModR,6);
      
      valarray<double> C8_UCHF_Deriv = ComputeC6_UCHF_Gradient(atnumA, atnumB, Coordination_Number[i], Coordination_Number[j], CNgrad[i], CNgrad[j]);
      C8_UCHF_Deriv = (-1.0*s_8*3.0*C8_UCHF_Deriv*damp8*sqrt(Q_A*Q_B))/pow(ModR,8);
      
      valarray<double> C6_CKS_Deriv = ComputeC6_CKS_Gradient(atnumA, atnumB, Coordination_Number[i], Coordination_Number[j], CNgrad[i], CNgrad[j]);
      C6_CKS_Deriv = (-1.0*C6_CKS_Deriv*damp)/pow(ModR,6);
      
      valarray<double> C8_CKS_Deriv = ComputeC6_CKS_Gradient(atnumA, atnumB, Coordination_Number[i], Coordination_Number[j], CNgrad[i], CNgrad[j]);
      C8_CKS_Deriv = (-1.0*s_8*3.0*C8_CKS_Deriv*damp8*sqrt(Q_A*Q_B))/pow(ModR,8);
      
      C6Grad += ((C6_CKS_Deriv + C8_CKS_Deriv) - (C6_UCHF_Deriv+C8_UCHF_Deriv)); 
      
      Gradient[3*i] += (argument_x/R_AB)*(XYZ[3*i]-XYZ[3*j]);
      Gradient[3*i+1] += argument_y/R_AB*(XYZ[3*i+1]-XYZ[3*j+1]);
      Gradient[3*i+2] += argument_z/R_AB*(XYZ[3*i+2]-XYZ[3*j+2]);
      Gradient[3*j] -= argument_x/R_AB*(XYZ[3*i]-XYZ[3*j]);
      Gradient[3*j+1] -= argument_y/R_AB*(XYZ[3*i+1]-XYZ[3*j+1]);
      Gradient[3*j+2] -= argument_z/R_AB*(XYZ[3*i+2]-XYZ[3*j+2]);
    }
  }
  
  double Total_Energy = ((E6CKS + E8CKS) - (0.5*Cos + 0.5*Css)*(E6UCHF + E8UCHF))*627.509;
  
  cout << "---------------------------------------" << endl;
  cout << "-   MP2D Dispersion correction v1.2   -" << endl;
  cout << "---------------------------------------" << endl;
  cout << "" << endl;
  
  cout << "Chandler Greenwell and Gregory Beran, 2018-2024." << endl;
  cout << "" << endl;
  
  cout << "Please cite:"<< endl;
  cout << "   MP2D: J. Rezac, C. Greenwell, and G. Beran." << endl;
  cout << "   J. Chem. Theory Comput. 14, 4711-4721 (2018).\n   DOI: 10.1021/acs.jctc.8b00548" << endl;
  cout << "" << endl;
  if (SCSMP2D_U) {
    cout << "   SCS-MP2D: C. Greenwell, J. Rezac, and G. Beran." << endl;
    cout << "   Phys. Chem. Chem. Phys. 24, 3695-3712 (2022).\n   DOI: 10.1039/D1CP04922D\n" << endl;
  }
  
  cout << "C6 coefficient files used:" << endl;
  cout << "   UCHF: " << UCHFPATH << endl;
  cout << "   CKS:  " << CKSPATH << endl;
  cout << "" << endl;
  if (SCSMP2D_U) {cout << "SCS-MP2D empirical parameters:" << endl;}
  else
    {cout << "MP2D empirical parameters:" << endl;}
  cout << "   s8    = " << s_8 << endl;
  cout << "   TT_a1 = " << a_one << endl;
  cout << "   TT_a2 = " <<  a_two << endl;
  cout << "   r_cut = " << rcut << endl;
  cout << "   w     = " <<  width << endl;
  if (SCSMP2D_U==true) {
    cout << "   Cos   = " << Cos << endl;
    cout << "   Css   = " <<  Css << endl;
  }
  
  // Print XYZ coordinates
  cout << "\nAtomic Coordinates (Angstroms): " << NumberOfAtoms << " atoms" << endl;
  for (int i =0; i<Ntot; i++) {
    cout << symbols[i];
    printf ("%*f %*f %*f \n", 15, XYZ[3*i]/1.8897, 15, XYZ[3*i+1]/1.8897, 15, XYZ[3*i+2]/1.8897);
  }
  cout << "\n";
  
  //cout << "MP2D Dispersion Energies (kcal/mol)" << endl;
  //printf("   UCHF Contribution:  %.4f \n", (0.5*Cos + 0.5*Css)*(E6UCHF + E8UCHF)*627.5096);
  //printf("   CKS Contribution:  %.4f \n", (E6CKS + E8CKS)*627.5095);
  //printf("   MP2D dispersion correction:  %.4f \n\n",  Total_Energy);
  
  if (SCSMP2D_U)  {
    cout << "SCS-MP2D Dispersion Energies (hartrees)" << endl;
    printf("   HF energy                           = %14.8f\n",HFT);
    printf("   MP2 same-spin correlation           = %14.8f\n",MP2_ss);
    printf("   MP2 opposite-spin correlation       = %14.8f\n",MP2_os);
    double scaled_mp2 = HFT + Cos*MP2_os + Css*MP2_ss;
    printf("   Scaled MP2 energy                   = %14.8f \n", scaled_mp2);
    printf("   Scaled UCHF dispersion contribution = %14.8f \n", (0.5*Cos + 0.5*Css)*(E6UCHF + E8UCHF));
    printf("   CKS dispersion contribution         = %14.8f \n",  (E6CKS + E8CKS));
    printf("   Scaled MP2D dispersion correction   = %14.8f \n", (E6CKS + E8CKS) - (0.5*Cos + 0.5*Css)*(E6UCHF + E8UCHF) );
    double SCSMP2D_Total_Energy = Total_Energy/627.5095 + HFT + Cos*MP2_os + Css*MP2_ss;
    printf("\n ! Total SCS-MP2D energy (hartrees)    = %14.8f \n", SCSMP2D_Total_Energy);
    
  }
  else {
    cout << "MP2D Dispersion Energies (hartrees)" << endl;
    printf("   UCHF dispersion contribution        = %14.8f \n", (0.5*Cos + 0.5*Css)*(E6UCHF + E8UCHF));
    printf("   CKS dispersion contribution         = %14.8f \n",  (E6CKS + E8CKS));
    printf(" ! MP2D dispersion correction          = %14.8f \n", (E6CKS + E8CKS) - (0.5*Cos + 0.5*Css)*(E6UCHF + E8UCHF) );
  }
  
  
  if (Gradient_U == true ) {
    cout << "MP2D Gradient (hartree/bohr)" << endl;
    for (int i=0; i<Ntot; i++ ) {
      printf ("%*f %*f %*f \n", 15, Gradient[3*i]+C6Grad[3*i], 15, Gradient[3*i+1]+ C6Grad[3*i+1], 15, Gradient[3*i+2]+ C6Grad[3*i+2]);
    }
  }
  
}











