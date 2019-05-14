#include<iostream>
#include<fstream>
#include<string>
#include<sstream>
#include<vector>
#include<stdio.h>
#include<stdlib.h>
#include<cmath>
#include<valarray>

 
using namespace std;

class MP2D {
    private:
        MP2D();
        ~MP2D();

        int Ntot;
        // Pointers for arrays to store info about each monomer. 
        int *Natoms_per_monomer;
        string *AtomicSymbols;
        string *symbols;
        int *AtomicNumbers;
        int *atomic_number;
        double *r2_r4;
        double *xyz;
        vector<double> XYZ;

        double a_one, a_two, rcut, width, s_8;
     //   double CutoffRadii[36][36];
       // double *r0AB;
        vector<vector<double> > r0AB;
        double **RcovAB;
        int C6_counter;
        int UCHF_counter;

        double AngToBohr;
        double HartreeToKcal;

        vector<vector<double> > matrix_Grimme;
        vector<vector<double> > matrix_UCHF;
        vector<double>  Coordination_Number;
        //vector<vector<double> > CNgrad;
        valarray<valarray<double> > CNgrad;
        vector<vector<double> > dC6_CKSGrad;
        vector<vector<double> > dC6_UCHFGrad;

        
 
    public:

        static MP2D& mp2d() {
            static MP2D themp2d;
            return themp2d;
        }

        // All the flags for the user defined inputs
        bool Gradient_U;
        bool Param_U;
        bool UCHF_U;
        bool CKS_U;
        bool w_U;
        bool rcut_U;
        bool TT_a1_U;
        bool TT_a2_U;
        bool s8_U;
        bool CKS_Coeffs;
        bool UCHF_Coeffs;

        double L_ij;

        // Variables for user defined inputs
        string Param;
        string UCHF_Path;
        string CKS_Path;
        string uchf_coeffs;
        string cks_coeffs;
        double Width;
        double Rcut;
        double TT_a1;
        double TT_a2;
        double s8;

        string Path;
        string GrimmePath;
        string UCHFPath;
        string CKSPATH;
        string UCHFPATH;

        // MP2D Functions
        void Initialize(ifstream &infile, int argc, char *argv[]);
        void GetUserParameters(int argc, char* argv[]);
        void GetParameters(ifstream& infile);
        int GetNumberOfAtoms(int iMon) {return Natoms_per_monomer[0];};
        int BeforeGetTotalNumberOfAtoms(ifstream& infile);
     //   int GetTotalNumberOfAtoms();
        int NumberOfAtoms;
        void GetCoordinatesDifferently(ifstream& infile);
        int *GetAtomicNumber();
        double *GetMultipoleExpectationValue();      
       // void GetCutoffRadii();
        void GetCovalentRadii();
        void GetC6Coefficients();
        double ComputeDistance(double a, double b, double c, double d, double e, double f);
        double Factorial(int n);
        double TangToennies(int n, double r0ab, double Rab);
        double TangToenniesDamp(int n, double r0ab, double Rab);
        double DoubleDamping(double R0ab, double R, int atnumA, int atnumB );
        double DoubleDampingDeriv(double R0ab, double R, int atnumA, int atnumB );
        void CoordinationNumber();
        void CoordinationNumberGradient();
        double ComputeC6_UCHF(int atnumA, int atnumB, double CN_A, double CN_B);
        double ComputeC6_CKS(int atnumA, int atnumB, double CN_A, double CN_B);
        valarray<double> ComputeC6_UCHF_Gradient(int atnumA, int atnumB, double CN_A, double CN_B, valarray<double> vectA, valarray<double> vectB);
        valarray<double> ComputeC6_CKS_Gradient(int atnumA, int atnumB, double CN_A, double CN_B, valarray<double> vectA, valarray<double> vectB);
        void Test_Function();

        
        void Rewind(ifstream& infile) {
            infile.clear();
            infile.seekg(0,ios::beg);
        }


};
