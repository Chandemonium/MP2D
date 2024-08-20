============================================
=                                          =
=    Dispersion-corrected MP2 Library      =
=    Chandler Greenwell & Gregory Beran    =
=                                          =
============================================


================
1. Introduction
================

MP2D is an open-source code (MIT license, see below) for calculating van der Waals dispersion corrections to second-order Moller-Plesset pertubation theory. It can compute the necessary corrections to convert an MP2 energy obtained from any electronic structure package to the MP2D or spin-component-scaled (SCS-MP2D) energies. 


========================
2. Theoretical Overview
========================
The calculation of the MP2D dispersion correction is similar to the MP2C model of Hesselman et al (DOI: 10.1063/1.2905808),

  E(MP2D) = E(MP2) - E(UCHF) + E(CKS)

The key difference is that MP2D employs Grimme's D3 dispersion model (DOI: 10.1063/1.3382344) to calculated the dispersion correction instead of using explicit intermolecular perturbation theory calculations. The CKS coefficients are identical to those used in Grimme's D3. The UCHF coefficients were calculated specifically for the MP2D models. Dispersion coefficients are currently available for the following elements:

  H, B, C, N, O, F, Ne, P, S, Cl, Ar, Br

MP2D employs 5 global, empirical parameters:
  s8 = scaling factor for the C8 dispersion contribution
  a1 = slope parameter for the Tang-Toennies damping function
  a2 = intercept parameter for the Tang-Toennies damping function
  rcut = distance below which the secondary short-range damping begins.
  w = width over which short-range doble damping heads to zero.


SCS-MP2D improves upon MP2 by re-scaling the residual MP2 correlation (i.e. the correlation energy that is not part of the dispersion correction).  The CKS dispersion remains unscaled.  The  SCS-MP2D energy expression is given by:

  E(SCS-MP2D) = E(HF) + Cos*E(MP2_os) + Css*E(MP2_ss) - Cos*E(UCHF_os) - Css*E(UCHF_ss) + E(CKS)

where "os" and "ss" refer to the opposite-spin and same-spin correlation contributions, respectively.


SCS-MP2D employs 7 global empirical parameters. These include the 5 parameters from MP2D plus 2 additional spin-component-scaling parameters:
  Cos = Opposite-spin energy scaling coefficient
  Css = Same-spin energy scaling coefficient



================
3. Installation
================

Installation has been tested on Linux systems.

1) Unpack the .tar.gz file in a directory of your choosing. Edit the makefile to set your C++ compiler as appropriate. Type "make" into the command line to compile. The executable is named MP2D. Make sure this executable is available in your PATH.

The code utilizes several key parmameter files:
   CKS_C6coeffs.dat    - contains the CKS dispersion coefficients
   UCHF_C6coeffs.dat   - contains the UCHF dispersion coefficients
   mp2d_parameters.dat - contains the empirical parameters and the 
			 locations of the CKS and UCHF coefficient files

2) Modify the mp2d_parameters.dat file to reflect the paths of the the dispersion coefficient files on your file system:

   CKS_Path: /the/directory/where/MP2D/was/unpacked/CKS_C6coeffs.dat
   UCHF_Path: /the/directory/where/MP2D/was/unpacked/UCHF_C6coeffs.dat

3) Set an environment variable MP2D_PARAM_PATH in your .cshrc or .bashrc that points to the directory containing mp2d_parameters.dat.

4) Test the code by running the sample jobs in the "examples/" folder.  See the Program Usage section below for details.


==================
4.1 Program Usage
==================

NOTE: This code only computes MP2 dispersion corrections, not the MP2 energy itself. MP2 energies must be computed separately with another quantum chemistry software package. 

For MP2D, the program requires only a Cartesian-coordinate XYZ file as input (using units of Angstroms), and it returns the dispersion correction energy which can be added to the MP2 energy obtained from another software package. To compute the dispersion correction for the intermolecular interaction between two molecules, the user would prepare three XYZ files: the dimer, monomer A (MonA), and monomer B (MonB).  Then run MP2D on each, combine them with their individual MP2 energies, and compute the interaction energy according to:
	E_Interaction = E_Dimer - E_MonA - E_MonB.
An error will be returned if atoms are present that do not currently have dispersion coefficients, or if the input file does not match standard XYZ format.

For SCS-MP2D, the user needs to input the XYZ file together with the HF, same-spin MP2 correlation, and opposite-spin MP2 correlation energies. The three energies are provided via command-line keywords as described in Section 4.3.

Default values of the MP2D and SCS-MP2D empirical parameters are stored in the text files "mp2d_parameters.dat" and "scsmp2d_parameters.dat", while the UCHF and CKS dispersion coefficients are stored in "UCHF_C6coeffs.dat" and "CKS_C6coeffs.dat".  The file paths to those coefficient files are specified in the parameter files.  Most users will never need to alter these parameters.  For expert users who wish to modify the model, it is possible to pass custom parameter sets or coefficient files to the program.

Basis set recommendations: The MP2D and SCS-MP2D dispersion coefficients were computed in a large augmented quadruple-zeta basis (the same basis sets used in Grimme's 2010 D3 paper) and the empirical MP2D/SCS-MP2D parameters were fitted using MP2 and CCSD(T) energies that had been extrapolated to the complete basis set limit (usually via aug-cc-pVTZ and aug-cc-pVQZ two-point extrapolation). For best results, the MP2D and SCS-MP2D dispersion corrections should be applied to large basis MP2 results, preferably those that have been extrapolated to the complete basis set limit.


==========================
4.2 Input XYZ File Format 
==========================
The geometry should be specified in a standard XYZ format, with atomic positions express in Cartesian coordinates in units of Angstroms.

---------------------------------------
<number of atoms>
<one line for comments or one blank line>
atom1 x1 y1 z1
atom2 x2 y2 z2
.
.
.
atomn xn yn zn
---------------------------------------


Example: Water dimer
---------------------------------------
6
Water dimer
O  -1.551007  -0.114520   0.000000
H  -1.934259   0.762503   0.000000
H  -0.599677   0.040712   0.000000
O   1.350625   0.111469   0.000000
H   1.680398  -0.373741  -0.758561
H   1.680398  -0.373741   0.758561
---------------------------------------

=============
4.3 Keywords
=============

Various keywords can be used to control aspects of the MP2D dispersion energy calculation. The keywords should be added after the input file, and be entered as given below.  If no keywords are specified, the code performs a standard MP2D calculation. To run SCS-MP2D, one must used the keywords described below.

--gradient      Compute and print the nuclear gradient of the dispersion energy (MP2D only).



--scsmp2d     Perform SCS-Mp2D instead. The Hartree-Fock, opposite-spin MP2 correlation energy, and
              same-spin MP2 correlation energies must also be specified with the --hf, --mp2os, and --mp2ss keywords.

--hf=<double>  Specify the Hartree-Fock energy (for SCS-MP2D only).

--mp2os=<double>  Specify the opposite-spin MP2 correlation energy (for SCS-MP2D only).

--mp2ss=<double>  Specify the same-spin MP2 correlation energy (for SCS-MP2D only).


Keywords for customizing the parameters:

--param=<path> Specify a path to a customized mp2d_parameters.dat and/or scsmp2d_parameters.dat file
               Format: /a/path/to/your/parameter_files/

--uchf=<path/filename>  Specify a path to a customized UCHF coefficient file
                        Format: /a/path/to/your/coefficient/file/UCHF_Coefficient_File.dat

--cks=<path/filename>   Specify a path to a customized CKS coefficient file
                        Format: /a/path/to/your/coefficient/file/CKS_Coefficient_File.dat

--TT_a1=<double> Specify the a1 parameter used by the Tang-Toennies damping function

--TT_a2=<double> Specify the a2 parameter use by the Tang-Toennies damping function

--w=<double>     Specify the width parameter used by the double damping function

--rcut=<double>  Specify the cutoff radius parameter used by the double damping function

--s8=<double>    Specify the s8 parameter used in the D3 energy expression.

--cos=<double>   Specify MP2 opposite-spin correlation scaling coeffient (SCS-MP2D only).

--css=<double>   Specify MP2 same-spin correlation scaling coeffient (SCS-MP2D only).



		
==========
Examples: 
==========
 1) Evaluate a basic MP2D energy correction:
	MP2D benzene_dimer.xyz
 2) Evaluate the MP2D energy correction and nuclear gradient
	MP2D benzene_dimer.xyz --gradient
 3) Evaluate the MP2D energy correction with a custom set of UCHF coefficients
	MP2D benzene_dimer.xyz --uchf=UCHF_Coefficient_File.dat
 4) Evaluate the SCS-MP2D calculation (passing in the HF and MP2 ss/os correlation energies in hartrees)
        MP2D benzene.xyz --SCSMP2D --hf=-230.7816242 --mp2ss=-0.2399095 --mp2os=-0.7228155
	(The HF, MP2ss, and MP2os energies should be set to the appropriate values for your system).



==========
Citations
==========
Please cite the following MP2D or SCS-MP2D papers in work resulting from use of this software: 

J. Rezac, C. Greenwell, G. Beran, "Accurate non-covalent interactions via dispersion-corrected second-order Moller-Plesset perturbation theory."  J. Chem. Theory Comput. 14, 4711-4721 (2018). DOI: 10.1021/acs.jctc.8b00548

C.Greenwell, J. Rezac, G. Beran, "Spin-component-scaled and dispersion-corrected second-order Moller-Plesset perturbation theory: A path toward chemical accuracy." Phys. Chem. Chem. Phys. 24, 3695-3712 (2022). DOI: 10.1039/D1CP04922D

Separately, a review of the MP2D models can be found in: 

G. Beran, C. Greenwell, C. Cook, and J. Rezac. "Improved Description of Intra- and Intermolecular Interactions Through Dispersion-Corrected Second-Order Moller-Plesset Perturbation Theory." Acc. Chem. Res. 56, 3525-3534 (2023). DOI: 10.1021/acs.accounts.3c00578



========
License
========

MIT License

Copyright (c) [2018-2024] [Chandler Scott Greenwell and Gregory J. O. Beran]

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.



