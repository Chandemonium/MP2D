Introduction
============

MP2D is an open-source code (MIT license, see below) for calculating 
van der Waals dispersion for second-order Moller-Plesset pertubation 
theory. The calculation of the MP2D dispersion correction is analogous 
to MP2C:

	E_MP2D = E_MP2 - E_UCHF + E_CKS 

The UCHF and CKS dispersion energies are calculated according to Grimme's D3 
method. The CKS coefficients identical to those used in D3. The UCHF 
coefficients were calculated specifically for this work. Dispersion 
coefficients are currently available for the following elements:

	H, B, C, N, O, F, Ne, P, S, Cl, Ar, Br

NOTE: This code does not compute MP2 itself, only the dispersion correction 
E_Dispersion = -E_UCHF + E_CKS.  You must compute the MP2 energy separately 
yourself with another quantum chemistry software package.

MP2D employs 5 global, empirical parameters:
        s8 = scaling factor for C8 dispersion contribution
        a1 = slope parameter for Tang-Toennies damping function
        a2 = intercept parameter for Tang-Toennies damping function
        rcut = distance below which the secondary short-range damping begins.
        w = width over which short-range doble damping heads to zero.

Default values of these parameters are specified in a global parameter file, 
but the values can be altered at run time as described below.

The MP2D dispersion coefficients were computed in a large augmented 
quadruple-zeta basis (the same as used in Grimme's 2010 D3 paper) and the 
empirical parameters fitted using extrapolated complete basis set MP2 and 
CCSD(T) results. Care should be taken when employing the correction to 
correct MP2 energies computed in smaller basis sets.


Citation
========
Please cite the following paper in any work resulting from use of this
software: 

J. Rezac, C. Greenwell, G. Beran, "Accurate non-covalent interactions via 
dispersion-corrected second-order Moller-Plesset perturbation theory." 
J. Chem. Theory Comput. 14, 4711-4721 (2018). DOI: 10.1021/acs.jctc.8b00548


Installation
============

1) Unpack the .tar.gz file in a directory of your choosing. Edit the 
makefile to set your C++ compiler as appropriate. Type "make" into the 
command line to compile. The executable is named MP2D. Make sure this 
executable is available in your PATH.

The code utilizes several key parmameter files:
   CKS_C6coeffs.dat    - contains the CKS dispersion coefficients
   UCHF_C6coeffs.dat   - contains the UCHF dispersion coefficients
   mp2d_parameters.dat - contains the empirical parameters and the 
			 locations of the CKS and UCH coefficient files

2) Modify the mp2d_parameters.dat file to set the locations where you 
have the dispersion coefficient files:

   CKS_Path: /the/directory/where/MP2D/was/unpacked/CKS_C6coeffs.dat
   UCHF_Path: /the/directory/where/MP2D/was/unpacked/UCHF_C6coeffs.dat

3) Set an environment variable MP2D_PARAM_PATH in your .cshrc or .bashrc
that points to the directory containing mp2d_parameters.dat.

1-3) Alternately, steps 1-3 above can be replaced by `git clone`,
`cmake` configure as outlined in CMakeLists.txt, `make`, and `make
install`. The `mp2d` executable built in this manner knows where its
parameter files are and requires no further editing or envvar setting.

4) Test the code by running the sample jobs in the examples/ folder.  
See the Usage section below for details.

Note) Continuous integration tests and other usual repository
paraphernalia found at https://github.com/MolSSI/QCEngine. A python
interface available in QCEngine v0.7.0. Conda package available
https://anaconda.org/psi4/mp2d.

Usage
======
When the program is executed a dispersion correction energy is 
returned. To compute the dispersion correction for the intermolecular 
interaction energy in a dimer, one would make three inputs; a dimer, 
monomer A (MonA), and monomer B (MonB). A dispersion energy would be 
calculated for each input. E_Interaction = E_Dimer - (E_MonA + E_MonB). 
The input geometry is specified according to standard XYZ Cartestian 
coordinate file format, with coordinates expressed in Angstroms. An 
error will be returned if atoms are present that do not currently have 
dispersion coefficients, or if the input differs from the xyz format.

Input File Format
=================
<number of atoms>
<a line for comments or a blank line>
atom1 x1 y1 z1
atom2 x2 y2 z2
.
.
.
atomn xn yn zn
==================

A number of keywords can be used to control aspects of the MP2D 
dispersion energy calculation. The keywords should come after the input 
file, and be entered as given below:

--PARAM=<path> Specify a path to a customized mp2d_parameters.dat file

--uchf=<path/filename>  Specify a path to a customized UCHF coefficient file
    Use this format: /a/path/to/your/coefficient/file/UCHF_Coefficient_File.dat

--cks=<path/filename>   Specify a path to a customized CKS coefficient file
    Use	this format: /a/path/to/your/coefficient/file/CKS_Coefficient_File.dat

--w=<double>   Specify an alternative value for the width parameter 
		used by the double damping function

--rcut=<double> Specify an alternative value for the cutoff radius parameter 
		 used by the double damping function

--TT_a1=<double>  Specify an alternative value for the a1 parameter used by the
		   Tang-Toennies damping function

--TT_a2=<double>  Specify an alternative value for the a2 parameter use by the
		   Tang-Toennies damping function

--s8=<double>  Specify an alternative value for the s8 parameter used in the 
		energy expression to compensate for the lack of higher oreder terms

--gradient  When specified, compute and print the nuclear gradient of the dispersion energy

Examples: 
	MP2D benzene_dimer.in --gradient
	MP2D benzene_dimer.in --uchf=UCHF_Coefficient_File.dat


License
=======

MIT License

Copyright (c) [2018] [Chandler Scott Greenwell]

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




