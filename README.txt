This file is part of GPCCA.

Copyright (c) 2020, 2019, 2018, 2017 Bernhard Reuter
with contributions of Marcus Weber and Susanna Röblitz

If you use this code or parts of it, cite the following reference:

Reuter, B., Weber, M., Fackeldey, K., Röblitz, S., & Garcia, M. E. (2018). Generalized
Markov State Modeling Method for Nonequilibrium Biomolecular Dynamics: Exemplified on
Amyloid β Conformational Dynamics Driven by an Oscillating Electric Field. Journal of
Chemical Theory and Computation, 14(7), 3579–3594. https://doi.org/10.1021/acs.jctc.8b00079

GPCCA is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published
by the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU Lesser General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
------------------------------------------------------------------------------------------
For questions or support contact

Bernhard Reuter
bernhard-reuter@gmx.de

or fill a pull request.
------------------------------------------------------------------------------------------

Download the code at: https://github.com/msmdev/gpcca

------------------------------------------------------------------------------------------

Please consider using the further developed Python version of GPCCA - *pyGPCCA* - 
especially, if you are working with large transition matrices with more than a few
thousand rows/columns. pyGPCCA is significantly more performant, if one uses the
optional support for sorted partial Schur decompositions utilizing the SPLEPc library.
Further, GPCCA is not designed to work with sparse matrices. If a sparse matrix is passed
to GPCCA, it is simply densified and thus the RAM consumption may increase dramatically.
On the contrary, pyGPCCA was designed to work equally well with dense and sparse matrices
and makes use of sparse data structures internally whenever possible, if supplied with a
sparse transition matrix.

*pyGPCCA* can be found here:
https://github.com/msmdev/pyGPCCA

Documentation is available here:
https://pygpcca.readthedocs.io/en/latest/

------------------------------------------------------------------------------------------
How to use GPCCA:

Firstly, you will need an installation of MATLAB. Since GPCCA was developed and tested 
with MATLAB version 2015b, it is advised to use this version but later versions should 
work as well.
To execute GPCCA, start MATLAB and run the script main.m. This script contains the 
relevant input parameters, which are preset for the example of butane (the count matrix of 
butane is stored in the file count.txt). The script executes the main program gpcca.m, 
which will ask you to supply the name of a input file containing the count matrix 
characterizing your system/process. This count matrix needs to be generated in advance 
from your data (e.g. molecular dynamics time series data) as described in detail in the 
above reference. There you can also find details about using GPCCA on protein 
conformational dynamics data gained from molecular dynamics simulations. You might also 
take a look at the flowchart supplied in the file Algorithm.png.

The program is pleasant to use interactively and automatically stores a whole range of 
interesting quantities and plots. The following procedure is advised:

1. If interactively selected: Use the minChi criterion for an interval of cluster numbers 
(e.g., 2 to 10) specified in the main.m script to obtain a graph of the minChi's.

2. Based on the minChi graph the interval of cluster numbers (e.g., 2 to 5), within which 
optimization is performed for each cluster number using Gauss-Newton, Nelder-Mead, or 
Levenberg-Marquardt, is selected and inputted. Here it is advisable to use Gauss-Newton, 
since even if it does not converge to the minimum it will come very close to it and is 
extremely fast. A graph of the crispness is outputted, which serves to identify the 
optimal cluster number for a final optimization, if one is desired.
The optimization loop can be executed in parallel (via parfor, by setting the option 
wk.parallel = 1) or serial (wk.parallel = 0). This saves hours and days if you have 
multiple cores and MATLABs Parallel Programming Toolbox.

3. If specified in the main.m script: Interactively, an predetermined optimal cluster 
number (or an interval of predetermined cluster numbers) is entered and the final 
optimization takes place - Nelder-Mead with a ( possibly very) large number of iterations 
has proven to be a good choice here.
The optimization loop can be executed in parallel (via parfor, by setting the option 
iopt.parallel = 1) or serial (iopt.parallel = 0).

There are two test suites that are very easy to use:

The unitTests Testsuite checks the functionality of each individual function, the 
regressionTest Testsuite checks the correct functioning of the entire program (i.e. of 
gpcca.m, which calls all other functions) based on known results for the butane count 
matrix. Extensive test results and reports are stored in folders created by the tests - 
these do not have to be checked manually and are only for archiving (can be deleted).

The test suites are called by run(unitTests) and run(regressionTest) respectively. Then 
the desired precision is chosen in which the test is performed: 'mp' (multiprecision) is 
only selectable if you have installed the Multiprecision Toolbox of Advanpix, 'double' is 
the Matlab standard.

Also for the overall program the precision can and must be determined. This happens in the 
main.m script and is by default set to 'double' there.
WARNING: Don't set a precision bellow 'double'!!! NEVER use single precision!!!
