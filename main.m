% This file is part of GPCCA.
%
% Copyright (c) 2018, 2017 Bernhard Reuter
%
% If you use this code or parts of it, cite the following reference:
%
% Reuter, B., Weber, M., Fackeldey, K., Röblitz, S., & Garcia, M. E. (2018). Generalized
% Markov State Modeling Method for Nonequilibrium Biomolecular Dynamics: Exemplified on
% Amyloid β Conformational Dynamics Driven by an Oscillating Electric Field. Journal of
% Chemical Theory and Computation, 14(7), 3579–3594. https://doi.org/10.1021/acs.jctc.8b00079
%
% GPCCA is free software: you can redistribute it and/or modify
% it under the terms of the GNU Lesser General Public License as published
% by the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU Lesser General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
% -------------------------------------------------------------------------
% main file for clustering of a row-stochastic matrix by GPCCA
% Written by Bernhard Reuter, Theoretical Physics II,
% University of Kassel, 2017
% -------------------------------------------------------------------------

%   set precision with miltiprecision toolbox for certain numerics
%   uncomment if you want to use multiprecision:
%mp.Digits(50) ;

%   uncomment if you want to use multiprecision:
%disp (['number of digits used in multiprecision numerics (mp): ' ...
    %int2str(mp.Digits)])

%   global variable to define precision to use
global class_t ;

%   set precision of variables or expressions wrapped by numeric_t(),
%   i.e. 'double' or 'mp' for multipresicion
class_t = 'double' ;
disp (' ')
disp (['precision to use in sensitive numerics ' ...
    '(i.e. Eigenvalue and Schur decomposition): ' class_t])

% -------------------------------------------------------------------------

%   Parameters for gpcca
kmin = 2                        % minimum number of clusters
kmax = 10                       % maximum number of clusters
wk.id = 'butane' ;
wk.schur = 1 ;                  % calculate Schurvectors (schur=1) 
                                % or use existing from file (schur=0)
wk.b = 0 ;                      % if b < 0 then -b blocks will be sorted,
                                % if b > 0 then  b or b+1 eigenvalues will 
                                % be sorted, depending on the sizes of 
                                % the blocks,
                                % if b = 0 then the whole Schur form will
                                % be sorted.
wk.init = 1 ;                   % if 1 use A=inv(EVS(index,:)) as starting
                                % guess, 
                                % if =0 read A from file.
wk.solver = 'nelder-mead' ;     % solver for unconstrained optimization 
                                % problem, either 'nelder-mead',
                                % 'levenberg-marquardt', 'gauss-newton'
wk.maxiter = -1 ;
wk.parallel = 0 ;
wk.tolx = 1e-8 ;
wk.tolfun = 1e-8 ;
iopt.init = 2 ;                 % If =1 use A=inv(EVS(index,:)) as starting 
                                % guess, 
                                % if =0 read A from file with identifier 
                                % interactively passed from the command 
                                % window,
                                % if =2 use the the optimized A matrices 
                                % from the  first optimization loop as 
                                % input for the final optimization.
iopt.solver = 'gauss-newton' ;   % solver for optional final optimization
iopt.maxiter = 10 ;
iopt.parallel = 0 ;

% -------------------------------------------------------------------------

%   read the count matrix from file, calculate the stochastic matrix,
%   and call gpcca rotine gpcca.m

%   load the count matrix from file
disp (' ')
COUNTMATRIX = input('Enter the name of the count-matrix file (IN QUOTES): ') ;
Tc = load_t(COUNTMATRIX,'-ascii',class_t) ;
dummy = (mod(Tc,1) ~= 0) ;
assert(~any(dummy(:)), ...
    'main:Tc_DataError','Tc doesnt seem to be a count matrix')
clearvars dummy
assert(isa(Tc,numeric_t),'main:Tc_DataTypeError', ...
    'Variable is type %s not %s',class(Tc),numeric_t)
assert(size(Tc,1)==size(Tc,2),'main:Tc_MatrixShapeError', ...
    'Matrix is not quadratic but %d x %d',size(Tc,1),size(Tc,2))
%   make sure there are now rows with zero rowsum in the count matrix
assert(~any(sum(Tc,2) < numeric_t('0.99')),'main:ZeroRowError', ...
    'Matrix has rows with rowsum zero')

%   calculate stochastic matrix P from the count matrix Tc
P = diag(numeric_t('1.0')./sum(Tc,2)) * Tc ;
assert(isa(P,numeric_t),'main:P_DataTypeError', ...
    'Variable is type %s not %s',class(P),numeric_t)
dummy = (sum(P,2) > numeric_t('0.0')) ;
assert(all(dummy(:)),'ZeroError', 'Not all rows of P are >0!')
clearvars dummy

%   calculate "initial distribution"
sd = sum(Tc,2) ; sd=sd/sum(sd) ;
assert(isa(sd,numeric_t),'main:sd_DataTypeError', ...
    'Variable is type %s not %s',class(sd),numeric_t)
assert(all(sd > numeric_t('0.0')),'ZeroError', 'Not all elements of sd are >0!')

%   perform GPCCA
[Pc, chi, A, wk, iopt] = gpcca(P, sd, kmin, kmax, wk, iopt) ;
disp('parameters for the first optimization procedure:')
disp(wk)
disp('parameters for the optional optimization procedure:')
disp(iopt)
