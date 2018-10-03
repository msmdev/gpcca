function [ badblocks, badblockrows ] = find_twoblocks( X, RR, fileid )
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
% This function checks the sorted part of the Schurform RR for 2x2-blocks. 
% If a 2x2-block (corresponding to two complex conjugate
% eigenvalues, that MUST NOT be splitted) at positions (rr_i;i, rr_i;i+1,
% rr_i+1;i, rr_i+1;i+1) is found, the element rr_i+1;i is saved and the
% row-index i of the first row of the 2x2-block is saved. This row-index is
% then passed to gpcca.m and is excluded from further analysis.
    % Input:
    %   X           (N,n)-matrix containing the ordered Schur-vectors 
    %               columnwise
    %   RR          (N,N) ordered Schur matrix
    %   fileid      ID-string for the naming of output files
    % Output:
    %   badblocks   (n-1,1)-vector of 0's and 1's indicating cluster
    %               numbers i to exclude (if badblocks(i)=1) from further 
    %               analysis
    %   badblockrows    vector of cluster numbers to exclude
    % Written by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2018

%       check RR for correctness
    assert(size(RR,1)==size(RR,2),'find_twoblocks:MatrixShapeError1', ...
        'RR Matrix is not quadratic but %d x %d',size(RR,1),size(RR,2))
%       check if the number of rows of X matches with the shape of RR
    assert( size(X,1)==size(RR,1), 'find_twoblocks:MatrixShapeError2', ...
        ['The number of rows size(X,1)=%d of X doesnt match with the ', ...
        'shape (%d,%d) of RR!'],...
        size(X,1), size(RR,1), size(RR,2) )

    b = size(X,2) ;
    RR_sorted = RR(1:b,1:b) ;
    subdiagonal = diag(RR_sorted,-1) ;
    name=strcat(fileid,'-subdiagonal-2-',int2str(b),'.txt') ;
    save(name,'subdiagonal','-ascii','-double')
    badblocks = (abs(subdiagonal) > 1E-12) ;
    badblocks = double(badblocks) ;
    name=strcat(fileid,'-badblocks-1-100.txt') ;
    save(name,'badblocks','-ascii','-double')
    badblockrows = find(badblocks) ;
    name=strcat(fileid,'-badblock-row-indices-1-',int2str(b),'.txt') ;
    save(name,'badblockrows','-ascii','-double')

end

