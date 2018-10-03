function A = fillA(A,EVS)
% This file is part of GPCCA.
%
% Copyright (c) 2018, 2017 Bernhard Reuter, Susanna R?blitz 
% and Marcus Weber
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
% make A feasible
    %
    % A = fillA( A, EVS )
    %
    % Input:
    %   A       (infeasible) (k,k)-matrix
    %   EVS     (N,k)-matrix with eigen-/Schurvectors of A (for eigenvalues
    %   close to one / for Schurvectors close to the unit matrix)
    %
    % Output:
    %   A       feasible (k,k)-matrix A
    %
    % Written by Susanna Roeblitz and Marcus Weber, Zuse Institute Berlin, 
    % 2012
    % Modified by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017

    [N,k]=size(EVS);
    
    assert( size(A,1)==size(A,2), 'fillA:MatrixShapeError1', ['A matrix '...
        'isnt quadratic!'] )
    assert( size(A,1)==k, 'fillA:MatchError', ['EVS second dimension '...
        'doesnt match A!'] )
    assert( size(A,1)>=2, 'fillA:MatrixShapeError2', ['A matrix must be '...
        'at least 2x2!'] )
    assert(N>=k,'fillA:MatrixShapeError3',['The Eigen- or Schurvector ', ...
        'matrix has more columns than rows. You cant cluster N data- ', ...
        'points into k>N clusters.'])
    assert(N~=k,'fillA:MatrixShapeError4',['The Eigen- or Schurvector ', ...
        'matrix has equal number of columns and rows. There is no ', ...
        'point in clustering N datapoints into k=N clusters.'])
    dummy = ( abs(EVS(:,1) - 1) < ( 100 * eps ) ) ;
    assert(all(dummy(:)), 'fillA:FirstColumnError', 'EVS(:,1) isnt equal 1!')

%       compute 1st column of A by row sum condition
    A(2:k,1) = -sum(A(2:k,2:k),2) ;

%       compute 1st row of A by maximum condition
    for j = 1:k
        A(1,j) = - EVS(1,2:k)*A(2:k,j) ;
        for l = 2:N
            dummy = - EVS(l,2:k)*A(2:k,j) ;
            if (dummy > A(1,j))
%                A(1,j) = max(dummy,0) ;
                A(1,j) = dummy ; % original definition
            end
        end
    end


%       reskale A to be in the feasible set
    A = A/sum(A(1,:)) ;

%       make sure, that there are no zero or negative elements in the first
%       row of A.
    if (any(A(1,:) == 0) == 1)
        A(1,:) % (B.R. 26.07.17)
        error('fillA:ZeroError','first row of A has elements =0!')
    end

    if (min(A(1,:)) < 0)
        A(1,:)
        error('fillA:NegativeError','first row of A has elements <0!')
    end
end
