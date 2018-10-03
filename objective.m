function optval = objective(alpha, EVS, optfile)
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
% compute objective function value
    %
    % optval = objective( alpha, EVS, optfile )
    %
    % Input:
    %   alpha               (k-1)^2-vector with current iterate
    %                       (it columnwise contains A(2:k,2:k))
    %   EVS                 (N,k)-matrix with eigen- or Schur-vectors 
    %                       column by column
    %   optfile             FileID of the output file
    %
    % Output:
    %    optval             current value of the objective function 
    %                       k-trace(S)
    %
    % written by Susanna Roeblitz and Marcus, Zuse Institute Berlin, 
    % Takustrasse 7, 14195 Berlin
    % Modified by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017

    [N,k] = size(EVS) ;
    A = zeros(k,k) ;

%       test for input errors
%   -----------------------------------------------------------------------
    assert(isa(EVS,'double'),'objective:DataTypeError1', ...
        'Variable EVS is type %s not %s.',class(EVS),'double')
    assert(isa(alpha,'double'),'objective:DataTypeError2', ...
        'Variable alpha is type %s not %s.',class(alpha),'double')
    assert(N>=k,'objective:MatrixShapeError1',['The Eigen- or ', ...
        'Schurvector matrix has more columns than rows. You cant ', ...
        'cluster N data-points into k>N clusters.'])
    assert(N~=k,'objective:MatrixShapeError2',['The Eigen- or ', ...
        'Schurvector matrix has equal number of columns and rows. ', ...
        'There is no point in clustering N data-points into k=N clusters.'])
    dummy = ( abs(EVS(:,1) - 1) < ( 100 * eps ) ) ;
    assert(all(dummy(:)), 'objective:FirstColumnError', ...
        'EVS(:,1) isnt equal 1!')
    assert((size(alpha,1) == 1 && size(alpha,2) == (k-1)^2), ...
        'objective:MatrixShapeError3',['alpha is not a (1 x (%d-1)^2) ' ...
        'matrix but (%d x %d).'], k, size(alpha,1), size(alpha,2))
%   -----------------------------------------------------------------------

%       rebuild transformation matrix A
    for i=1:k-1
        for j=1:k-1
            A(i+1,j+1)=alpha(j + (i-1)*(k-1)) ;
        end
    end


%       make A feasible
    A = fillA(A, EVS) ;


%       compute value of objective function
    optval = k - trace(diag(1 ./ A(1,:)) * (A' * A)) ;

    % Note: other choices are possible here, e.g.:
    % optval=-trace(log((diag(1./A(1,:))*(A'*A))));  
    %[White/Shalloway, 2009]

%       save actual function value
    fprintf(optfile, '%.16e\n', optval) ;

end
