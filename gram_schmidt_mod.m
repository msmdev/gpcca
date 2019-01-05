function Q = gram_schmidt_mod(EVS,sd)
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
% function to orthonormalize vectors - modified numerically stable version
    %
    % Q = gram_schmidt_mod( EVS, sd )
    %
    % Input:
    %   EVS     (N,N)-matrix with eigen-/Schurvectors of 
    %           Pd=diag(sqrt(sd))*P*diag(1.0./sqrt(sd))
    %   sd      "initial distribution"
    %
    % Output:
    %   Q       (N,N)-matrix with orthonormalized eigen-/Schur-
    %           vectors of Pd and sqrt(sd) on the first column
    %
    % Written by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017,2018,2019
%   -----------------------------------------------------------------------
    class_t1 = class(EVS);
    class_t2 = class(sd) ;
    if (strcmpi(class_t1,class_t2))
        class_t = class_t1 ;
    else
        error('gram_schmidt_mod:EVS_sd_DataTypeError', ...
        ['class(EVS) is not equal class(sd)! this will lead to numeric' ... 
        'presicion errors!'])
    end

    function r = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), r = mp(expression) ;
            else
                if isnumeric(expression)
                    r = expression ;
                else
                    r = eval(expression) ;
                end
            end
        else
            r = class_t ;
        end
    end  % num_t
%   -----------------------------------------------------------------------
    [m,n] = size(EVS) ;
    
    assert(m > 1,'gram_schmidt_mod:MatrixShapeError1', ...
        'The given matrix has only 1 row.') ;
    assert(n > 1,'gram_schmidt_mod:MatrixShapeError2', ...
        'The given matrix has only 1 column.') ;

    Q = num_t(zeros(m,n)) ;

    R = num_t(zeros(n,n)) ;
    
%       search for the constant eigen-/Schurvector, if explicitly present
    max_vsum = 0.0 ;
    for i = 1:n
        vsum = sum(EVS(:,i)) ;
        if ( abs(vsum) > max_vsum )
            max_vsum = abs(vsum) ;
            max_i = i ;
        end
    end

%       keep copy of the original eigen-/Schurvectors for later sanity check
    cEVS = EVS

%       shift non-constant first eigen-/Schurvector to the right
    EVS(:,max_i) = EVS(:,1) ;
%       set first eigen-/Schurvector equal sqrt(sd) 
%       (in do_schur.m the Q-matrix, orthogonalized by gram_schmidt_mod.m,
%       will be multiplied with 1.0./sqrt(sd) - so the first eigen-/Schur-
%       vector will become the unit vector 1!
    EVS(:,1) = num_t(sqrt(sd)) ;

%       sd-orthogonalize
    for j=1:n
        v = EVS(:,j) ;
        for i=1:j-1
            R(i,j) = Q(:,i)' * v ;
            v = v - R(i,j) * Q(:,i) ;
        end
        R(j,j) = num_t(norm(v)) ;
        Q(:,j) = v / R(j,j) ;
    end

%       assert that the subspace didn't change!
    assert(subspace(EVS,cEVS) < (num_t('10000')*eps(num_t)),'gram_schmidt_mod:SubspaceError', ...
        'The derived subspace doesnt match the original one!') ;

end
