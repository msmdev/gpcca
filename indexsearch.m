function index = indexsearch( evs )
% This file is part of GPCCA.
%
% Copyright (c) 2018, 2017 Bernhard Reuter and Marcus Weber
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
% find a simplex structure in the data
    %
    % index = indexsearch( evs )
    %
    % Input:
    %   evs     (N,k)-matrix with eigen or Schur-vectors columnwise
    %
    % Output:
    %   index   k-vector with indices of objects that build the simplex
    %           vertices
    %
    % Written by Marcus Weber, Zuse Institute Berlin, Takustrasse 7, 
    % 14195 Berlin
    % Modified by Bernhard Reuter, Theoretical Physics II, 
    % University of Kassel, 2017
%   -----------------------------------------------------------------------
    class_t = class(evs);

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
    [N,k] = size(evs) ;
    
    assert(N>=k,'indexsearch:MatrixShapeError',['The matrix has more ', ...
        'columns than rows. You cant get a k-dimensional simplex ', ...
        'from N<k datavectors'])
    
    maxdist = num_t('0.0') ;
    OrthoSys = evs ;
    index = num_t(zeros(1,k)) ;

%       first vertex: row with largest norm
    for l = 1:N
        dist = norm(OrthoSys(l,:)) ;
        if (dist > maxdist)
            maxdist = dist ;
            index(1) = l ;
        end
    end


    OrthoSys = OrthoSys - num_t(ones(N,1)) * evs(index(1),:) ;

%       all further vertices as rows with maximum distance to existing 
%       subspace
    for j = 2:k
        maxdist = num_t('0.0') ;
        temp = OrthoSys(index(j-1),:) ;
        for l = 1:N
            sclprod = OrthoSys(l,:) * temp' ;
            OrthoSys(l,:) = OrthoSys(l,:) - sclprod*temp ;
            distt = norm(OrthoSys(l,:)) ;
            if distt > maxdist  % && ~ismember(l,index(1:j-1))
                maxdist = distt ;
                index(j) = l ;
            end 
        end
        OrthoSys = OrthoSys/maxdist ;
    end
end
