function [ Chi, indic ] = cluster_by_isa( Evs )
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
% classification of dynamic data based on NoOfClus orthonormal eigen- or 
% Schur-vectors of the (stochastic) transition matrix.
% Hereby NoOfClus determines the number of clusters to cluster the data
% into. Applied method is the Inner Simplex Algorithm.
% Constraint: Evs matrix needs to contain at least NoOfClus
% eigen-/Schur-vectors.
%
% [ Chi, indic ] = cluster_by_isa ( Evs ) 
%
% Input:
%   Evs         orthonormal (N,k)-matrix with eigen- or Schur-vectors 
%               (column by column), with the first column constantly equal
%               1
%
% Output:
%   Chi         (N,k)-matrix with membership values
%   indic       minChi indicator, see [1]
%
% References:
%   [1] (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
%   Application to Markov State Models and Data Classification.
%   Adv. Data Anal. Classif. 2013, 7 (2), 147?179.
%
% Written by Marcus Weber, Zuse Institute Berlin
% Modified by Bernhard Reuter, Theoretical Physics II,
% University of Kassel, 2017
%   -----------------------------------------------------------------------
    class_t = class(Evs) ;

    function r = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), r = mp(expression) ;
            else
                if isnumeric(expression)
                    r = expression ;
                else
                    r = eval(expression) ;
                end ;
            end ;
        else
            r = class_t ;
        end ;
    end  % num_t
%   -----------------------------------------------------------------------
    [N,NoOfClus] = size(Evs) ;
    
    dummy = ( abs(Evs(:,1) - 1) < ( numeric_t('100') * eps(numeric_t) ) ) ;
    assert(all(dummy(:)), 'cluster_by_isa:FirstColumnError', ...
        'Evs(:,1) isnt equal 1!')

%       special cases
    if NoOfClus < 2 % l2
        Chi = ones(size(Evs,1),1) ;
        indic = 0 ;
        %RotMat = 1/Evs(1,1) ;
        %cF= ones(size(Evs,1),1) ;
        %ind = [] ;
    elseif NoOfClus == N % l3
        Chi = eye(size(Evs,1)) ;
        indic = 0 ;
        %RotMat = pinv(Evs) ;
        %cF = [1:size(Evs,1)] ;
        %ind = [] ;
    elseif NoOfClus > N % l4
        error('cluster_by_isa:MatrixShapeError',['The Eigen- or ', ...
            'Schurvector matrix has more columns than rows.'])
    else % l1
%           ISA-Algorithm
        C = Evs(:, 1:NoOfClus) ;
        OrthoSys = C ;
        maxdist = num_t('0.0') ;
%           first two representatives with maximum distance
        for i = 1:size(Evs,1) % l5        
            if norm(C(i,:)) > maxdist % l6
                maxdist = norm(C(i,:)) ;
                ind(1)=i ;
            end
        end
        for i = 1:size(Evs,1) % l7
            OrthoSys(i,:) = OrthoSys(i,:) - C(ind(1),:) ;
        end ;
    
%           further representatives through Gram-Schmidt orthogonalization
        for k = 2:NoOfClus % l8
            maxdist = num_t('0.0') ;
            temp = OrthoSys(ind(k-1),:) ;
            for i = 1:size(Evs,1) % l9
                OrthoSys(i,:) = OrthoSys(i,:) - ...
                    (temp * transpose(OrthoSys(i,:))) * temp ;
                distt = norm(OrthoSys(i,:)) ;
                if distt > maxdist % l10
                    maxdist = distt ;
                    ind(k) = i ;
                end
            end
            OrthoSys = OrthoSys ./ norm(OrthoSys(ind(k),:)) ;
        end
        
%           linear transformation of the eigen-/Schur-vectors
        size_of_C = size(C(ind,:)) % print 
        condition_of_C = cond(C(ind,:)) % print 
        RotMat = pinv(C(ind,:)) ;
        Chi = C * RotMat ;
        Chi = double(Chi) ;
%           determine indicator
        indic = min(min(Chi)) ;
%           defuzzification of the membership functions
%        [minVal, cF] = max(transpose(Chi)) ;
    end
end
