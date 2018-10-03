function A = initialize_A( evs, init )
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
% initialize the A (n,n)-matrix
    %
    % A = initialize_A( evs, init )
    %
    % Input: 
    %   init    = 0: A will be read from file
    %           = 1: A will be calculated by inverting of the n first 
    %           eigen- or Schurvectors
    %   evs     the n first eigen- or Schurvectors
    %
    % Output: 
    %   A       the A (n,n)-matrix
    %
    % Written by Bernhard Reuter, Theoretical Physics II, 
    % University of Kassel, 2017
%   -----------------------------------------------------------------------
    class_t = class(evs) ;

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
    [N,k]=size(evs) ;
    
    dummy = ( abs(evs(:,1) - 1) < ( num_t('100') * eps(num_t) ) ) ; 
    assert(all(dummy(:)), 'initialize_A:FirstColumnError', ...
        'evs(:,1) isnt equal 1!')
    assert(N>=k,'initialize_A:MatrixShapeError1',['The Eigen- or ', ...
        'Schurvector matrix has more columns than rows. You cant ', ...
        'cluster N datapoints into k>N clusters.'])
    assert(N~=k,'initialize_A:MatrixShapeError2',['The Eigen- or ', ...
        'Schurvector matrix has equal number of columns and rows. ', ...
        'There is no point in clustering N datapoints into k=N clusters.'])

    if init == 1
        
%               search start simplex vertices ('inner simplex algorithm')
        index = indexsearch(evs) ;

%               compute transformation matrix A as initial guess for local
%               optimization (maybe not feasible)
        A = evs(index,:) ;
        dummy = cond(A) ;
        assert((dummy < (num_t('1')/eps(num_t))), ...
            'initialize_A:A_ConditionError', ...
            ['The condition number of the initial guess for A is too ', ...
            'high: cond(A)=%.15e'],double(dummy))
        if (dummy > 1e4)
            disp(['Warning for ' num2str(k) ' clusters: The condition ' ...
                'number of the initial guess for A is > 1e4'])
        end
        disp(' ')
        disp(['condition number of EVS(index,:) for ' num2str(k) ...
            ' clusters: ' num2str(dummy)])
        %disp(num2str(dummy))
        clearvars dummy
        disp(' ')
        A = pinv(A) ;
        A = double(A) ;
            
    elseif init == 0
        
        oldfileid = inputT(['Enter the fileid (IN QUOTES) of the files ' ...
            'containing the initial A to optimize: '],'oldfileid') ;
        name = strcat(oldfileid,'-','n=',int2str(k),'-','A','.mat') ;
        load(name, 'A') ;
            
    else
        
        error('initialize_A:InputError', ...
            'invalid input: init has to be either 0 or 1')
            
    end

end
