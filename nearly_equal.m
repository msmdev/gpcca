function [ a, b, c ] = nearly_equal(a, b)
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
% Checks if two floats a, b are equal with respect to the precision.
%   The precision can bei either mp or double.
%   Based on code from http://www.http://floating-point-gui.de/errors/comparison/
%   Written by Bernhard Reuter, Theoretical Physics II, University of
%   Kassel, 2017
%   -----------------------------------------------------------------------   
%       catch inappropriate use
    assert( (isa(a,'mp') || isa(a,'double')), ...
        'nearly_equal:a_DataTypeError', 'a is neither mp nor double')
    assert( (isa(b,'mp') || isa(b,'double')), ...
        'nearly_equal:b_DataTypeError', 'b is neither mp nor double')    
%   -----------------------------------------------------------------------
    class_t1 = class(a);
    class_t2 = class(b);
    if (strcmp(class_t1,class_t2))
        class_t = class_t1;
    else
        disp(['class(a)=' class_t1 ' is not equal class(b)=' class_t2 '!'])
        disp('will convert to datatype double before comparison ')
        a = double(a) ;
        b = double(b) ;
        class_t = 'double' ;
    end
    
    function r = num_t(expression)
        if (nargin > 0)
            if(strcmp(class_t,'mp')), r = mp(expression);
            else
                if isnumeric(expression)
                    r = expression;
                else
                    r = eval(expression);
                end;
            end;
        else
            r = class_t;
        end;
    end  % num_t
%   -----------------------------------------------------------------------
    
    absA = abs(a) ;
    absB = abs(b) ;
	diff = abs(a - b) ;

    if (a == b) % shortcut, handles infinities
        c = true ;
        return
    elseif (a == num_t('0.0') || b == num_t('0.0') || diff < realmin(num_t))
        % a or b is zero or both are extremely close to it
        % relative error is less meaningful here.
        % In case of standard double floats diff has to be smaller than
        % eps*realmin (as suggested by
        % http://www.http://floating-point-gui.de/errors/comparison/).
        if strcmp(num_t,'double')
            c = ( diff < (eps(num_t) * realmin(num_t)) ) ;
        % In case of multiprecision diff<eps*realmin doesn't work, since
        % realmin('mp') is really the smallest positive floating point
        % number and not the "smallest positive NORMALIZED floating point
        % number": denormal values like 0.00001*realmin!=0 are possible in
        % double floating point precision but not in multiprecision mode -
        % there 0.00001*realmin==0 ALWAYS!
        elseif strcmp(num_t,'mp')
            c = ( diff < realmin(num_t) ) ;
        end
        return
    else % use relative error
        c = (diff / min((absA + absB), realmax(num_t)) < eps(num_t)) ;
        return
    end
end
