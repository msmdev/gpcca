function r = numeric_t(expression)
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
% flexible wrapper function to unify work with different numeric types
    % 
    % function r = numeric_t(expression)
    %
    % Input:
    %   expression      An expression like 'pi/2' or a float like '2.0'
    %                   (Quotes are important!). '2' and '2.0' are 
    %                   equivalent.
    %                   Also float matrices (user defined or elementary  
    %                   ones like eye(), ones(), zeros() can be passed but 
    %                   WITHOUT QUOTES!
    %
    % Output:
    %   r               Either the result of the input expression or the 
    %                   passed float or Matrix of floats but in any case of 
    %                   datatype class_t (as defined in main.m!).
    %
    % Written by Bernhard Reuter, Theoretical Physics II, 
    % University of Kassel, 2017 -- based on example code by Pavel 
    % Holoborodko
    % from  https://www.advanpix.com/2016/07/21/how-to-write-precision ...
    % -independent-code-in-matlab/

    global class_t;

    if (nargin > 0)
        if(strcmpi(class_t,'mp')), r = mp(expression);
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
end
