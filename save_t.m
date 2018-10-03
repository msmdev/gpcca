function save_t( name, var, fmt )
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
% flexible wrapper function to unify save with different numeric types
    % syntax is basically like in normal matlab save()
    %
    % save_t( name, var, fmt )
    %
    % Input:
    %   name     string containing the filename with ending i.e. 'dummy.txt'
    %   var      variable to be saved i.e. T (NOT in quotes!)
    %   fmt      format specifier i.e. '-ascii'
    %
    % Written by Bernhard Reuter, Theoretical Physics II, 
    % University of Kassel, 2017

    class_t = class(var) ;

    if (strcmpi(class_t,'mp'))
        mp.write(var,name)
    else
        precision = strcat('-',class_t) ;
        save(name,'var',fmt,precision)
    end
end
