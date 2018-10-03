function r = load_t( name, fmt, prec )
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
% flexible wrapper function to unify load with different numeric types
    % syntax is basically like in normal matlab load()
    % 
    % r = load_t( name, fmt, prec )
    %
    % Input:
    %   name        string containing the filename with ending
    %   fmt         format specifier '-ascii' - CAUTION: NOT FOR .mat
    %               FILES!
    %   prec        precision to use: 'double' or in case of multiprecision
    %               'mp'
    %
    % Output:
    %   r           loaded data
    %
    % Written by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017

    if (strcmpi(prec,'mp'))
        r = mp.read(name) ;
    elseif strcmpi(prec,'double') && strcmp(fmt, '-ascii')
        r = load(name,fmt) ;
    else
        error('load_t:InputError',['Invalid format or precision given ', ...
            'fmt: ',fmt,', prec=',prec])
    end ;
end
