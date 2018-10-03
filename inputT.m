function inp = inputT( inputstring, varname )
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
% Get input from file in testmode.
%
%   inp = inputT( inputstring, varname )
%
%   Get input from keyboard as usual, if program isn't used in testmode.
%   In case the program is in testmode, get the input from file.
%
%   Input:
%       inputstring     the string displayed on screen to invite for input
%       varname         the name-string of the variable to input
%
%   Output:
%       inp             variable with the input from keyboard or file
%
%   Written by Bernhard Reuter, Theoretical Physics II,
%   University of Kassel, 2017
    global testmode
    
%       if testmode=='minChiON', get input from file 
%       'testInput_minChiON.mat': minChi will be used    
    if strcmpi(testmode,'minChiON')
        dummy = load('testInput_minChiON',varname) ;
        inp = getfield(dummy,varname) ;
%       if testmode=='minChiOFF', get input from file 
%       'testInput_minChiOFF.mat': minChi won't be used
    elseif strcmpi(testmode,'minChiOFF')
        dummy = load('testInput_minChiOFF',varname) ;
        inp = getfield(dummy,varname) ;
%       in normal mode just read input from keyboard as usual
    else
        inp = input(inputstring) ;
    end

end

