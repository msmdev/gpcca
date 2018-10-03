function [ iopt, fileid, final_opt ] = initialize_optional( iopt, id )
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
% initialize the workspace parameters for the optional final optimization 
% procedure.
    %
    % [ iopt, fileid, final_opt ] = initialize_optional( iopt, id )
    %
    % Input:
    %   iopt        structure of optional optimization parameters
    %   id          input ID-string to construct the ID-string for naming
    %               of the output files
    %
    % Output:
    %   iopt        (initialized) structure
    %   fileid      ID-string for naming of the output files
    %   final_opt   flag to indicate, if the optional final optimization
    %               shall be performed (1=yes, 0=no)
    % 
    % Written by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017
    
    if isfield(iopt,'solver')
    
        solver = iopt.solver ;
        final_opt = true ;
        [init,iopt] = getopt(iopt,'init',1) ;
        [parallel,iopt] = getopt(iopt,'parallel',0) ;
        
        if (strcmpi(solver,'gauss-newton'))
        
            [maxiter,iopt] = getopt(iopt,'maxiter',100) ;
%               no Display currentlly implemented
            [display,iopt] = getopt(iopt,'display',0) ; 
%               relative (!) tolerance in x (real tolerance related to
%               XSCALE*XTOL)
            [xtol,iopt] = getopt(iopt,'xtol',1e-4) ;
            if xtol > 1e-2
                xtol = 1.e-2 ;
                iopt.xtol = xtol ;
            end
%               factor for the scaling vector 
%               xscale = xscale * ones(size(xguess))
            [xscale,iopt] = getopt(iopt,'xscale',1e-6) ;
            if xscale > xtol
                xscale = xtol ;
                iopt.xscale = xscale ;
            end
            fileid=strcat(id,'-',solver,'-','maxiter','_', ...
                num2str(maxiter),'-','xscale','_',num2str(xscale), ...
                '-','xtol','_',num2str(xtol),'-',numeric_t) ;
        
        elseif (strcmpi(solver,'nelder-mead'))
        
%               maxiter=-1 means off (internal maximum number of iterations
%               of the algorithm is used);
%               typically maxiter equals 2000-100000 for nelder-mead
            [maxiter,iopt] = getopt(iopt,'maxiter',2000) ;
%               if display=1, iterative output (and plot - in case of 
%               Nelder-Mead) is shown;
%               if display=0, then not
            [display,iopt] = getopt(iopt,'display',0) ;
%               typically tolfun=1e-4 for nelder-mead
            [tolfun,iopt] = getopt(iopt,'tolfun',1e-8) ;
%               typically tolx=1e-4 for nelder-mead
            [tolx,iopt] = getopt(iopt,'tolx',1e-8) ;
            fileid=strcat(id,'-',solver,'-','maxiter','_', ...
                num2str(maxiter),'-','tolfun','_',num2str(tolfun), ...
                '-','tolx','_',num2str(tolx),'-',numeric_t) ;
        
        elseif (strcmpi(solver,'levenberg-marquardt'))
        
%               maxiter=-1 means off (internal maximum number of iterations
%               of the algorithm is used);
%               typically maxiter ~500 for levenberg-marquardt
            [maxiter,iopt] = getopt(iopt,'maxiter',500) ;
%               if display=1, iterative output (and plot - in case of 
%               Nelder-Mead) is shown;
%               if display=0, then not
            [display,iopt] = getopt(iopt,'display',0) ;
%               typically tolfun=1e-8 for levenberg-marquardt
            [tolfun,iopt] = getopt(iopt,'tolfun',1e-8) ;
%               typically tolx=1e-10 for levenberg-marquardt
            [tolx,iopt] = getopt(iopt,'tolx',1e-10) ;
            fileid=strcat(id,'-',solver,'-','maxiter','_', ...
                num2str(maxiter),'-','tolfun','_',num2str(tolfun), ...
                '-','tolx','_',num2str(tolx),'-',numeric_t) ;
        
        else
        
            final_opt = 0 ;
            error('initialize_optional:SolverError', ...
                'Unknown solver for optimization problem.')
        
        end
        
    else
        final_opt = false ;
        fileid = '' ;
    end

end
