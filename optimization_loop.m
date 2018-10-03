function [ Pc, A_cell, chi, val_vec, opt_vec ] = optimization_loop( P, ...
    sd, EVS, A_cell, kmin, kmax, parameters, fileid )
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
%loop over opt_soft() for cluster numbers kt in [kmin,kmax].
%   [ Pc, A_cell, chi, val_vec, opt_vec ] = optimization_loop( P, sd, ...
%       EVS, A_cell, kmin, kmax, parameters, fileid )
%
%   Input:
%       P               (N,N)-matrix to be clustered (row-stochastic)
%       sd              "initial distribution".
%       EVS             (N,k)-matrix with eigen- or Schur-vectors (column 
%                       by column).
%       A_cell          cell vector of initial (k,k)-tranformation 
%                       matrices A.
%       kmin            minimum number of clusters.
%       kmax            maximum number of clusters.
%       parameters      structure with parameters (see gpcca()).
%       fileid          Identification string for naming of the output 
%                       files.
%
%   Output:
%       Pc              Coarse grained stochastic transition matrix as
%                       returned from the last optimization.
%       A_cell          cell vector of (k,k)-tranformation matrices A.
%       chi             (N,k)-matrix with membership values.
%       val_vec         vector of optimality values.
%       opt_vec         vector of crispnesses.
%   Written by Bernhard Reuter, Theoretical Physics II, 
%   University of Kassel, 2017
%   -----------------------------------------------------------------------

    assert(isa(A_cell,'cell'),'optimization_loop:DataTypeError1', ...
        'Variable is type %s not %s', class(A_cell),'cell')
    assert(isa(EVS,'double'),'optimization_loop:DataTypeError2', ...
        'Variable is type %s not %s', class(EVS),'double')
    assert(ischar(fileid),'optimization_loop:DataTypeError3', ...
        'Variable is type %s not %s', class(fileid),'char')
    assert(isstruct(parameters),'optimization_loop:DataTypeError4', ...
        'Variable is type %s not %s', class(parameters),'struct')
    dummy = ( abs(EVS(:,1) - 1) < ( 100 * eps ) ) ;
    assert(all(dummy(:)), 'optimization_loop:FirstColumnError', ...
        'EVS(:,1) isnt equal 1!')
    clearvars dummy

    opt_vec = NaN(kmax-kmin+1,2) ;
    opt_vec(:,1) = kmin:kmax ;
    val_vec = NaN(kmax-kmin+1,2) ;
    val_vec(:,1) = kmin:kmax ;

%       compute optimal solution for every desired number of clusters
     
    if parameters.parallel == 0 %l2
    
        for kt = kmin:kmax %l4
            
            if isempty(A_cell{kt}) == false

                evs = EVS(:,1:kt) ;

%                   pass cluster number to opt_soft for naming the optimization
%                   output plots
                counter = kt ;

%                   call to optimization routine
                [chi, A_cell{kt}, val] = opt_soft(A_cell{kt}, evs, counter, ...
                    fileid, parameters) ;
%               -----------------------------------------------------------
%                   put the actual optimality value 'val' and crispness in two
%                   vectors         
                opt_vec(kt-kmin+1,2) = (kt-val)/kt ;
                val_vec(kt-kmin+1,2) = val ;

%                   postprocess (display, pot, save, calculate) the relevant 
%                   output quantities        
                Pc = postprocess( kt, val, chi, sd, A_cell{kt}, P, fileid ) ;
                
            else
                
                dummy = ['No optimization will be performed for ', ...
                    num2str(kt),' clusters, '] ;
                disp(' ')
                disp(dummy)
                clearvars dummy
                disp('since this would split a 2x2-block of eigenvalues!')
                disp(' ')
                
            end

        end    
        
    elseif parameters.parallel == 1 %l3
        
        disp(' ')
        disp('You are performing the first optimization loop in')
        disp('parallel. Hence the iterations are executed in non-')
        disp('deterministic order. Therefore the output will also')
        disp('occur in non-deterministic order. This is no error')
        disp('and no reason to worry.')
        disp(' ')
%           delete active parallel pool, if there is one
        delete(gcp('nocreate')) ;
%           get the number of physical cores of the machine
        numCores = feature('numcores') ;
%           create a new pool with #workers=#cores
        parpool(numCores) ;
%           create extra variables needed for parfor loop   
        shift_pf = kmin-1 ;
        Pc_cell = cell([1,kmax]) ;
        chi_cell = cell([1,kmax]) ;
        evs_cell = cell([1,kmax]) ;
        for kt = kmin:kmax %l5
            evs_cell{kt} = EVS(:,1:kt) ;
        end
%           perform parallel loop iteration        
        parfor kt = kmin:kmax %l6
            
%               little hack to ensure that, if kmin is an invalid cluster
%               number (leading to splitting of a 2x2-block), the
%               associated opt and val values are NaN. If the following two
%               lines are removed, the opt and val values will be 0 in such
%               a case, which is not intended. This seems to be a MATLAB
%               bug only appearing in parallel by unknown reasons...
            opt_vec(kt-shift_pf,2) ;
            val_vec(kt-shift_pf,2) ;
            
            if isempty(A_cell{kt}) == false
        
                try %l7

                    evs = evs_cell{kt} ;

%                       pass cluster number to opt_soft for naming the 
%                       optimization output plots
                    counter = kt ;

%                       call to optimization routine
                    [chi_cell{kt},A_cell{kt},val]=opt_soft(A_cell{kt}, ...
                        evs, counter, fileid, parameters) ;
%               -----------------------------------------------------------
%                       put the actual optimality value 'val' and crispness in
%                       two vectors         
                    opt_vec(kt-shift_pf,2) = (kt-val)/kt ;
                    val_vec(kt-shift_pf,2) = val ;

%                       postprocess (display, pot, save, calculate) the  
%                       relevant output quantities        
                    Pc_cell{kt} = postprocess( kt, val, chi_cell{kt}, ...
                        sd, A_cell{kt}, P, fileid ) ;
                catch ME %l8
                    disp(' ')
                    disp(['Error in iteration ' num2str(kt) ...
                        ' of the parfors loop in gpcca.m!'])
                    disp('ME.identifier:')
                    disp(ME.identifier)
                    disp('ME.message:')
                    disp(ME.message)
                    rethrow(ME)
                end
            
            else
                
                dummy = ['No optimization will be performed for ', ...
                    num2str(kt),' clusters, '] ;
                disp(' ')
                disp(dummy)
                disp('since this would split a 2x2-block of eigenvalues!')
                disp(' ')
                
            end

        end
        
        Pc = Pc_cell{kmax} ;
        clearvars Pc_cell
        chi = chi_cell{kmax} ;
        clearvars chi_cell
        
    else %l1
        error('optimization_loop:PoolInitializationError',['Couldnt ' ...
            'initialize parallel pool since parameters.parallel was ' ...
            'neither =0 nor =1!'])
    end

end

