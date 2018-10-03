function Pc = postprocess( k, val, chi, sd, A, P, fileid )
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
% Postprocess (calculate, display, store) the relevant output quantities.
    %
    % Pc = postprocess( k, val, chi, sd, A, P, fileid )
    %
    % The relevant output quantities like membership vectors chi, starting
    % distribution sd, coarse grained matrix Pc etc. are postprocessed here:
    % They are displayed on screen, plotted and saved to files. Some
    % quantities (like Pc) are even calculated here.
    %
    % Input:      
    %   k           number of clusters
    %   val         optimality value
    %   chi         membership vectors (as (N,k)-matrix)
    %   sd          "initial" distribution
    %   A           (k,k)-tranformation matrix A (chi = EVS*A)
    %   P           stochastic (fine) transition matrix
    %   fileid      ID-string for naming of the output files
    %
    % Output:
    %   Pc          Coarse grained stochastic transition matrix
    %
    % Written by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017
%   -----------------------------------------------------------------------
%       Global variable 'figureVisible' defines, if figures are visible.
%       Even if figures are invisible (figureVisible=false), the figures
%       are still generated and printed to files.
    global figureVisible
    if figureVisible == false
        % all subsequent figures invisible
        set(0,'DefaultFigureVisible','off');  
    end
%   ----------------------------------------------------------------------- 
%       assertions
    assert(size(P,1)==size(P,1), 'postprocess:MatrixShapeError1', ...
        'P isnt quadratic') ;
    N = size(P,1) ;
    assert(mod(k,1) == 0, 'postprocess:InputError1', ...
        'k is not an integer value') ;
    assert(isa(val,'double'), 'postprocess:InputError2', ...
        'val is not of class double') ;
    assert(isa(chi,'double'), 'postprocess:InputError3', ...
        'chi is not of class double') ;
    assert(isa(sd,'double'), 'postprocess:InputError5', ...
        'sd is not of class double') ;
    assert(isa(A,'double'), 'postprocess:InputError6', ...
        'A is not of class double') ;
    assert(isa(P,'double'), 'postprocess:InputError7', ...
        'P is not of class double') ;
    assert(isa(fileid,'char'), 'postprocess:InputError8', ...
        'fileid is not of class char') ;
    assert(all(size(chi)==[N,k]), 'postprocess:MatrixShapeError2', ...
        'chi is of shape (%d,%d) but should be of shape (%d,%d)', ...
        size(chi,1), size(chi,2), N, k) ;
    assert(all(size(sd)==[N,1]), 'postprocess:MatrixShapeError4', ...
        'sd is of shape (%d,%d) but should be of shape (%d,1)', ...
        size(sd,1), size(sd,2), N) ;
    assert(all(size(A)==[k,k]), 'postprocess:MatrixShapeError5', ...
        'A is of shape (%d,%d) but should be of shape (%d,%d)', ...
        size(A,1), size(A,2), k, k) ;
%   -----------------------------------------------------------------------
    crispness = (k-val)/k ;
    disp(' ')
    disp(['For ' int2str(k) ' clusters : crispness = ' num2str(crispness)])
    disp(' ')
    name = strcat(fileid,'-','n=',int2str(k),'-','val','.txt') ;
    save_t(name,val,'-ascii')
    name = strcat(fileid,'-','n=',int2str(k),'-','crispness','.txt') ;
    save_t(name,crispness,'-ascii')

    disp(' ')
    disp('Statistical sd-weights of these clusters')
    sd_weights = chi'*sd ;
    for i = 1:length(sd_weights)
        disp(['For ' int2str(i) '-th cluster : Weight = ' ...
            num2str(sd_weights(i))])
    end
    disp(' ')
    name = strcat(fileid,'-','n=',int2str(k),'-','sd_weights','.txt') ;
    save_t(name,sd_weights,'-ascii')
        
    name = strcat(fileid,'-','n=',int2str(k),'-','chi','.txt') ;
    save_t(name,chi,'-ascii')
    name = strcat(fileid,'-','n=',int2str(k),'-','A','.txt') ;
    save_t(name,A,'-ascii')
    name = strcat(fileid,'-','n=',int2str(k),'-','A','.mat') ;
    save(name,'A')
        
    dummy = cond(chi' * diag(sd) * chi) ;
    assert(dummy < ( 1.0 / eps ),'postprocess:Pc_ConditionError', ...
    'The condition number of chi.transpose*diag(sd)*chi is too high: ', ...
    'cond(A)=%.15e',dummy)
    if (dummy > 1e4)
        disp(['Warning: The condition number of ', ...
            'chi.transpose*diag(sd)*chi is > 1e4'])
    end
    disp(' ')
    disp('condition number of chi.transpose*diag(sd)*chi')
    disp(num2str(dummy))
    clearvars dummy
    disp(['If the condition is large better be careful since ', ...
        'chi.transpose*diag(sd)*chi is inverted!!'])
    disp(' ')
    disp('Coarse-grainded sd-normed transition matrix:')
    Pc = pinv(chi'*diag(sd)*chi)*(chi'*diag(sd)*P*chi)
    disp(' ')
    name = strcat(fileid,'-','n=',int2str(k),'-','Pc','.txt') ;
    save_t(name,Pc,'-ascii')
%       plot Pc matrix with logarithmic scaling
    name = strcat(fileid,'-','n=',int2str(k),'-','Pc') ;
    plotmatrix( Pc, name )
%       plot chi matrix with logarithmic scaling
    name = strcat(fileid,'-','n=',int2str(k),'-','chi') ;
    plotmatrix( chi, name )
        
%       plot membership values chi for all objects
    fig3 = figure ;
    title(['n=' int2str(k)])
    plot(chi)
    xlabel('object')
    ylabel('membership')
    name = strcat(fileid,'-','n=',int2str(k),'-','chi-figure') ;
    savefig(fig3,strcat(name,'.fig'),'compact')
    print(fig3,strcat(name,'.pdf'),'-dpdf')
%       figure(3) might look chaotic if consecutively numbered states are
%       assigned to different clusters (e.g. if you load P_pentan);
%       Therefore: sort states according to cluster assignment and plot 
%       again:

%       sort the states according to the cluster assignment
    [~, numclus] = size(chi) ;
    [~,cluster] = max(chi,[],2) ;
    idx_vec = [] ;
    ticks = [] ;
    ticklabels = {} ;
    for k = 1:numclus
        idx = find(cluster == k) ;
        tick1 = length(idx_vec)+1 ;
        idx_vec = [idx_vec;idx] ;
        tick2 = length(idx_vec) ;
        if tick2 > tick1
            ticks = [ticks,tick1,tick2] ;
            ticklabels = [ticklabels,num2str(idx(1)), ...
                num2str(idx(length(idx)))] ;
        elseif ~isempty(idx)
            ticks = [ticks,tick1] ;
            ticklabels = [ticklabels,num2str(idx(1))] ;
        end
    end

    name = strcat(fileid,'-','n=',int2str(k),'-','idx_vec','.txt') ;
    save_t(name,idx_vec,'-ascii')

    fig4 = figure ;
    plot(chi(idx_vec,:))
    ax = gca ;
    ax.XTickMode = 'manual' ;
    ax.XTick = ticks ;
    ax.XTickLabelMode = 'manual' ;
    ax.XTickLabel = ticklabels ;
    xlabel('object')
    ylabel('membership')
    title(['n=' int2str(k) ', states ordered according to cluster assignment:'])
    name = strcat(fileid,'-','n=',int2str(k),'-', ...
        'chi-ordered-figure') ;
    savefig(fig4,strcat(name,'.fig'),'compact')
    print(fig4,strcat(name,'.pdf'),'-dpdf')
    
    %set(0,'DefaultFigureVisible','on');  % all subsequent figures visible

end
