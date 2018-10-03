function minChi_vec = use_minChi( EVS, kmin, kmax, fileid )
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
% calculate minChi for cluster numbers between kmin and kmax 
    %
    % use_minChi( EVS, kmin, kmax, fileid )
    %
    % decide for a good interval of clusternumbers to use in the next 
    % search step: the optimal value for minChi can be at most zero
    % therefore you should pick a interval were minChi is closest to zero
    %
    % Input:
    %   EVS         eigen-/Schur-vectors
    %   kmin        minimal number of clusters
    %   kmax        maximum number of clusters
    %   fileid      ID-string for naming of the output files
    %
    % Output:
    %   minChi_vec  vector of minChi values for numbers of clusters from
    %               kmin till kmax
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
%       make sure kmin and kmax are defined
    assert(~isempty(kmin),'use_minChi:EmptyInput1', ...
        'Variable kmin is empty') ;
    assert(~isempty(kmax),'use_minChi:EmptyInput2', ...
        'Variable kmax is empty') ;
%       make sure kmin, kmax are integer values and kmin<=kmax
    assert(mod(kmin,1) == 0, 'use_minChi:InputError1', ...
        'kmin is not an integer value') ;
    assert(mod(kmax,1) == 0, 'use_minChi:InputError2', ...
        'kmax is not an integer value') ;
    assert(kmax >= kmin, 'use_minChi:InputError3', 'kmax !>= kmin') ;
%   -----------------------------------------------------------------------
    disp (' ')
    disp ('Calculate minChi vs. number of clusters')
    disp ('=============================================')
    disp (' ')

    minChi_vec=zeros(kmax-kmin+1,2) ;

    for kt=kmin:kmax
        [ chi, minChi ] = cluster_by_isa( EVS(:,1:kt) ) ;
        name=strcat(fileid,'-','n=',num2str(kt),'-','chi-by-ISA','.txt') ;
        save_t(name,chi,'-ascii')
        minChi_vec(kt-kmin+1,1)=kt ;
        minChi_vec(kt-kmin+1,2)=minChi ;
        fig1 = figure ;
        plot(chi)
        title(['n=' int2str(kt) '; minChi=' num2str(minChi)])
        xlabel('number of clusters')
        ylabel('Chi')
        name=strcat(fileid,'-','n=',num2str(kt),'-', ...
            'chi-by-ISA-figure') ;
        savefig(fig1,strcat(name,'.fig'),'compact')
        print(fig1,strcat(name,'.pdf'),'-dpdf')
    end

    name=strcat(fileid,'-','minChi_vec','-','n=',num2str(kmin), ...
        '-',num2str(kmax),'.txt') ;
    save_t(name,minChi_vec,'-ascii')

    fig2 = figure ; 
    plot(kmin:kmax,minChi_vec(:,2),'-x')
    xlabel('number of clusters')
    ylabel('Chi')
    name=strcat(fileid,'-','minChi-n=',num2str(kmin),'-', ...
        num2str(kmax)) ;
    h = gcf ; % Current figure handle
    set(h,'PaperPositionMode','manual') ;
    set(h,'PaperUnits','inches') ;
    set(h,'PaperPosition',[0 0 3.33 2.63]) ;
    set(h,'PaperUnits','inches') ;
    set(h,'PaperSize',[3.33 2.63]) ;
    savefig(fig2,strcat(name,'.fig'),'compact')
    print(fig2,strcat(name,'.pdf'),'-dpdf')
    print(fig2,strcat(name,'.png'),'-dpng','-r600')
    
    %set(0,'DefaultFigureVisible','on');  % all subsequent figures visible

    disp (' ')
    disp(['Choose a interval were the value in Figure 2 is as ' ...
        'close to zero as possible!'])
    disp (' ')
        
end
