function plotmatrix( M, name )
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
%Plot a matrix with logarithmic scaling.
%   The function plots the input matrix M as checkerboard color plot, with
%   logarithmic color scaling. Negative elements are colored red, zero
%   elements are colored white.
%   The plot is outputted as .fig and .pdf file.
%
%   plot_matrix( M, name )
%   
%   Input:
%       M       Matrix to plot
%       name    string for naming of the output files
%
%   Written by Bernhard Reuter, Theoretical Physics II, University of
%   Kassel, 2017
%   -----------------------------------------------------------------------
%       Global variable 'figureVisible' defines, if figures are visible.
%       Even if figures are invisible (figureVisible=false), the figures
%       are still generated and printed to files.
    global figureVisible
    if figureVisible == false
            % all subsequent figures invisible
        set(0,'DefaultFigureVisible','off') ;  
    end
%   -----------------------------------------------------------------------
%       prepare matrix for logarithmic surface plotting
    [m,n] = size(M) ;
    xstep = ceil(n/4) ;
    ystep = ceil(m/4) ;
%       minimum of M
    Mmin = min(M(:)) ;
%       get indices of matrix elements >0
    idx = (M>0) ;
%       get indices of matrix elements =0
    idx1 = (M==0) ;
%       get indices of matrix elements <0
    idx2 = (M<0) ;
%       get minimum of the nonnegative matrix elements of M
    Mminpos = min(min(M(M~=0))) ;
%       initialize ML
    ML = M ;
%       fill ML with log10(m_ij) were m_ij>0
    ML(idx) = log10(M(idx)) ;
%       get minimum of ML
    MLmin = min(min(ML(idx))) ; 
%       get maximum of ML
    MLmax = max(max(ML(idx))) ;
%       if Mminpos<min(ML) set the minimum of the colorbar ticks to
%       floor(Mminpos), else set it to floor(min(ML))
    if Mminpos < MLmin
        ticksmin = floor(Mminpos) ;
    else 
        ticksmin = floor(MLmin) ;
    end
%       set the maximum of the colorbar ticks to ceil(max(ML))
    ticksmax = ceil(MLmax) ;
%       if ticksmin and ticksmax ==0, set ticksmin=-1 and set the lower
%       limit of the colorbar scala to limit=-1.
%       The purpose of this is to account for matrices with only elements
%       between 1 and 0, since log10(1)=0.
%       Else set the lower limit
%       of the colorbar to limit=ticksmin-(ticksmax-ticksmin)/(2^16-1).
%       The purpose of this is to account for m_ij=0, which shall be
%       colored white later.
    if (ticksmax == 0 && ticksmin == 0)
        ticksmin = -1 ;
        limit = -1 ;
    else
        limit = ticksmin - (ticksmax-ticksmin)/(2^16-1) ;
    end
%       set the matrix elements m_ij=0 equal the lower limit of the color
%       scala.
    ML(idx1) = limit ;
%       if there are negative matrix elements, set the lower colorbar limit
%       even deeper
    if any(idx2(:))
        limit = limit - (ticksmax-limit)/(2^16-1) ;
%           set matrix elements m_ij<0 equal the new deeper lower limit of 
%           the color scala
        ML(idx2) = limit ;
    end
%       determine 3 colorbarticks between ticksmin and ticksmax, include
%       ticksmin and ticksmax to get 5 colorbar ticks
    tickstep = (ticksmax-ticksmin)/4 ;
    ticks = ticksmin:tickstep:ticksmax ;
    if ticks(length(ticks)) ~= ticksmax
        ticks = [ ticks, ticksmax ] ;
    end
%       round the tick values to two significant digits
    ticks = round(ticks,2,'significant') ;
%       set the ticklabels to decimal values
    ticklabelsvec = 10.^ticks ;
%       convert the numeric ticklabels to a cell array of strings with
%       scientific format and two significant digits (1 digit behind the
%       decimal point)
    ticklabels = {} ;
    for i = 1:length(ticklabelsvec)
        ticklabels{i} = num2str(ticklabelsvec(i),'%.1e') ;
    end 
%       plot Pc matrix
    f = figure('units','inches') ;
    imagesc(ML) ;
    set(gca,'fontsize',10,'XTick',[1:xstep:n-1 n],'YTick',[1:ystep:m-1 m]) ;
%       define a colormap with 16bit colors scaling
    c = colormap(winter(2^16)) ;
%       set the lowest color scale to white
    if ~any(idx2(:))
        c(1,:) = 1 ;
%       if there are negative matrix elements, set the lowest color scale
%       to red (representing negative matrix elements) and the second 
%       lowest to white (representing zero matrix elements)
    else
        c(1,:) = [ 1, 0, 0 ] ;
        c(2,:) = 1 ;
    end
    colormap(c) ;
    caxis([limit,ticksmax]) ;
%       set the colorbar limits 
    cb = colorbar('LimitsMode','manual','Limits',[ticksmin,ticksmax]) ;
    %cb = colorbar ;
    set(cb,'Ticks',ticks) ;
    set(cb,'LineWidth',1.0) ;
    set(cb,'TickLength',0.02) ;
    set(cb,'Ticklabels',ticklabels) ;
    set(cb,'FontSize',10) ;
%       set the location of the colorbar to the top outside the axes
    set(cb,'Location','northoutside') ;
    %pos = get(gcf,'pos') ;
    h = gcf ; % Current figure handle
    %set(h,'Resize','off');
    set(h,'PaperPositionMode','manual') ;
    set(h,'PaperUnits','inches') ;
    set(h,'PaperPosition',[0 0 3.33 3.8]) ;
    set(h,'PaperUnits','inches') ;
    set(h,'PaperSize',[3.33 3.8]) ;
    title(cb, ['minimum matrix element: ' num2str(min(M(:)),'%.3e')])
    savefig(f,strcat(name,'.fig'),'compact')
    print(f,strcat(name,'.pdf'),'-dpdf')
    
end
