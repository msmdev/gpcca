function [ chi, A, fopt ] = opt_soft( A, evs, counter, fileid, wk )
% This file is part of GPCCA.
%
% Copyright (c) 2018, 2017 Bernhard Reuter, Susanna R?blitz 
% and Marcus Weber
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
% computation of membership matrix by PCCA+
    %
    % [ chi, A, fopt ] = opt_soft( A, evs, counter, fileid, wk )
    %
    % Input:
    %   A           initial (k,k)-tranformation matrix A
    %   evs         (N,k)-matrix with eigen- or Schur-vectors (column by 
    %               column)
    %   counter     k, actual cluster number
    %   fileid      Identification string for naming of the output files.
    %   wk          structure. See declaration of gpcca().
    %
    % Output:
    %   chi         (N,k)-matrix with membership values
    %   A           (k,k)-tranformation matrix A s.th. chi=evs*A
    %   fopt        optimal function value
    %
    % cite:
    %
    % [1] P. Deuflhard, M. Weber: Robust Perron Cluster Analysis in 
    %     Conformation Dynamics. Lin. Alg. App. 2005, 398c, Special 
    %     issue on matrices and mathematical biology, 161-184.
    % [2] S. Roeblitz and M. Weber: Fuzzy Spectral Clustering by PCCA+.
    %     In: Classification and Clustering: Models, Software and 
    %     Applications, WIAS Report No. 26, Berlin 2009.
    %
    % written by Bernhard Reuter, Theoretical Physics II, 
    % University of Kassel, 2017,
    % based on algorithms and code by Susanna Roeblitz and Marcus Weber, 
    % Zuse Institute Berlin, Takustrasse 7, 14195 Berlin
%   -----------------------------------------------------------------------

    [N,k]=size(evs) ;

    assert(isa(A,'double'),'opt_soft:DataTypeError1', ...
        'Variable is type %s not %s', class(A),'double')
    assert(isa(evs,'double'),'opt_soft:DataTypeError2', ...
        'Variable is type %s not %s', class(evs),'double')
    assert(ischar(fileid),'opt_soft:DataTypeError3', ...
        'Variable is type %s not %s', class(fileid),'char')
    assert(isstruct(wk),'opt_soft:DataTypeError4', ...
        'Variable is type %s not %s', class(wk),'struct')
    assert( size(A,1)==size(A,2), 'opt_soft:MatrixShapeError1', ...
        'A matrix isnt quadratic!')
    assert( size(A,1)==k, 'opt_soft:MatchError', ...
        'EVS second dimension doesnt match A!')
    assert( size(A,1)>=2, 'opt_soft:MatrixShapeError2', ...
        'A matrix must be at least 2x2!')
    assert(N>=k,'opt_soft:MatrixShapeError3',['The Eigen- or ', ...
        'Schurvector matrix has more columns than rows. You cant ', ...
        'cluster N datapoints into k>N clusters.'])
    assert(N~=k,'opt_soft:MatrixShapeError4',['The Eigen- or ', ...
        'Schurvector matrix has equal number of columns and rows. ', ...
        'There is no point in clustering N datapoints into k=N clusters.'])
    dummy = ( abs(evs(:,1) - 1) < ( 100 * eps ) ) ;
    assert(all(dummy(:)), 'opt_soft:FirstColumnError', ...
        'evs(:,1) isnt equal 1!')
        
%       first, overwrite eventually existing file with values of the
%       objective function gained during optimization
    name = strcat(fileid,'-','n=',num2str(counter),'-','objective','.txt') ;
    optfile = fopen(name, 'w') ;

%       reduce optimization problem to size (k-1)^2
    alpha = zeros(1,(k-1)^2) ;
    for i = 1:k-1 % l3
        for j = 1:k-1 % l4
            alpha(j + (i-1)*(k-1)) = A(i+1,j+1) ;
        end
    end


%       perform optimization
    switch lower(wk.solver)

        case 'nelder-mead' % l5
            
            if wk.maxiter == -2
                wk.maxiter = (50 * (k-1)^2) ;
            end

            if ((wk.maxiter > 0) && (wk.display == false)) % l6+l7

                options=optimset('maxfunevals',wk.maxiter,'MaxIter', ...
                    wk.maxiter,'TolFun',wk.tolfun,'TolX',wk.tolx) ;

            elseif ((wk.maxiter <= 0) && (wk.display == true)) % l8+l9

                options = optimset('Display','iter','PlotFcns', ...
                    @optimplotfval,'TolFun', wk.tolfun,'TolX',wk.tolx) ;

            elseif ((wk.maxiter > 0) && (wk.display == true)) % l10+l11

                options = optimset('Display','iter','PlotFcns', ...
                    @optimplotfval, 'maxfunevals',wk.maxiter,'MaxIter', ...
                    wk.maxiter,'TolFun',wk.tolfun,'TolX',wk.tolx) ;

            else % l5

                options = optimset('TolFun',wk.tolfun,'TolX',wk.tolx) ;

            end
%               optimization
            [alpha,fopt,exitflag,output] = fminsearch(@(alpha) ...
                objective(alpha,evs,optfile),alpha,options) ;
%               save optimization plot, if option wk.display==true
            if wk.display == true % l12
                name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                    'optimization','.fig') ;
                savefig(name)
            end
%               save options structure to .mat file 
            name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                'optimization-options','.mat') ;
            save(name,'options')
%               save output, exitflag and certain options to text file
            name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                'optimization-output','.txt') ;
            T = struct2table(output) ;
            writetable(T,name)
            fid = fopen(name,'a') ;
            fprintf(fid,'exitflag: %d\n\n',exitflag) ;
            fprintf(fid,'%s\n','options:') ;
            fprintf(fid,'Display: %s\n',options.Display) ;
            fprintf(fid,'MaxFunEvals: %d\n',options.MaxFunEvals) ;
            fprintf(fid,'MaxIter: %d\n',options.MaxIter) ;
            fprintf(fid,'TolFun: %.4e\n',options.TolFun) ;
            fprintf(fid,'TolX: %.4e\n',options.TolX) ;
            if isa(options.PlotFcns,'function_handle')
                fprintf(fid,'PlotFcns: %s\n',func2str(options.PlotFcns)) ;
            else
                fprintf(fid,'PlotFcns: %s\n',options.PlotFcns) ;
            end
            fclose(fid) ;
            

        case 'levenberg-marquardt' % l13

%            MAXFUNEVALS = 200*2*((k-1)*(k-1)) ;
            if ((wk.maxiter > 0) && (wk.display == false)) % l14+l15

                options=optimoptions('lsqnonlin','Algorithm', ...
                    'levenberg-marquardt','MaxIter',wk.maxiter, ...
                    'TolFun',wk.tolfun,'TolX',wk.tolx) ;
                options.ScaleProblem = 'jacobian' ;

            elseif ((wk.maxiter <= 0) && (wk.display == true)) % l16+l17

                options = optimoptions('lsqnonlin','Algorithm', ...
                    'levenberg-marquardt','Display','iter','TolFun', ...
                    wk.tolfun,'TolX',wk.tolx) ;
                options.ScaleProblem = 'jacobian' ;

            elseif ((wk.maxiter > 0) && (wk.display == true)) % l18+l19

                options = optimoptions('lsqnonlin','Algorithm', ...
                    'levenberg-marquardt','Display','iter','MaxIter', ...
                    wk.maxiter,'TolFun',wk.tolfun,'TolX',wk.tolx) ;
                options.ScaleProblem = 'jacobian' ;
%                options.MaxFunEvals = MAXFUNEVALS ;

            else % l13

                options = optimoptions('lsqnonlin','Algorithm', ...
                    'levenberg-marquardt', 'TolFun',wk.tolfun, ...
                    'TolX',wk.tolx) ;
                options.ScaleProblem = 'jacobian' ;

            end
%               optimization                    
            [alpha,~,fopt,exitflag,output] = lsqnonlin(@(alpha) ...
                objective(alpha,evs,optfile),alpha,[],[],options) ;
%               save options structure to .mat file 
            name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                'optimization-options','.mat') ;
            save(name,'options')
%               save output, exitflag and certain options to text file
            name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                'optimization-output','.txt') ;
            T = struct2table(output,'AsArray',true) ;
            writetable(T,name)
            fid = fopen(name,'a') ;
            fprintf(fid,'\nexitflag: %d\n\n',exitflag) ;
            fprintf(fid,'%s\n','options:') ;
            fprintf(fid,'Algorithm: %s\n',options.Algorithm) ;
            fprintf(fid,'Display: %s\n',options.Display) ;
            %fprintf(fid,'MaxFunEvals: %d\n',options.MaxFunEvals) ;
            fprintf(fid,'MaxIter: %d\n',options.MaxIter) ;
            fprintf(fid,'TolFun: %.4e\n',options.TolFun) ;
            fprintf(fid,'TolX: %.4e\n',options.TolX) ;
            fprintf(fid,'ScaleProblem: %s\n',options.ScaleProblem) ;
            fclose(fid) ;

        case 'gauss-newton' % l20

            if k > 2 % l21
                par.evs = evs ;
                problem = 'problem_pcca_nlscon' ;
                alpha = main_nlscon(alpha, par, problem, wk.maxiter, ...
                    wk.xscale, wk.xtol, fileid, counter) ;
                fopt = feval(problem,alpha, '', par) ;
            else % l20
                options = optimset('MaxIter',2000,'TolFun',1e-8, ...
                    'TolX',1e-8) ;
                [alpha,fopt,exitflag,output] = fminsearch(@(alpha) ...
                    objective(alpha,evs,optfile),alpha,options) ;
%               save options structure to .mat file 
                name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                    'optimization-options','.mat') ;
                save(name,'options')
%                   save output, exitflag and certain options to text file
                name = strcat(fileid,'-','n=',num2str(counter),'-', ...
                    'optimization-output','.txt') ;
                T = struct2table(output) ;
                writetable(T,name)
                fid = fopen(name,'a') ;
                fprintf(fid,'exitflag: %d\n\n',exitflag) ;
                fprintf(fid,'%s\n','options:') ;
                fprintf(fid,'Display: %s\n',options.Display) ;
                fprintf(fid,'MaxFunEvals: %d\n',options.MaxFunEvals) ;
                fprintf(fid,'MaxIter: %d\n',options.MaxIter) ;
                fprintf(fid,'TolFun: %.4e\n',options.TolFun) ;
                fprintf(fid,'TolX: %.4e\n',options.TolX) ;
                if isa(options.PlotFcns,'function_handle')
                    fprintf(fid,'PlotFcns: %s\n',func2str(options.PlotFcns)) ;
                else
                    fprintf(fid,'PlotFcns: %s\n',options.PlotFcns) ;
                end
                fprintf(fid,'%s\n',['Gauss-Newton method does not work '...
                    'for k=2.']) ;
                fprintf(fid,'%s\n',['Instead, the Nelder-Mead algorithm '...
                    'was used.']) ;
                fclose(fid) ;
            end

        otherwise % l1

            error('opt_soft:SolverError', ...
                'Unknown solver for optimization problem.')

    end


%       complete A to meet constraints (positivity, partition of unity)
    for i = 1:k-1 % l22
        for j = 1:k-1 % l23
            A(i+1,j+1) = alpha(j + (i-1) * (k-1)) ;
        end
    end
    A = fillA(A, evs) ;
    chi = (evs * A) ;

%       close the file with values of the
%       objective function gained during optimization
    fclose(optfile) ;

end
