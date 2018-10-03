function [final_par] = main_nlscon( xguess, par, problem, MAXITER, ...
    XSCALE, XTOL, FILEID, counter )
% initialize user-customized input parameters for NLSCON and call NLSCON
%
% [ final_par ] = main_nlscon( xguess, par, problem, MAXITER, XSCALE, XTOL,
%                              FILEID, counter )
%
% input to nlscon driver routine:
% -------------------------------
% 	xguess      initial parameter guess
% 	par         measurement timepoints and values
% 	problem     file name of model description 
%               (rhs evaluation and jacobian calculation)
%   MAXITER     Number of Gauss-Newton-iterations
%   XSCAL       factor for the scaling vector 
%               xscale = XSCALE * ones(size(xguess))
%   XTOL        relative (!) tolerance in x 
%               (real tolerance related to XSCALE*XTOL)
%   FILEID      file-ID for the naming of the output files
%   counter     gives the actual number of clusters
%
% ouput:
% ------
%   final_par   resulting vector of iteration run
%
%     ____________________________________________________________ 
%
%
%
%*  Written by        L. Weimann 
%*  Purpose           Testexample for code NLSCON
%*  Version           2.3
%*  Revision          June 2004
%*  Latest Change     June 2004
%*  Library           CodeLib
%*  Code              Matlab 6.0
%*  Environment       Standard Matlab 6.0 environment on PC's,
%                     workstations and hosts.
%
%   modified by C. Stoetzel, T. Dierkes, 2011
%
%   modified by Bernhard Reuter, Theoretical Physics II, 
%   University of Kassel, 2017
%     ____________________________________________________________
%

    fidout = fopen(strcat(FILEID,'-','n=',num2str(counter),'-', ...
        'nlscon.m.out'),'w'); 
    fiddat = fopen(strcat(FILEID,'-','n=',num2str(counter),'-', ...
        'nlscon.m.dat'),'w');
    fprintf(1,'%s\n',' monitor: nlscon.m.out , data: nlscon.m.dat') ;
  
    xtol = XTOL ;
  
%       Execution mode: 0=Standard Mode, 1=Stepwise mode
    iopt.mode = 0 ;
  
%       Jacobian: 0=(same as value 3)
%                 1=supplied by user routine JAC
%                 2=computed by numerical differentation (no feedback) 
%                 3=computed by numerical differentation (with feedback)
    iopt.jacgen =3 ;
  
%       A posteriori statistical analysis: 0 = no, 1 = yes
    iopt.qstat = 1 ;
    iopt.mprsta = 1 ;
  
%       Problem classification:
%       1 = linear , 2 = mildly nonlinear,  3 = highly nonlinear, 
%       4 = extremely nonlinear
    iopt.nonlin = 4 ;

%       Bounded damping strategy switch:
%       =0 means currently always IBDAMP = off
%       (but may depend on the settings of other
%       options in future versions)
%       =1 means always IBDAMP = on
%       =2 means always IBDAMP = off
    iopt.ibdamp = 1 ;
  
%       Broyden updates: 0 = inhibit, 1=allow
    iopt.qrank1 = 0 ;
  
%       Set output level on various units
    iopt.mprerr = 3 ;
    iopt.mprmon = 6 ; 
    iopt.mprsol = 2 ; 
    iopt.mprtim = 1 ; 
  
%       Set output units
    iopt.luerr = fidout ;
    iopt.lumon = fidout ;
    iopt.lusol = fiddat ;
    iopt.lutim = fidout ;
  
%       Automatic row scaling of linear system (constraints)
%       and user scaling (measurements):
%       0 = allowed , 1 = inhibited
    iopt.norowscal = 0 ;

%       Minimal allowed damping factor
    wk.fcmin = 1e-14 ; 

%       Damping factor for first Gauss-Newton iteration
%       - overrides option NONLIN, if set
%    wk.fcstart = 1e-6;

%       Number of Gauss-Newton-iterations
    wk.nitmax = MAXITER ;

%       Initial estimate of parameters
    x = xguess(:) ;

%       Data obtained by measurements 
%       (or here the value of val=k-trace(S) aimed at)
    fobs=0 ;
%    fobs = par(:,2) ;
%    tpoints = par(:,1);

%    fprintf(fidout,' Simulated  experimential  data  :%s\n',' ') ;
%    fprintf(fidout,' %10.7f\n',fobs) ;
%    fprintf(fidout,'%s\n\n\n',' ') ;

%       User scaling (lower threshold) of the
%       iteration vector X(N)
    xscal = ones(size(x)) ;
    xscal = XSCALE .* xscal ;

%       User weighting vector of measurements (=1 is suitable here...)
    fscal = ones(size(fobs)) ;

%       initialize the model (objective) function to get m
    f = feval(problem,x,'',par) ;
%       Sum of number of measurement data and equality constraints
    m = length(f) ;
  
    stime = cputime ;
  
    [x,info,wk] = nlscon(m,problem,x,xscal,fobs,fscal,xtol,iopt,par,wk) ;

    %   new extra output to file from here...
    infout = fopen(strcat(FILEID,'-','n=',num2str(counter),'-', ...
        'nlscon.m.info.out'),'w') ;
    fprintf(infout, 'xscal\n') ;
    fprintf(infout, '%.15e\n', info.xscal) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'rtol\n') ;
    fprintf(infout, '%.15e\n', info.rtol) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'ierr = %d\n', info.ierr) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'precision\n') ;
    fprintf(infout, '%.15e\n', info.precision) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'dampingfc\n') ;
    fprintf(infout, '%.15e\n', info.dampingfc) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'natlevel\n') ;
    fprintf(infout, '%.15e\n', info.natlevel) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'simlevel\n') ;
    fprintf(infout, '%.15e\n', info.simlevel) ;
    fprintf(infout, '\n') ;
    fprintf(infout, 'stdlevel\n') ;
    fprintf(infout, '%.15e\n', info.stdlevel) ;
    fclose(infout) ;

    etime = cputime ;
    cptime = etime-stime ;
    if cptime ~= 0
        fprintf(fidout, 'Time used = %9.3f Sec\n\n',cptime) ;
        %fprintf('Time used = %9.3f Sec\n\n',cptime) ;
    end
 
    final_par = x(:) ;
   
    fclose(fidout) ;
    fclose(fiddat) ;
  
end