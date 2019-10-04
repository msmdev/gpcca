function [x,info,wk] = nlscon(m,fcn,x,xscal,fi,fscal,rtol,iopt,par,wk,info)
%     ------------------------------------------------------------
%
%*  Title
%
%     Numerical solution of nonlinear (NL) least squares (S)
%     problems with nonlinear constraints (CON), especially
%     designed for numerically sensitive problems.
%
%*  Written by        U. Nowak, L. Weimann 
%*  Purpose           Solution of highly nonlinear, optionally 
%                     constrained least squares problems
%*  Method            Damped affine invariant Gauss-Newton method
%                     (see references below)
%*  Category          K1b2b. - Nonlinear least squares approxi-
%                     mation with nonlinear constraints
%*  Keywords          Nonlinear least squares problems, 
%                     Gauss-Newton methods
%*  Version           2.3.2
%*  Revision          December 1993
%*  Latest Change     August 2001
%*  Library           CodeLib
%*  Code              Matlab 6.0
%*  Environment       Standard Matlab 6.0 environment on PC's,
%                     workstations and hosts.
%*  Copyright     (c) Konrad-Zuse-Zentrum fuer
%                     Informationstechnik Berlin (ZIB)
%                     Takustrasse 7, D-14195 Berlin-Dahlem
%                     phone : + 49/30/84185-0
%                     fax   : + 49/30/84185-125
%*  Contact           Lutz Weimann
%                     ZIB, Division Scientific Computing, 
%                          Department Numerical Analysis and Modelling
%                     phone : + 49/30/84185-185
%                     fax   : + 49/30/84185-107
%                     e-mail: weimann@zib.de
%
%*    References:
%
%     /1/ P. Deuflhard:
%         Newton Methods for Nonlinear Problems. -
%         Affine Invariance and Adaptive Algorithms.
%         Series Computational Mathematics 35, Springer (2004)
%
%     /2/ U. Nowak, L. Weimann:
%         A Family of Newton Codes for Systems of Highly Nonlinear
%         Equations - Algorithm, Implementation, Application.
%         ZIB, Technical Report TR 90-10 (December 1990)
%
%  ---------------------------------------------------------------
%
%* Licence
%    You may use or modify this code for your own non commercial
%    purposes for an unlimited time. 
%    In any case you should not deliver this code without a special 
%    permission of ZIB.
%    In case you intend to use the code commercially, we oblige you
%    to sign an according licence agreement with ZIB.
%
%* Warranty 
%    This code has been tested up to a certain level. Defects and
%    weaknesses, which may be included in the code, do not establish
%    any warranties by ZIB. ZIB does not take over any liabilities
%    which may follow from acquisition or application of this code.
%
%* Software status 
%    This code is under care of ZIB and belongs to ZIB software class 1.
%
%     ------------------------------------------------------------
%
%*    Summary:
%     ========
%     Damped Gauss-Newton-algorithm with rank strategy for highly 
%     nonlinear least squares approximation problems (optionally
%     constrained, over- and underdetermined problems) - 
%     due to Ref.(1).
%
%     (The iteration is done by function NCINT currently. NLSCON
%      itself does some house keeping and builds up workspace.)
%
%     The problem solved by this program looks as follows:
%
%     Denoting below the n-dimensional real space with IR(n),
%     the number of parameters to be estimated with N,
%     the number of measurements (data to be fitted) with MFIT, and
%     the number of equality constraints with MCON,
%     M := MCON + MFIT ,
%     let   F : IR(N) --> IR(MFIT) ,  G : IR(N) --> IR(MCON) 
%     be nonlinear functions,
%     FI     in IR(MFIT) the vector of measurement data and 
%     FSCAL  in IR(MFIT) the vector of measurement weights.
%
%     For M >= N, find a parameter vector X in IR(N), which 
%     minimizes  Sum (j=1,...,MFIT) ((Fj(X)-FIj)/FSCALj)**2 and
%     satisfies   G(X) = 0.
%     For M < N, find the parameter vector X in IR(N) with the 
%     smallest possible euclidian norm, 
%     which satisfies  F(X) = 0  and  G(X) = 0  . 
%
%     Jacobian approximation by numerical differences or user
%     supplied function FCN, called with 'jacobian'-flag.
%
%     The numerical solution of the arising linear least squares
%     problem is done by means of the functions DECCON and SOLCON
%     (QR decomposition with subcondition estimation, rank decision
%     and computation of the rank-deficient pseudoinverse) .
%     For special purposes these routines may be substituted.
%
%     A statistical a posteriori analysis of the parameter estimate
%     is optionally available.
%
%     This is a driver routine for the core solver NCINT.
%
%     ------------------------------------------------------------
%
%*    Parameters list description (* marks inout parameters)
%     ======================================================
%
%*    External functions (to be supplied by the user)
%     =================================================
% 
%     [FOUT,FAIL] =FCN(X,FLAG,PAR) String  Name of problem- and Jacobian function
%       X()           Float   Vector of unknowns (input)
%       FLAG          String  Operation flag - the following values must be
%                             supported:
%                             '' (empty string) : FOUT must return the problem-
%                                                 function value f(X) ;
%                             'jacobian' : FOUT must return the associated 
%                                          Jacobian Jac(x)
%       PAR            AnyType  A (required!) user parameter
%
%       FOUT()         Float  Vector of returned function values or Jacobian matrix
%       FAIL           Int    evaluation-failure indicator. (output)
%                             On output: Indicates failure of FCN eval-
%                                uation, if having a value <= 2.
%                             If <0: NLEQ1 will be terminated with 
%                                    error code = 82, and FAIL stored
%                                    to wk.ifail.
%                             If =1: A new trial Newton iterate will
%                                    computed, with the damping factor
%                                    reduced to its half.
%                             If =2: A new trial Newton iterate will
%                                    computed, with the damping factor
%                                    reduced by a reduct. factor, which
%                                    must be output through F(1) by FCN,
%                                    and it is value must be >0 and < 1.
%                             Note, that if FAIL = 1 or 2, additional
%                             conditions concerning the damping factor,
%                             e.g. the minimum damping factor or the
%                             bounded damping strategy may also influ-
%                             ence the value of the reduced damping 
%                             factor.
%
%
%*    Input parameters of NLSCON
%     ==========================
%
%     M              Int    Sum of number of measurement data and
%                           equality constraints
%     X(N)           Dble   Initial estimate of parameters
%     XSCAL(N)       Dble   User scaling (lower threshold) of the 
%                           iteration vector X(N)
%     FI(MFIT)       Dble   Data obtained by measurements
%     FSCAL(MFIT)    Dble   User weighting vector of measurements
%     RTOL           Dble   Required relative precision of
%                           solution components -
%                           RTOL >= EPMACH*TEN*N
%     IOPT(50)       Int    Structure of run-time options. Set the 
%                           fields to zero or omit them to get
%                           default values (details see below)
%     PAR            AnyTyp This parameter will be passed on every
%                           call of the user function FCN as the
%                           third input argument. 
%     WK             Struct A structured workspace. Normally an empty
%                           array should be passed to this argument.
%                           However, on successor steps in onestep mode,
%                           the output parameter WK from the previous step
%                           call must be passed to this parameter.
%
%*    Output parameters of NLSCON
%     ===========================
%
%     X(N)           Dble   Solution parameters ( or final values,
%                           respectively )
%     INFO           Struct A structure, holding additional info
%                           obtained from the iteration run. The
%                           fields are as follows:
%       XSCAL(N)     Dble   After return with IERR >= 0, it contains
%                           the latest internal scaling vector used
%                           After return with IERR == -1 in onestep-
%                           mode it contains a possibly adapted 
%                           (as described below) user scaling vector:
%                           if XSCAL(I) <  SMALL) XSCAL(I) = SMALL ,
%                           if XSCAL(I) >  GREAT) XSCAL(I) = GREAT .
%                           For SMALL and GREAT, see section machine
%                           constants below  and regard note 1.
%       RTOL         Dble   Finally achieved (relative) accuracy
%                           The estimated absolute error of component i
%                           of x_out is approximately given by
%                             abs_err(i) = RTOL * XSCAL_out(i) ,
%                           where (approximately)
%                             XSCAL_out(i) = 
%                                max(abs(X_out(i)),XSCAL_in(i)).
%                           Note that RTOL_out may be greater than
%                           RTOL_in, but NLSCON claims 'solution found'
%                           - see IOPT(36).
%       IERR         Int    Return value parameter
%                           =-1 sucessfull completion of one iteration
%                               step, subsequent iterations are needed 
%                               to get a solution. (stepwise mode only)
%                           = 0 successfull completion of iteration
%                           > 0 see list of error messages below
%        XITER(N,:)   Dble   An array holding all iterates of the
%                            Newton iteration run, stored as column
%                            vectors. 
%        NATLEVEL(:)  Dble   The sequence of natural levels of the 
%                            newton corrections over the iteration steps
%        SIMLEVEL(:)  Dble   The sequence of natural levels of the simplified 
%                            newton corrections over the iteration steps 
%        STDLEVEL(:)  Dble   The sequence of standard levels over the
%                            iteration steps
%        PRECISION(:) Dble   The sequence of acheived precisions over
%                            the iteration steps.
%        DAMPINGFC(:) Dble   The sequence of accepted damping factors
%                            over the iteration steps.
%        NITER        Int    Number of Newton-iterations done
%        NCORR        Int    Number of corrector steps done
%        NREJR1       Int    Number of rejected Newton iteration steps
%                             done with a rank-1 approximated Jacobian
%        NJAC         Int    Number of Jacobian generations or JAC-calls done
%        NFCN         Int    Number of FCN-evaluations done
%        NFCNJ        Int    Number of FCN-evaluations for Jacobian
%                             approximation done
%
%     Note 1.
%        The machine dependent values SMALL, GREAT and EPMACH are
%        gained from calls of the machine constants function D1MACH.
%        As delivered, this function is adapted to use constants 
%        suitable for all machines with IEEE arithmetic. If you use
%        another type of machine, you have to change the DATA state-
%        ments for IEEE arithmetic in D1MACH into comments and to 
%        uncomment the set of DATA statements suitable for your machine.
%
%*   Options IOPT.name:
%    ==================
%
%          Name   Default  Meaning
%
%          mode   0        =0 Standard mode initial call:
%                             Return when the required accuracy for the
%                             iteration vector is reached. User defined
%                             parameters are evaluated and checked.
%                             Standard mode successive call:
%                             If NLSCON was called previously with
%                             MODE=1, it performs all remaining 
%                             iteration steps.
%                          =1 Stepwise mode:
%                             Return after one Gauss-Newton
%                             iteration step.
%          jacgen 0        Method of Jacobian generation
%                          =0 Standard method is JACGEN=2
%                          =1 User supplied function FCN will be 
%                             called to generate Jacobian matrix
%                          =2 Jacobian approximation by numerical
%                             differentation (see function NCJAC)
%                          =3 Jacobian approximation by numerical
%                             differentation with feedback control
%                             (see function NCJCF)
%          iscal  0        Determines how to scale the iterate-vector:
%                          =0 The user supplied scaling vector XSCAL is
%                             used as a (componentwise) lower threshold
%                             of the current scaling vector
%                          =1 The vector XSCAL is always used as the
%                             current scaling vector
%          mprerr 0        Print error messages
%                          =0 No output
%                          =1 Error messages
%                          =2 Warnings additionally
%                          =3 Informal messages additionally
%          luerr  1        Logical unit number for error messages
%          mprmon 0        Print iteration Monitor
%                          =0 No output
%                          =1 Standard output
%                          =2 Summary iteration monitor additionally
%                          =3 Detailed iteration monitor additionally
%                          =4,5,6 Outputs with increasing level addi-
%                             tional increasing information for code
%                             testing purposes. Level 6 produces
%                             in general extremely large output!
%          lumon  1        Logical unit number for iteration monitor
%          mprsol 0        Print solutions
%                          =0 No output
%                          =1 Initial values and solution values
%                          =2 Intermediate iterates additionally
%          lusol  1        Logical unit number for solutions
%          mprtim 0        Output level for the time monitor
%                          = 0 : no time measurement and no output
%                          = 1 : time measurement will be done and
%                                summary output will be written -
%                                regard note 4a.
%          lutim  1        Logical output unit for time monitor
%          qstat  0        Statistical Analysis of the final 
%                          least squares estimate:
%                          = 0 : Analysis will not be done
%                          = 1 : Analysis will be done,
%                                and certain results are stored
%                                to the RWK array (for details, see
%                                RWK description below)
%          mprsta 0        Printing of statistical Analysis for
%                          the final least squares estimate:
%                          = 0 : Printing will not be done
%                          = 1 : Printing will be done (and this
%                                implies QSTAT to be set to 1)
%          nonlin 3        Problem type specification
%                          =1 Linear problem
%                             Warning: If specified, no check will be
%                             done, if the problem is really linear, and
%                             NLSCON terminates unconditionally after
%                             one Gauss-Newton-iteration step.
%                          =2 Mildly nonlinear problem
%                          =3 Highly nonlinear problem
%                          =4 Extremely nonlinear problem
%          qrank1 0        =0 (0) Rank-1 updates by Broyden-
%                             approximation are inhibited.
%                          =1 (1) Rank-1 updates by Broyden-
%                             approximation are allowed.
%          qnscal 0        Inhibit automatic row scaling: 
%                          =0 (0) Automatic row scaling of
%                             the linear system is activ: 
%                             Rows i=1,...,N will be divided by
%                             max j=1,...,N (abs(a(i,j))) 
%                          =1 (1) No row scaling of the linear
%                             system. Recommended only for well row-
%                             scaled nonlinear least squares problems.
%          iterm  0        Determines the iteration termination cri-
%                          terium to be chosen:
%                          =0 Iteration is terminated, if one of the
%                             stopping criteria used for ITERM=1 and
%                             ITERM=2 is satisfied (see below).
%                          =1 Iteration is terminated, if, from a
%                             statistical point of view, a resonable
%                             precision is achieved, i.e.
%                             simplified GN-correction < RTOL, and an
%                             estimate of the accuracy is available.
%                             Recommended to be used for incompatible
%                             problems.
%                          =2 Iteration is terminated, if 
%                             GN-correction < RTOL . Using this option
%                             may force 'to much' precision from
%                             the statistical point of view.
%          ibdamp          Bounded damping strategy switch:
%                          =0 means currently always IBDAMP = off 
%                             (but may depend on the settings of other
%                              options in future versions)
%                          =1 means always IBDAMP = on 
%                          =2 means always IBDAMP = off 

%     Note 3:
%         If NLSCON terminates with IERR=2 (maximum iterations)
%         or  IERR=3 (small damping factor), you may try to continue
%         the iteration by increasing NITMAX or decreasing FCMIN
%         (see WK) and setting WK.QSUCC to 1.
%
%*   Optional input/output in WK:
%    ============================
%
%          Name          Meaning
%
%          niter  IN/OUT Number of Gauss-Newton-iterations
%          ncorr  IN/OUT Number of corrector steps
%          nfcn   IN/OUT Number of FCN-evaluations
%          njac   IN/OUT Number of Jacobian generations or
%                        JAC-calls
%          nfcnj  IN/OUT Number of FCN-evaluations for Jacobian
%                        approximation
%          nrejr1 IN/OUT Number of rejected Gauss-Newton iteration
%                        steps done with a rank-1 computed Jacobian
%          idcode IN/OUT Output: The 8 decimal digits program identi-
%                        fication number ppppvvvv, consisting of the
%                        program code pppp and the version code vvvv.
%                        Input: If containing a negative number,
%                        it will only be overwritten by the identi-
%                        fication number, immediately followed by
%                        a return to the calling program.      
%          ifail  OUT    Set in case of failure of NCFACT (IERR=80),
%                        N2SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
%                        to the nonzero IFAIL value returned by the 
%                        routine indicating the failure .
%          nitmax IN     Maximum number of permitted iteration
%                        steps (default: 50)
%          irank  IN     Initial rank
%                        =0 : full rank N
%                        =k with 0 < k < N : deficient rank assumed
%                           for the Jacobian in the starting point
%          new    IN/OUT Count of consecutive rank-1 updates
%          ifccnt INTERN Count of consecutive special steps before 
%                        convergence stop of iteration
%          conv   OUT    The maximum norm of the latest ordinary 
%                        respective simplified (scaled) Gauss-Newton
%                        correction.
%          sumx   OUT    Natural level (not Normx of printouts)
%                        of the current iterate, i.e. Sum(DX(i)**2),
%                        where DX = scaled Gauss-Newton correction.
%          dlevf  OUT    Standard level (not Normf of printouts)
%                        of the current iterate, i.e. Norm2(F(X)),
%                        where F =  nonlinear model function.
%          fcbnd  IN     Bounded damping strategy restriction factor
%                        (Default is 10)
%          fcstrt IN     Damping factor for first Gauss-Newton iteration
%                        - overrides option NONLIN, if set (see note 5)
%          fcmin  IN     Minimal allowed damping factor (see note 5)
%          sigma  IN     Broyden-approximation decision parameter
%                        Required choice: SIGMA >= 1. Increasing this
%                        parameter make it less probable that the algo-
%                        rith performs Broyden steps.
%                        Rank1 updates are inhibited, if 
%                        SIGMA > 1/FCMIN is set. (see note 5)
%          cond   IN     Maximum permitted subcondition for rank-
%                        decision by linear solver
%                        (Default: 1/epmach, epmach: relative
%                         machine precision) 
%          ajdel  IN     Jacobian approximation without feedback:
%                        Relative pertubation for components
%                        (Default: sqrt(epmach*10), epmach: relative
%                         machine precision) 
%          ajmin  IN     Jacobian approximation without feedback:
%                        Threshold value (Default: 0.0)
%                          The absolute pertubation for component k is
%                          computed by 
%                          DELX := AJDEL*max(abs(Xk),AJMIN)
%         etadif  IN     Jacobian approximation with feedback:
%                        Target value for relative pertubation ETA of X
%                        (Default: 1.0e-6)
%         etaini  IN     Jacobian approximation with feedback:
%                        Initial value for denominator differences
%                        (Default: 1.0e-6)
%         prec    OUT    An estimate for the achieved relative accuracy.
%                        This number is only available, if IERR=0 or 
%                        IERR=1 and an estimate for the incompatibility
%                        factor kappa (SKAP=RWK(32)) can be made. If no
%                        meaningful information is at hand, a -1.0 is
%                        stored.
%         skap    OUT    An estimate of the incompatibility factor of
%                        the given nonlinear least squares problem.
%                        This number is only available, if IERR=0 or 
%                        IERR=1 and a certain further condition holds.
%                        If no meaningful information is at hand, a -1.0
%                        is stored.
%                        A small value of SKAP indicates that the given
%                        measurement data can be well fitted by the
%                        model, whereas a value SKAP >= 0.5 may give a 
%                        hint to an insufficient modelling of the
%                        experiment from which the measurements
%                        originate.
%         sigma2  OUT    Holds the estimated variance of the residual
%                        on final exit of NLSCON, if IOPT(21)=1 is set.
%         xl(1:n) OUT    Holds the left bounds of the final parameters
%                        confidence intervals; 
%                        on final exit of NLSCON, if IOPT(21)=1 is set.
%         xr(1:n) OUT    Holds the right bounds of the final parameters
%                        confidence intervals; 
%                        on final exit of NLSCON, if IOPT(21)=1 is set.
%         vcv(1:n,1:n)
%                 OUT    Holds the columnwise stored correlation-matrix
%                        on final exit of NLSCON, if IOPT(21)=1 is set.
%
%     Note 5:
%       The default values of the internal parameters may be obtained
%       from the monitor output with at least IOPT field MPRMON set to 2.
%
%*   Error messages:
%    ===============
%
%      1    Termination at stationary point (rank deficient Jacobian
%           and termination criterion fulfilled)
%      2    Termination after NITMAX iterations ( as indicated by
%           input parameter NITMAX=IWK(31) )
%      3    Termination, since damping factor became to small and
%           Jacobian rank was already reduced to zero
%     20    Bad or inconsistent input to the dimensional parameter M
%     21    Nonpositive value for RTOL supplied
%     22    Negative scaling value via vector XSCAL supplied
%     30    One or more fields specified in IOPT are invalid
%           (for more information, see error-printout)
%     80    Error signalled by linear solver routine NCFACT,
%           for more detailed information see IFAIL-value
%           stored to IWK(23)
%           (not used with standard routine NCFACT)
%     81    Error signalled by linear solver routine NCFIT,
%           for more detailed information see IFAIL-value
%           stored to IWK(23)
%           (not used with standard routine NCFIT)
%     82    Error signalled by user routine FCN when called
%           with ''-flag (blank-flag) (Nonzero value
%           returned via IFAIL-flag; stored to WK.ifail )
%     83    Error signalled by user routine FCN when called 
%           with 'jacobian'-flag (Nonzero value
%           returned via IFAIL-flag; stored to WK.ifail )
%     180,182,183
%           see error codes 80,82,83, but the failure of NCFACT or 
%           FCN occured when preparing the call of the statistics
%           function STACON for the final iterate of NLSCON.
%
%     Note 6 : in case of failure:
%        -    use non-standard options
%        -    use another initial guess
%        -    or reformulate model
%        -    or turn to general optimization routine
%
%*    Machine dependent constants used:
%     =================================
%
%     DOUBLE PRECISION EPMACH  in  NLSCON, NCPCHK, NCINT
%     DOUBLE PRECISION GREAT   in  NLSCON, NCPCHK
%     DOUBLE PRECISION SMALL   in  NLSCON, NCPCHK, NCINT, NCSCAL
%
%*    Functions called: NCPCHK, NCINT, D1MACH
%
%     ------------------------------------------------------------
%*    End Prologue
%
%*    Summary of changes of the underlying Fortran code:
%     ==================================================
%      
%     2.2.1   91, June  3   Time monitor included
%     2.2.2   91, June  3   Bounded damping strategy implemented
%     2.2.3   91, July 26   AJDEL, AJMIN as RWK-options for JACGEN == 2,
%                           ETADIF, ETAINI as RWK-opts. for JACGEN == 3
%                           FCN-count changed for anal. Jacobian
%             91, Sept.     DECCON with new fail exit, for the case that
%                           the square root of a negative number will
%                           appear during pseudo inverse computation.
%                           (Occured, although theoretical impossible!)
%     2.2.6  91, Sept.  17  Damping factor reduction by FCN-fail imple-
%                           mented
%                Sept.  24  Meaning of option ITERM modified.
%     2.3    91, Dec.   20  New Release for CodeLib
%            92, June    2  Level of accepted simplified correction
%                           stored to RWK(IRWKI+4)
%     2.3.1  92, Oct.   13  Corrected errors in subr. NCJAC, NCJCF:
%                           dimensions of FX(M), FU(M),
%                           description of FCN in these functions
%     2.3.2  93, Nov.   24  Optional call of statistical analysis
%                           routine STACON to get statistics about
%                           the quality of the parameter estimate
%                           which has been computed by NLSCON.
%            00, July   12  RTOL output-value bug fixed
%            05, Sept.  12  FCREDU setting bug near 3.6.1 fixed
%     2.3.3  06, Nov.   23  Spurious number output fixed, and made
%                           several changes recommended by M-Lint
%                           tool of MATLAB version 7.2.0.
% 
%     ------------------------------------------------------------
%
%     Further WK positions (only internally used)
%
%               Name     Meaning
%
%               FCA      Previous damping factor
%               SUMXS    natural level of accepted simplified correction
%
%     Internal arrays stored in WK (see routine NCINT for descriptions)
%
%               Array         Type   Remarks
%
%               AA(M,N)       Perm
%               DX(N)         Perm  
%               DXQ(N)        Perm 
%               XA(N)         Perm
%               FMODEL(M)     Perm
%               F(M)          Perm
%               FA(M)         Perm
%               FW(M)         Perm
%               ETA(N)        Perm   Only used for JACGEN=3
%               A(M,N)        Temp
%               QA(N,N)       Temp
%               XW(N)         Temp
%               DXQA(N)       Temp
%               QU(M)         Temp
%               RQ(M)         Temp
%               RQKEEP(M)     Temp
%               DELXQ(N)      Temp
%     
%      
      persistent xiter sumxall dlevfall sumxqall tolall fcall ;
      one  =  1.0 ;
      ten  = 10.0 ;
      zero =  0.0 ;
      iver = 22112323 ;
%
%     Version: 2.3.3             Latest change:
%     -----------------------------------------
%
      chgdat  = 'November 23, 2006   ' ;
      prodct  =  'NLSCON  '            ;
%*    Begin
      n = length(x) ;
      mfit = length(fi) ;
      epmach = d1mach(3) ;
      great  = sqrt(d1mach(2)/ten) ;
      small  = d1mach(6) ;
      tolmin = epmach*ten*n ;
      ierr = 0 ;
      if isfield(wk,'iver')
        qvchk = wk.iver < 0 ;
        wk.iver = iver ;
        if qvchk
          huch=qvchk
          return ;
        end ;
      else
        wk.iver = iver ;
      end

%        Print error messages?
      [mprerr,iopt] = getopt(iopt,'mprerr',0) ;
      [luerr, iopt] = getopt(iopt,'luerr', 1) ;
      if luerr  ==  0
        luerr = 1 ;
        iopt.luerr = luerr ;
      end
%        print iteration monitor?
      [mprmon,iopt] = getopt(iopt,'mprmon',0) ;
      [lumon ,iopt] = getopt(iopt,'lumon', 1) ;
      if lumon  <=  0  ||  lumon  >  99
        lumon = 1 ;
        iopt.lumon = lumon ;
      end
%        print intermediate solutions?
      [mprsol,iopt] = getopt(iopt,'mprsol',0) ;
      [lusol, iopt] = getopt(iopt,'lusol', 1) ;
      if lusol  ==  0
        lusol = 1 ;
        iopt.lusol=lusol ;
      end
%        print time summary statistics?
      [mprtim,iopt] = getopt(iopt,'mprtim',0) ;
      [lutim, iopt] = getopt(iopt,'lutim', 1) ;
      if lutim  ==  0
        lutim = 1 ;
        iopt.lutim=lutim ;
      end
      [qsucc,wk] = getopt(wk,'qsucc',0) ;
      qinimo = mprmon >= 1 &&  ~ qsucc ;
%     Print NLSCON heading lines
      if qinimo
        fprintf(lumon,'%s\n\n%s\n\n','    N L S C O N  *****  V e r s i o n  2 . 3 . 2 ***', ...
        ' Gauss-Newton-Method for the solution of nonlinear least squares problems') ;
      end
%     Check input parameters and options
      [ierr,iopt,rtol,xscal,fscal] = ncpchk(n,m,mfit,x,xscal,fi,fscal,rtol,iopt,wk) ;
%     Exit, if any parameter error was detected till here
      if ierr ~= 0
          huch1=ierr
        return ;
      end
      if ~ qsucc
        xiter    = [] ;
        sumxall  = [] ;
        dlevfall = [] ;
        sumxqall = [] ;
        tolall   = [] ;
        fcall    = [] ;
      end
      mcon = m-mfit ;
      m1 = m ;
      m2 = m ;
      maxmn = max(m,n) ;
      [jacgen,iopt] = getopt(iopt,'jacgen',0) ;
      if jacgen == 0
        jacgen=2 ;
      end
      iopt.jacgen = jacgen ;
%     WorkSpace: WK
      wk = iniopt(wk,'xl',    zeros(n,1)) ;
      wk = iniopt(wk,'xr',    zeros(n,1)) ;
%      wk = iniopt(wk,'vcv',   zeros(n,n)) ;
%      wk = iniopt(wk,'aa',    zeros(m2,n)) ;
%      wk = iniopt(wk,'a',     zeros(m1,n)) ;
%      wk = iniopt(wk,'qa',    zeros(n,n)) ;
      wk = iniopt(wk,'vcv',   0.0) ;
      wk = iniopt(wk,'aa',    0.0) ;
      wk = iniopt(wk,'a',     0.0) ;
      wk = iniopt(wk,'qa',    0.0) ;
      wk = iniopt(wk,'dx',    zeros(n,1)) ;
      wk = iniopt(wk,'dxq',   zeros(n,1)) ;
      wk = iniopt(wk,'xa',    zeros(n,1)) ;
      wk = iniopt(wk,'xwa',   zeros(n,1)) ;
      wk = iniopt(wk,'fmodel',zeros(m,1)) ;
      wk = iniopt(wk,'f',     zeros(m,1)) ;
      wk = iniopt(wk,'fa',    zeros(m,1)) ;
      wk = iniopt(wk,'eta',   zeros(n,1)) ;
      wk = iniopt(wk,'xw',    zeros(n,1)) ;
      wk = iniopt(wk,'fw',    zeros(m,1)) ;
      wk = iniopt(wk,'dxqa',  zeros(n,1)) ;
      wk = iniopt(wk,'qu',    zeros(m,1)) ;
      wk = iniopt(wk,'rq',    zeros(m,1)) ;
      wk = iniopt(wk,'rqkeep',zeros(m,1)) ;
      wk = iniopt(wk,'delxq', zeros(n,1)) ;

      wk = iniopt(wk,'fcmin',0.0) ;
      wk = iniopt(wk,'sigma',0.0) ;
      wk = iniopt(wk,'fca',0.0) ;
      wk = iniopt(wk,'conv',0.0) ;
      wk = iniopt(wk,'sumx',0.0) ;
      wk = iniopt(wk,'sumxs',0.0) ;
      wk = iniopt(wk,'dlevf',0.0) ;
      wk = iniopt(wk,'niter',0) ;
      wk = iniopt(wk,'ncorr',0) ;
      wk = iniopt(wk,'nfcn',0) ;
      wk = iniopt(wk,'njac',0) ;
      wk = iniopt(wk,'nfcnj',0) ;
      wk = iniopt(wk,'nrejr1',0) ;
      wk = iniopt(wk,'new',0) ;
      wk = iniopt(wk,'iconv',0) ;
      wk = iniopt(wk,'ifccnt',0) ;
%
      iopt = iniopt(iopt,'norowscal',0) ;
      if qinimo
        fprintf(lumon,'\n %s : %4i\n %s : %4i\n %s : %4i\n\n %s : %10.2e\n', ...
                      'Number of parameters to be estimated (N)',n, ...
                      'Number of data to fitted, e.g. observations (MFIT)',mfit, ...
                      'Number of equality constraints (MCON) : ',mcon, ...
                      'Prescribed relative precision (RTOL) : ',rtol) ;
        if jacgen == 1
          jacg = 'a user function' ;
        elseif jacgen == 2
          jacg = 'numerical differentiation (without feedback strategy)' ;
        elseif jacgen == 3
          jacg = 'numerical differentiation (feedback strategy included)' ;
        end
        fprintf(lumon,'\n The Jacobian is supplied by %s\n',jacg) ;
        if iopt.norowscal == 1
          rsmode = 'inhibited' ;
        else
          rsmode = 'allowed' ;
        end
        fprintf(lumon,' Automatic row scaling of the jacobian is %s\n',rsmode) ;
      end
      [qrank1,iopt] = getopt(iopt,'qrank1',0) ;
      [nonlin,iopt]=getopt(iopt,'nonlin',3) ;
      iopt = iniopt(iopt,'boundeddamp',0) ;
      if iopt.boundeddamp == 0
        qbdamp = nonlin == 4 ;
      elseif iopt.boundeddamp == 1
        qbdamp = 1 ;
      elseif iopt.boundeddamp == 2
        qbdamp = 0 ;
      end
      wk = iniopt(wk,'fcbnd',0.0) ;
      if qbdamp
        if wk.fcbnd < 1.0
           wk.fcbnd = 10.0 ;
        end
      end
      if qrank1 && m > n
        if mprerr >= 2
          fprintf(luerr,'\n\n Warning: Broyden-steps set to inhibited for the overdetermined system') ;
        end ;
        qrank1 = 0 ;
        iopt.qrank1 = 0 ;
      end
      if qinimo
        if qrank1
          rk1mode = 'allowed' ;
        else
          rk1mode = 'inhibited' ;
        end
        fprintf(lumon,'\n Rank-1 updates are %s\n',rk1mode) ;
        if nonlin == 1
          nonlinmode = 'linear' ;
        elseif nonlin == 2
          nonlinmode = 'mildly nonlinear' ;
        elseif nonlin == 3
          nonlinmode = 'highly nonlinear' ;
        elseif nonlin == 4
          nonlinmode = 'extremely nonlinear' ;
        end
        fprintf(lumon,' Problem is specified as being %s\n',nonlinmode) ;
        if qbdamp
          fprintf(lumon, ...
                  ' Bounded damping strategy is active\n Bounding factor is %10.3e\n', wk.fcbnd) ;
        else
          fprintf(lumon, ' Bounded damping strategy is off\n') ;
        end
      end
%     Maximum permitted number of iteration steps
      [nitmax,wk]=getopt(wk,'nitmax',50) ;
      if nitmax <= 0
        nitmax=50 ;
      end
      wk.nitmax=nitmax ;
      if qinimo
        fprintf(lumon,' Maximum permitted number of iteration steps : %6i\n',nitmax) ;
      end
%     Initial damping factor for highly nonlinear problems
      wk = iniopt(wk,'fcstart',0.0) ;
      qfcstr = wk.fcstart > 0.0 ;
      if  ~ qfcstr
        wk.fcstart=1.0e-2 ;
        if nonlin == 4
          wk.fcstart=1.0e-4 ;
        end
      end
%     Minimal permitted damping factor
      wk = iniopt(wk,'fcmin',0.0) ;
      if wk.fcmin <= 0.0
        wk.fcmin=1.0e-4 ;
        if nonlin == 4
          wk.fcmin=1.0e-8 ;
        end
      end
      fcmin=wk.fcmin ;
%     Broyden-update decision parameter SIGMA
      wk = iniopt(wk,'sigma',0.0) ;
      if wk.sigma < 1.0
        wk.sigma=2.0 ;
      end
      if  ~ qrank1
        wk.sigma=10.0/fcmin ;
      end
%     Starting value of damping factor (FCMIN <= FC <= 1.0)
      if nonlin <= 2 &&  ~ qfcstr
%       for linear or mildly nonlinear problems
        fc = 1.0 ;
      else
%       for highly or extremely nonlinear problems
        fc = wk.fcstart ;
      end
      wk.fcstart = fc ;
%     Initial rank
      wk = iniopt(wk,'irank',0) ;
      irank = wk.irank ;
      minmn = min(m,n) ;
      if irank <= 0 || irank > minmn
        wk.irank = minmn ;
      end
%     Maximum permitted sub condition number of matrix A
      [cond,wk] = getopt(wk,'cond',0.0) ;
      if cond < one
        cond = one/epmach ;
      end
      wk.cond = cond ;
      if mprmon >= 2 && ~ qsucc
        fprintf(lumon,'\n\n%s\n\n%s%9.2e\n%s%9.2e\n%s%9.2e\n%s%6i\n%s%9.2e\n',    ...
                      ' Internal parameters:',                                    ...
                      ' Starting value for damping factor FCSTART = ',wk.fcstart, ...
                      ' Minimum allowed damping factor FCMIN = ',fcmin,           ...
                      ' Rank-1 updates decision parameter SIGMA = ',wk.sigma,     ...
                      ' Initial Jacobian pseudo-rank IRANK =',wk.irank,           ...
                      ' Maximum permitted subcondition COND = ',cond) ;
      end
      if  ~ qsucc
         for i=1:mfit
           if fscal(i) >= small && fscal(i) <= great
             wk.fw(mcon+i) = one/fscal(i) ;
           else
             wk.fw(mcon+i) = one ;
             if fscal(i) ~= zero && mprerr >= 2
               fprintf(luerr,'\n Warning: Bad scaling value fscal(%5i) = %10.2e replaced by 1.0\n', ...
                             i, fscal(i) ) ;
             end
           end
         end
      end
%
%       Initialize and start time measurements monitor
%
      dummy = 1 ;
      if  (~ qsucc)  &&  mprtim ~= 0 
        monini(' NLSCON',lutim) ;
        mondef(0,'NLSCON') ;
        mondef(1,'FCN') ;
        mondef(2,'Jacobi') ;
        mondef(3,'Lin-Fact') ;
        mondef(4,'Lin-Sol') ;
        mondef(5,'Output') ;
        monsrt(dummy) ;
      end
      
%
      ierr = -1 ;
%     If IERR is unmodified on exit, successive steps are required
%     to complete the Gauss-Newton iteration


    [x,xscal,fi,rtol,wk.irank,ierr,wk.xl,wk.xr,wk.vcv,                         ...
        wk.aa,wk.a,wk.qa,wk.dx,wk.dxq,wk.xa,wk.fmodel,                           ...
        wk.f,wk.fa,wk.eta,wk.fw,wk.xw,wk.dxqa,wk.qu,                             ...
        wk.rq,wk.rqkeep,wk.delxq,wk.fcstart,wk.fcmin,wk.sigma,wk.fca,cond,       ...
        wk.conv,wk.sumx,wk.sumxs,wk.dlevf,                                       ...
        xiter,sumxall,dlevfall,sumxqall,tolall,fcall,                            ...
        wk.niter,wk.ncorr,wk.nfcn,wk.njac,wk.nfcnj,wk.nrejr1,wk.new,wk.ifccnt,   ...
        wk.qsucc] =                                                              ...
      ncint(n,m,mfit,fcn,par,x,xscal,fi,rtol,nitmax,nonlin,                      ...
        wk.irank,iopt,ierr,wk,m1,m2,maxmn,                                       ...
        wk.aa,wk.a,wk.qa,wk.dx,wk.dxq,wk.xa,wk.fmodel,                           ...
        wk.f,wk.fa,wk.eta,wk.fw,wk.xw,wk.dxqa,wk.qu,                             ...
        wk.rq,wk.rqkeep,wk.delxq,wk.fcstart,wk.fcmin,wk.sigma,wk.fca,cond,       ...
        wk.conv,wk.sumx,wk.sumxs,wk.dlevf,                                       ...
        xiter,sumxall,dlevfall,sumxqall,tolall,fcall,                            ...
        tolmin,mprerr,mprmon,mprsol,luerr,lumon,lusol,wk.niter,                  ...
        wk.ncorr,wk.nfcn,wk.njac,wk.nfcnj,wk.nrejr1,wk.new,wk.ifccnt,qbdamp) ;
%
      if mprtim ~= 0 && ierr ~= -1 && ierr ~= 10
        monend(dummy) ;
      end
%     Print statistics
      if mprmon >= 1 && ierr ~= -1 && ierr ~= 10
        fprintf(lumon,'\n%s%s%s\n%s%7i%s\n%s%7i%s\n%s%7i%s\n%s%7i%s\n%s%7i%s\n%s%7i%s\n%s\n\n', ...
                '   ******  Statistics * ',prodct, ' *******', ...
                '   ***  Gauss-Newton iter.: ',wk.niter,  '  ***', ...
                '   ***  Corrector steps   : ',wk.ncorr,  '  ***', ...
                '   ***  Rejected rk-1 st. : ',wk.nrejr1, '  ***', ...
                '   ***  Jacobian eval.    : ',wk.njac,   '  ***',  ...
                '   ***  Function eval.    : ',wk.nfcn,   '  ***', ...
                '   ***  ...  for Jacobian : ',wk.nfcnj,  '  ***', ...
                '   *************************************') ;
      end
      info.xscal = xscal ;
      if ierr == -1
        info.rtol = tolall(wk.niter) ;
      else
        info.rtol = rtol ;
      end
      info.ierr = ierr ;
      info.xiter = xiter ;
      info.natlevel = sumxall ;
      info.simlevel = sumxqall ;
      info.stdlevel = dlevfall ;
      info.precision = tolall ;
      info.dampingfc = fcall ;
      info.niter = wk.niter ;
      info.ncorr = wk.ncorr ;
      info.nrejr1 = wk.nrejr1 ;
      info.njac = wk.njac ;
      info.nfcn = wk.nfcn ;
      info.nfcnj = wk.nfcnj ;
      
%     End of function NLSCON

function [ierr,iopt,rtol,xscal,fscal] = ncpchk(n,m,mfit,x,xscal,fi,fscal,rtol,iopt,wk)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     N C P C H K : Checking of input parameters and options
%                   for NLSCON.
%
%*    Parameters:
%     ===========
%
%     See parameter descriptions in driver routine and NCINT.
%
%*    functions called: D1MACH
%
%*    Machine dependent constants used:
%     =================================
%
%     EPMACH = relative machine precision
%     GREAT = squareroot of maxreal divided by 10
%     SMALL = squareroot of "smallest positive machine number
%             divided by relative machine precision"
%
%     ------------------------------------------------------------
      numopt=50 ;
%
%      DATA IOPTL = [0,0,0,0,0,0,0,0,0,0,0,1,0,1,0,1,0,0,0,1, ...
%                    0,0,0,0,0,0,0,0,0,0,0,0,0,0,0, ...
%                    0,0,0,0,0,0,0,0,0,0, ...
%                    -9999999,-9999999,-9999999,-9999999,-9999999]
%      DATA IOPTU = [1,1,3,0,0,0,0,0,1,0,3,99,6,99,3,99,0,0,1,99, ...
%                    1,1,0,0,0,0,0,0,0,0,4,1,0,0,1, ...
%                    2,0,2,0,0,0,0,0,0,0, ...
%                    9999999,9999999,9999999,9999999,9999999]
%
      epmach = d1mach(3) ;
      great  = sqrt(d1mach(2)/10.0) ;
      small  = d1mach(6) ;
      ierr = 0 ;
%        print error messages?
      mprerr = iopt.mprerr ;
      luerr  = iopt.luerr  ;
      if luerr  <=  0  ||  luerr  >  99
        luerr = 1 ;
        iopt.luerr=luerr ;
      end
%
%     Checking dimensional parameters N,M and MFIT
      if  n <= 0  ||  m <= 0  ||  mfit < 0  ||  mfit > m 
        if mprerr >= 1
          fprintf(luerr,'\nError: %s\n        %s\n        %s%5i%s%5i%s%5i\n', ...
                  'Bad or inconsistent input to dimensional parameters supplied', ...
                  'choose N and M positive, and MFIT <= M, MFIT nonnegative', ...
                  'your input is: N = ',n,' M = ',m,' MFIT = ',mfit) ;
        end
        hier1=n
        hier2=m
        hier3=mfit
        ierr = 20 ;
      end
%
%     Problem type specification by user
      [nonlin,iopt] = getopt(iopt,'nonlin',0) ;
      if nonlin == 0
        nonlin = 3 ;
      end
      iopt.nonlin = nonlin ;
%
%     Checking and conditional adaption of the user-prescribed RTOL
      if rtol <= 0.0
        if mprerr >= 1
          fprintf(luerr,'\n%s\n',' Error: Nonpositive RTOL supplied') ;
        end
        ierr = 21 ;
      else
        tolmin = epmach*10.0*n ;
        if rtol < tolmin
          rtol = tolmin ;
          if mprerr >= 2
            fprintf(luerr, ...
            '\n warning: user prescribed rtol %s to reasonable %s value rtol = %11.2e\n', ...
            'increased','smallest',rtol ) ;
          end
        end
        tolmax = 1.0e-1 ;
        if rtol > tolmax
          rtol = tolmax ;
          if mprerr >= 2 
            fprintf(luerr, ...
            '\n warning: user prescribed rtol %s to reasonable %s value rtol = %11.2e\n', ...
            'decreased','largest',rtol) ;
          end
        end
      end
%     
%     Test user prescribed accuracy and scaling on proper values
      if n <= 0
        return ;
      end
      if nonlin >= 3
        defscl = rtol ;
      else
        defscl = 1.0 ;
      end
      for i=1:n
        if xscal(i) < 0.0
          if mprerr >= 1
            fprintf(luerr,' error: negative value in xscal(%5i) supplied\n', i) ;
          end
          ierr = 22 ;
        end
        if xscal(i) == 0.0
          xscal(i) = defscl ;
        end
        if  xscal(i) > 0.0  &&  xscal(i) < small 
          if mprerr >= 2
            fprintf(luerr,' warning: xscal(%5i) = %9.2e too small, increased to %9.2e\n', ...
                    i,xscal(i),small ) ;
          end
          xscal(i) = small ;
        end
        if xscal(i) > great
          if mprerr >= 2
            fprintf(luerr,' warning: xscal(%5i) = %9.2e too big, decreased to %9.2e\n', ...
                    i,xscal(i),great ) ;
          end
          xscal(i) = great ;
        end
      end
%     Checks options
%      DO 20 I=1,30
%        if IOPT(I) < IOPTL(I)  ||  IOPT(I) > IOPTU(I)
%          IERR=30
%          if MPRERR >= 1
%            WRITE(LUERR,20001) I,IOPT(I),IOPTL(I),IOPTU(I)
%20001       FORMAT(' Invalid option specified: IOPT(',I2,')=',I12,';',
%     $             /,3X,'range of permitted values is ',I8,' to ',I8)
%          end
%        end
%20    CONTINUE
%     End of function NCPCHK

function [x,xscal,fi,rtol,irank,ierr,xl,xr,vcv,                            ...
          aa,a,qa,dx,dxq,xa,fmodel,f,fa,eta,fw,xw,dxqa,qu,rq,rqkeep,delxq, ...
          fc,fcmin,sigma,fca,cond,conv,sumx,sumxs,dlevf,                   ...
          xiter,sumxall,dlevfall,sumxqall,tolall,fcall,                    ...
          niter,ncorr,nfcn,njac,nfcnj,nrejr1,new,ifccnt,qsucc] =           ...
ncint(n,m,mfit,fcn,par,x,xscal,fi,rtol,nitmax,nonlin,                      ...
      irank,iopt,ierr,wk,m1,m2,maxmn,                                      ...
      aa,a,qa,dx,dxq,xa,fmodel,f,fa,eta,fw,xw,dxqa,qu,rq,rqkeep,delxq,     ...
      fc,fcmin,sigma,fca,cond,conv,sumx,sumxs,dlevf,                       ...
      xiter,sumxall,dlevfall,sumxqall,tolall,fcall,                        ...
      tolmin,mprerr,mprmon,mprsol,luerr,lumon,lusol,niter,ncorr,nfcn,njac, ...
      nfcnj,nrejr1,new,ifccnt,qbdamp)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     N C I N T : Core routine for NLSCON .
%     Damped Gauss-Newton-algorithm with rank-strategy for highly
%     nonlinear least squares estimation problems especially
%     designed for numerically sensitive problems.
%
%*    Parameters:
%     ===========
%
%       N          Int    Number of parameters to be fitted
%       M          Int    Number of equality constraints plus
%                            number of measurements
%       MFIT       Int    Number of measurements
%       FCN,X,XSCAL,FI,RTOL   
%                         See parameter description in driver routine
%
%       NITMAX     Int    Maximum number of allowed iterations
%       NONLIN     Int    Problem type specification
%                         (see IOPT-field NONLIN)
%       IRANK      Int    Initially proposed (in) and final (out) rank
%                         of Jacobian
%       IOPT       Int    See parameter description in driver routine
%       IERR       Int    See parameter description in driver routine
%       WK         Struct Workspace structure array
%       M1         Int    Leading dimension of Jacobian array A
%                         for full case Jacobian: M
%                         (other matrix types are not yet implemented)
%       M2         Int    Leading dimension of Jacobian array AA
%                         for full case Jacobian: M
%       MAXMN      Int    Max(M,N) - length of some temporary
%                         workspace
%       VCV(N,N)   Dble   Correlation matrix, if IOPT.qstat=1 or
%                         IOPT.mprsta=1 (out)
%       XL(N)     Dble    Confidence interval left sides,
%                         if IOPT.qstat=1 or IOPT.mprsta=1 (out)
%       XR(N)     Dble    Confidence interval right sides,
%                         if IOPT.qstat=1 or IOPT.mprsta=1 (out)
%       AA(M2,N)   Dble   Holds the originally computed Jacobian
%       A(M1,N)    Dble   Holds the Jacobian matrix (decomposed form
%                         after call of linear decomposition
%                         routine)
%       QA(N,N)    Dble   Holds the pseudo inverse in case of rank-
%                         deficiency
%       DX(N)      Dble   Current Gauss-Newton correction
%       DXQ(N)     Dble   Simplified Gauss-Newton correction J(k)*X(k+1)
%       XA(N)      Dble   Previous Gauss-Newton iterate
%       F(M)       Dble   Function (FCN) value minus measured values
%                         vector FI of current iterate
%       FA(M)      Dble   Holds the previous values of vector F(M)
%       ETA(N)     Dble   Jacobian approximation: updated scaled
%                         denominators
%       FW(M)      Dble   Scaling factors for rows of the system
%       XW(N)      Dble   Scaling factors for iteration vector
%       DXQA(N)    Dble   Previous simplified Gauss-Newton 
%                         correction J(k-1)*X(k)
%       QU(M)      Dble   Savespace for right hand side belonging
%                         to upper triangular linear system
%       RQ(M)      Dble   Gets the residuum of the linear problems
%                         solution in each iteration step :
%                         JAC(X(i-1))*DXQ(i)+F(i) for the i-th iterate
%       RQKEEP(M)  Dble   Keeps a copy of RQ(M) for restoring it in
%                         case of Jacobian recomputation (M <= N ,
%                         rejected rank-1 step, and Jacobian of
%                         previous iterate was rank-deficient).
%       DELXQ(N)   Dble   Gets the vector of parameters best fitting
%                         the linear problems residuum of the previous
%                         iterate
%       FC         Dble   Current Gauss-Newton iteration damping factor.
%       FCMIN      Dble   Minimum permitted damping factor. If
%                         FC becomes smaller than this value, one
%                         of the following may occur:
%                         a.    Recomputation of the Jacobian
%                               matrix by means of difference
%                               approximation (instead of Rank1
%                               update), if Rank1 - update
%                               previously was used
%                         b.    Rank reduction of Jacobi
%                               matrix ,  if difference
%                               approximation was used previously
%                               and Rank(A) ~= 0
%                         c.    Fail exit otherwise
%       SIGMA      Dble   Decision parameter for rank1-updates.
%       FCA        Dble   Previous Gauss-Newton iteration damping factor.
%       COND       Dble   Maximum permitted subcondition for rank-
%                         decision by linear solver.
%       CONV       Dble   Scaled maximum norm of the Gauss-Newton-
%                         correction. Passed to RWK-field on output.
%       SUMX       Dble   Square of the natural level (see equal-
%                         named IOPT-output field)
%       SUMXS       Dble   Square of the "simplified" natural level
%                          (see equal-named RWK-internal field)
%       DLEVF      Dble   Square of the standard level (see equal-
%                         named IOPT-output field)
%       MPRERR,MPRMON,MPRSOL,LUERR,LUMON,LUSOL :
%                         See description of equal named IOPT-fields
%                         in the driver function
%       NITER,NCORR,NFCN,NJAC,NFCNJ,NREJR1,NEW :
%                         See description of equal named WK-fields
%                         in the driver function
%       QBDAMP     Logic  Flag, that indicates, whether bounded damping
%                         strategy is active:
%                         1  = bounded damping strategy is active
%                         0 = normal damping strategy is active
%
%
%*    Internal variables
%     ==================
%
%       AJDEL    See RWK(26) (num. diff. without feedback)
%       AJMIN    See RWK(27) (num. diff. without feedback)
%       COND1    Gets the subcondition of the linear system
%                as estimated by the linear solver (NCFACT)
%       CONVA    Holds the previous value of CONV .
%       DELQ     Gets the projection defect in case of rank-
%                deficiency.
%       DMUE     Temporary value used during computation of damping 
%                factors predictor.
%       EPDIFF   sqrt(10*epmach) (num. diff. with feedback)
%       ETADIF   See description of RWK(28) (num. diff. with feedback)
%       ETAINI   Initial value for all ETA-components (num. diff. fb.)
%       ETAMAX   Maximum allowed pertubation (num. diff. with feedback)
%       ETAMIN   Minimum allowed pertubation (num. diff. with feedback)
%       FCDNM    Used to compute the denominator of the damping 
%                factor FC during computation of it's predictor,
%                corrector and aposteriori estimate (in the case of
%                performing a Rank1 update) .
%       FCK2     Aposteriori estimate of FC.
%       FCMIN2   FCMIN^2 . Used for FC-predictor computation.
%       FCMINH   sqrt(FCMIN).
%                Used in rank decision logical expression.
%       FCNUMP   Gets the numerator of the predictor formula for FC.
%       FCNMP2   Temporary used for predictor numerator computation.
%       FCNUMK   Gets the numerator of the corrector computation 
%                of FC .
%       SENS1    Gets the sensitivity of the Jacobian as
%                estimated by the linear solver NCFACT.
%       SKAP     In case of termination at stationary point:
%                incompatibility factor
%       SUMXA    Natural level of the previous iterate.
%       TH       Temporary variable used during corrector- and 
%                aposteriori computations of FC.
%
%
%     IFAIL      Gets the return value from functions called from
%                NCINT (NCFACT, NCFIT, FCN, JAC)
%     ISCAL      Holds the scaling option from the IOPT-field ISCAL      
%     MODE       Matrix storage mode (see IOPT-field MODE) 
%     NRED       Count of successive corrector steps
%
%
%     QGENJ      Jacobian updating technique flag:
%                =1  : Call of analytical function Jacobian or
%                           numerical differentiation
%                =0 : rank1- (Broyden-) update
%     QINISC     Iterate initial-scaling flag:
%                =1  : at first call of NCSCAL
%                =0 : at successive calls of NCSCAL
%     QSUCC      See description of IOPT-field QSUCC.
%     QJCRFR     Jacobian refresh flag:
%                set to 1 if damping factor gets too small
%                and Jacobian was computed by rank1-update. 
%                Indicates, that the Jacobian needs to be recomputed
%                by function JAC or numerical differentation.
%     QLINIT     Initialization state of linear solver workspace:
%                =0 : Not yet initialized
%                =1  : Initialized - NCFACT has been called at
%                           least one time.
%     QREPET     Operation mode flag for linear solver:
%                =0 : Normal operation (full rank matrix)
%                =1  : Special operation in rank deficient case:
%                           Compute rank-deficient pseudo-inverse,
%                           partial recomputation when solving the
%                           linear system again prescribing a lower
%                           rank as before.
%     QNEXT      Set to 1 to indicate success of the current
%                Gauss-Newton-step, i.e. : sucessfull monotonicity-test.
%     
%     QREDU      Set to 1 to indicate that rank-reduction (or
%                refreshment of the Jacobian) is needed - if the
%                computed damping factor gets too small.
%     QSCALE     Holds the value of  ~ QNSCAL. See description
%                of IOPT-field QNSCAL.
%
%*    Functions called:
%     ===================
%
%       NCFACT, NCFIT,  NCJAC,  NCJCF, NCLVLS, NCSCRF, NCPRJN,
%       NCRNK1, NCSOUT, NCPRV1, NCPRV2, NCSCAL,
%       MONON,  MONOFF, D1MACH, WNORM
%
%*    Machine constants used
%     ======================
%
%       EPMACH,SMALL
% 
%     ------------------------------------------------------------
      one=1.0 ; zero=0.0 ;  half=0.5 ;  ten=10.0 ;
      epmach = d1mach(3) ;
      small  = d1mach(6) ;
%*    Begin
%       ----------------------------------------------------------
%       1 Initialization
%       ----------------------------------------------------------
        xl = [] ;
        xr = [] ;
        vcv = [] ;
%       1.1 Control-flags and -integers
        qsucc = wk.qsucc ; 
        qscale = iopt.norowscal ~= 1 ; 
        qrank1 = iopt.qrank1 ; 
        qmixio = lumon == lusol  &&  mprmon ~= 0  &&  mprsol ~= 0 ; 
        [iscal,iopt] = getopt(iopt,'iscal',0) ;
        [mode,iopt]  = getopt(iopt,'mode' ,0) ;
        [iterm,iopt] = getopt(iopt,'iterm',0) ; 
        jacgen = iopt.jacgen ; 
        mprtim = iopt.mprtim ; 
%       ----------------------------------------------------------
%       1.2 Derivated dimensional parameters
        mcon = m-mfit ;
%       ----------------------------------------------------------
%       1.3 Derivated internal parameters
        minmn = min(m,n) ;
        fcmin2 = fcmin*fcmin ;
        fcminh = sqrt(fcmin) ;
        rsmall = sqrt(ten*rtol) ;
%       ----------------------------------------------------------
%       1.4 Adaption of input parameters, if necessary
        if fc < fcmin
          fc = fcmin ;
        end
        if fc > one
          fc = one ;
        end
%       ----------------------------------------------------------
%       1.5 Initial preparations
        qjcrfr = 0 ;
        qlinit = 0 ;
        qiter = 1 ;
        qrepet = 0 ;
        ifail = 0 ;
        d = zeros(n,1) ;
        pivot = zeros(n,1) ;
        fcbnd = zero ;
        if qbdamp
          fcbnd = wk.fcbnd ;
        end
%       ----------------------------------------------------------
%       1.5.1 Numerical differentation related initializations
        if jacgen == 2
          [ajdel,wk] = getopt(wk,'ajdel',0.0) ;
          if ajdel <= small
            ajdel = sqrt(epmach*ten) ;
          end
          [ajmin,wk] = getopt(wk,'ajmin',0.0) ;
        elseif jacgen == 3
          [etadif,wk] =  getopt(wk,'etadif',0.0) ;
          if etadif  <=  small
            etadif = 1.0e-6 ;
          end
          [etaini,wk] =  getopt(wk,'etaini',0.0) ;
          if etaini  <=  small
            etaini = 1.0e-6 ;
          end
          epdiff = sqrt(epmach*ten) ;
          etamax = sqrt(epdiff) ;
          etamin = epdiff*etamax ;
        end
%       ----------------------------------------------------------
%       1.5.2 Miscellaneous preparations of first iteration step
        if  ~ qsucc
          niter = 0 ;
          wk.niter = niter ;
          xiter(1:n,1) = x ;
          rq = zeros(m,1) ;
          rqkeep = zeros(m,1) ;
          ncorr = 0 ;
          nrejr1 = 0 ;
          nfcn = 0 ;
          njac = 0 ;
          nfcnj = 0 ;
          ifccnt = 0 ;
          iranka = irank ;
          qgenj = 1 ;
          qinisc = 1 ;
          fca = fc ;
          fck2 = fc ;
          conv = zero ;
          sumxk = zero ;
          if jacgen == 3
            eta = etaini*ones(n,1) ;
          end
          xa = x ;
%         ------------------------------------------------------
%         1.6 Print monitor header
          if mprmon >= 2 && ~ qmixio
            fprintf(lumon,'\n\n\n%1s',' ') ;
            fprintf(lumon,'%s\n\n%s\n', ...
                  '  **************************************************************************' , ...
                  '        It       Normf              Normx      Damp.Fct.   New      Rank'        ) ;
          end
%         --------------------------------------------------------
%         1.7 Startup step
%         --------------------------------------------------------
%         1.7.1 Computation of the residual vector
          if mprtim ~= 0
             monon(1) ;
          end
          [fmodel,ifail] = feval(fcn,x,'',par) ;
          if mprtim ~= 0
            monoff(1) ;
          end
          nfcn = nfcn+1 ;
%     Exit, if ...
          if ifail ~= 0
            ierr = 82 ;
            qiter = 0 ;
          else
            f = zeros(m,1) ;
            f(mcon+1:mcon+mfit) = fmodel(mcon+1:mcon+mfit)-fi(1:mfit) ;
            f(1:mcon)=fmodel(1:mcon) ;
          end
        else
          qinisc = 0 ;
        end
        dummy = 1 ;
%
%       Main iteration loop
%       ===================
%
%       Repeat
        while qiter
%         --------------------------------------------------------
%         2 Startup of iteration step
          if  ~ qjcrfr
%           ------------------------------------------------------
%           2.0 Scaling of variables X(N)
            xw = ncscal(n,x,xa,xscal,iscal,qinisc) ;
            qinisc = 0 ;
            if niter ~= 0
%             Preliminary pseudo-rank
              iranka = irank ;
              dxqa = dxq ;
%             ----------------------------------------------------
%             2.1 Computation of linear residuum RQ(M)
              if m > n || iranka ~= m
                rq(1:m) = zeros(m,1) ;
                for l1=mcon+1:m
                  rq(l1)=sum(aa(l1,1:n).*dxqa(1:n)') + f(l1) ;
                end
                rqkeep = rq ;
                mfout=mfit ;
              else
                mfout=0 ;
              end
%             ----------------------------------------------------
%             2.2 Aposteriori estimate of damping factor
              fcnump = sum((dx./xw).^2) ;
              th = fc-one ;
              fcdnm = sum(((dxqa+th*dx)./xw).^2) ;
              fch = fc*fc*half*sqrt(fcnump/fcdnm) ;
%             ------------------------------------------------------
%             2.2.1 Computation of the numerator of damping
%                   factor predictor
              fcnmp2 = sum((dxqa./xw).^2) ;
              fcnump = fcnump*fcnmp2 ;
%             ----------------------------------------------------
%             2.2.2 Decision criterion for Jacobian updating
%                   technique
              qgenj = fc < fca && new > 0 || fch < fc*sigma || iranka ~= n || m > n || ~ qrank1 ;
              fca = fc ;
              fck2 = fca ;
              irank = minmn ;
              if nonlin > 1
                fc = min(one,fch) ;
              end
%             ----------------------------------------------------
%             2.2.3 Denominator for kappa (SKAP) estimate
              sumxk = sum((dx./xw).^2) ;
            end
          end
          qjcrfr =0 ;
%         --------------------------------------------------------
%         2.3 Jacobian matrix (stored to array AA(M2,N))
%         --------------------------------------------------------
%         2.3.1 Jacobian generation by routine JAC or
%               difference approximation (If QGENJ == 1)
%               - or -
%               Rank-1 update of Jacobian (If QGENJ == 0)
          if qgenj
            new = 0 ;
            if jacgen == 1
               if mprtim ~= 0
                 monon(2) ;
               end
               [aa,ifail] = feval(fcn,x,'jacobian',par) ;
               if mprtim ~= 0
                 monoff(2) ;
               end
            else
                if mprtim ~= 0
                  monon(2) ;
                end
                if jacgen == 3
                  [aa,eta,nfcnj,ifail] = ncjcf(fcn,par,n,m,x,fmodel,xw,eta, ...
                                         etamin,etamax,etadif,conv,nfcnj) ;
                end
                if jacgen == 2 
                  [aa,nfcnj,ifail] = ncjac(fcn, par, n, m, x, fmodel, xw, ...
                                           ajdel, ajmin, nfcnj) ;
                end
                if mprtim ~= 0
                  monoff(2) ;
                end
            end
            njac = njac + 1 ;
%     Exit, If ...
            if jacgen == 1  &&  ifail < 0
              ierr = 83 ;
              break ;
            end
            if jacgen ~= 1  &&  ifail ~= 0
              ierr = 82 ;
              break ;
            end
          else
            new = new+1 ;
            aa = ncrnk1(n,m,xw,dx,f,fa,aa,fca) ;
          end
%         --------------------------------------------------------
%         2.3.2 Copy matrix to Jacobian work array A(M2,N)
          a(1:m2,1:n)=aa(1:m2,1:n) ;
%         --------------------------------------------------------
%         2.4 Prepare solution of the linear system
%         --------------------------------------------------------
%         2.4.1 internal column scaling of matrix A
          for k=1:n
            a(1:m,k) = -a(1:m,k)*xw(k) ;
          end
%         ------------------------------------------------------
%         2.4.2 Row scaling of matrix A(M,N)
          if qscale
            [a,fw] = ncscrf(n,m,mcon,a,fw) ;
          else
            fw = ones(m,1) ;
          end
%         --------------------------------------------------------
%         2.4.3 Save and scale values of F(M)
          fa = f ;
          t3 = f .* fw ;
%         --------------------------------------------------------
%         2.4.4 Scaling of the linear residuum RQ(M)
          if m > n || iranka ~= m
             rq(1:m) = rq(1:m) .* fw(1:m) ;
          end
          qnext = 0 ;
          qredrnk = 1 ;
%         --------------------------------------------------------
%         3 Central part of iteration step
%
%         Pseudo-rank reduction loop
%         ==========================
          while qredrnk
%         DO (Until)
%3        CONTINUE
%           --------------------------------------------------------
%           3.1 Solution of linear (MFIT,N)-least squares problem
%               with MCON equality constraints
%           --------------------------------------------------------
%           3.1.1 Decomposition of (N,N)-matrix A
            wk.niwla(1) = qrepet ;
            cond1 = cond ;
            mconh = mcon ;
            if mprtim ~= 0
               monon(3) ;
            end
            [a,d,cond1,qa,pivot,irank,wk,ifail] = ncfact(n,m1,mconh,n,a,d,pivot,cond1,irank,iopt,wk) ;
            if mprtim ~= 0
               monoff(3) ;
            end
            qlinit = 1 ;
%     Exit Repeat If ...
            if ifail ~= 0
              ierr = 80 ;
              qiter = 0 ;
              break ;
            end
            irankc = wk.niwla(2) ;
            sens1  = wk.nrwla(1) ;
            condco = wk.nrwla(2) ;
%           --------------------------------------------------------
%           3.1.2 Solution of linear (N,N)-system
            if mprtim ~= 0
              monon(4) ;
            end
            [t3,t2,ifail] = ncfit(n,m1,mconh,n,a,qa,t3,irank,d,pivot,wk) ;
            if mprtim ~= 0
              monoff(4) ;
            end
%     Exit Repeat If ...
            if ifail ~= 0
              ierr = 81 ;
              qiter = 0 ;
              break ;
            end
            if  ~ qrepet && irank ~= 0
               qu(1:m)=t3(1:m) ;
            end
%           --------------------------------------------------------
%           3.2 Evaluation of scaled natural level function SUMX
%               scaled maximum error norm CONV
%               evaluation of (scaled) standard level function
%               DLEVF ( DLEVF only, if MPRMON >= 2 )
%               and computation of ordinary Gauss-Newton corrections DX(
%               N)
            [dx,conv,sumx,dlevf] = nclvls(n,m,t2,xw,f,mprmon) ;
            xa =x ;
            sumxa = sumx ;
            dlevxa = sqrt(sumxa/n) ;
            conva = conv ;
            dxanrm = wnorm(n,dx,xw) ;
            sumxall(niter+1) = dlevxa ;  dlevfall(niter+1) = dlevf ;
%           --------------------------------------------------------
%           3.3 A - priori estimate of damping factor FC
            qredu = 0 ;
            if niter ~= 0 && nonlin ~= 1
              if new == 0 || qrepet
%               ------------------------------------------------------
%               3.3.1 Comp. of the denominator of a-priori estimate
                if m > n || iranka ~= m
                  if mprtim ~= 0
                    monon(4) ;
                  end
                  [rq,delxq,ifail] = ncfit(n,m1,mconh,n,a,qa,rq,irank,d,pivot,wk) ;
                  %delxq=0;
                  %norm(delxq)
                  if mprtim ~= 0
                    monoff(4) ;
                  end
%     Exit Repeat If ...
                  if ifail ~= 0
                    ierr = 81 ;
                    qiter = 0 ;
                    break ;
                  end
                  fcdnm = sum(((dx-dxq)./xw-delxq).^2) ;
                else
                  fcdnm = sum(((dx-dxq)./xw).^2) ;
                end
                if irank ~= n
%                 ------------------------------------------------
%                 3.3.2 Rank-deficient case (if previous rank
%                           was full) computation of the projected
%                       denominator of a-priori estimate
                  t1 =dxqa ./xw ;
%                 Norm of projection of reduced component T1(N)
                  delq = ncprjn(n,irank,t1,d,qa,pivot) ;
                  fcdnm = fcdnm-delq ;
                end
                fcdnm = fcdnm*sumx ;
%               ------------------------------------------------------
%               3.3.3 New damping factor
                if fcdnm > fcnump*fcmin2
                  dmue = fca*sqrt(fcnump/fcdnm) ;
                  fc = min(dmue,one) ;
                  if dmue > ten
                    ifccnt = ifccnt + 1 ;
                  end
                else
                  fc = one ;
%$test-begin
                  dmue = -1.0 ;
%$test-end
                  if fca >= one
                    ifccnt = ifccnt+1 ;
                  end
                end
%$Test-begin
                if mprmon >= 5
                   fprintf(lumon,'%s\n%s%10i          %s%18.10e\n%s%18.10e  %s%18.10e\n%s%18.10e  %s%18.10e\n%s\n', ...
                      ' +++ apriori estimate +++',  ...
                      '  ifccnt = ', ifccnt, ...
                      '  fc     = ', fc,     ...
                      '  fca    = ', fca,    ...
                      '  dmue   = ', dmue,   ...
                      '  fcnump = ', fcnump, ...
                      '  fcdnm  = ', fcdnm,  ...
                      ' ++++++++++++++++++++++++');
                end
%$Test-end 
                if qbdamp
                  fcbh = fca*fcbnd ;
                  if fc > fcbh
                    fc = fcbh ;
                    if mprmon >= 4 
                       fprintf(lumon,'%s\n',' *** incr. rest. act. (a prio) ***') ;
                    end
                  end
                  fcbh = fca/fcbnd ;
                  if fc < fcbh
                    fc = fcbh ;
                    if mprmon >= 4
                       fprintf(lumon,'%s\n',' *** decr. rest. act. (a prio) ***') ; 
                    end
                  end
                end
              end
              qredu = fc < fcmin ;
            end
            qrepet = 0 ;
            if  ~ qredu
%             --------------------------------------------------------
%             3.4 Save natural level for later computations of
%                 corrector and print iterate
              fcnumk = sumx ;
              if mprmon >= 2
                if mprtim ~= 0
                  monon(5) ;
                end
                ncprv1(dlevf,dlevxa,fca,niter,new,irank,mprmon,lumon,qmixio) ;
                if mprtim ~= 0
                  monoff(5) ;
                end
              end
              nred = 0 ;
              qred = 1 ;
%             Damping-factor reduction loop
%             ================================
              while qred
%             DO (Until)
%34           CONTINUE
%               ------------------------------------------------------
%               3.5 Preliminary new iterate
                x = xa + dx*fc ;
                fcall(niter+1) = fc ;
%               -----------------------------------------------------
%               3.5.2 Exit, if problem is specified as being linear
%     Exit Repeat If ...
                if  nonlin == 1 
                  ierr = 0 ;
                  qiter = 0 ;
                  break ;
                end
%               ------------------------------------------------------
%               3.6.1 Computation of the residual vector
                if mprtim ~= 0
                  monon(1) ;
                end
                [fmodel,ifail] = feval(fcn,x,'',par) ;
                if mprtim ~= 0
                  monoff(1) ;
                end
                nfcn = nfcn+1 ;
%     Exit, if ...
                if ifail < 0
                  ierr = 82 ;
                  qiter = 0 ;
                  break ;
                end
                if ifail == 1  ||  ifail == 2
                  if ifail == 1
                    fcredu = half ;
                  else
                    fcredu = fmodel(1) ;
%     exit, if ...
                    if fcredu <= 0  ||  fcredu >= 1
                      ierr = 83 ;
                      qiter = 0 ;
%                       break ;
                    end
                  end
                  if mprmon >= 2
                    fprintf(lumon,'%8s%2i%41s%5.3f%4s%2i      %4i\n',          ...
                                 '        ',niter,                            ...
                                 ' fcn could not be evaluated              ', ...
                                 fc,'    ',new ,irank) ;
                  end
                  fch = fc ;
                  fc = fcredu*fc ;
                  if fch > fcmin
                    fc = max(fc,fcmin) ;
                  end
                  if qbdamp
                    fcbh = fch/fcbnd ;
                    if fc < fcbh
                      fc = fcbh ;
                      if mprmon >= 4
                        fprintf(lumon,' *** decr. rest. act. (fcn redu.) ***\n') ;
                      end
                    end
                  end
                  if fc < fcmin
                    ierr = 3 ;
                    qiter = 0 ;
                    break ;
                  end  
%     Break DO (Until) ...
                  break ;
                end
                f(mcon+1:mcon+mfit) = fmodel(mcon+1:mcon+mfit)-fi(1:mfit) ;
                f(1:mcon) = fmodel(1:mcon) ;
                t3(1:m) = f(1:m) .* fw(1:m) ;
%               ------------------------------------------------------
%               3.6.2 Solution of linear (MFIT,N)-least squares problem
%                     with MCON equality constraints
                wk.niwla(1) = qrepet ;
                if mprtim ~= 0
                  monon(4) ;
                end
                [t3,t2,ifail] = ncfit(n,m1,mconh,n,a,qa,t3,irank,d,pivot,wk) ;
                if mprtim ~= 0
                  monoff(4) ;
                end
%     Exit Repeat If ...
                if ifail ~= 0
                  ierr = 81 ;
                  qiter = 0 ;
                  break ;
                end
%               ------------------------------------------------------
%               3.6.3 Evaluation of scaled natural level function
%                     SUMX
%                     scaled maximum error norm CONV and evaluation
%                     of (scaled) standard level function DLEVF
                [dxq,conv,sumx,dlevf] = nclvls(n,m,t2,xw,f,mprmon) ;
                sumxqall(niter+1) = sqrt(sumx/n) ; 
                dxnrm = wnorm(n,dxq,xw) ;
%               ------------------------------------------------------
%               3.6.4 Convergence test
%     Exit Repeat If ...
                tolall(niter+1) = dxnrm ;
                if iterm == 0
                  if  (dxnrm <= rtol  &&  ifccnt >= 3)  ||  dxanrm <= rtol
                    ierr = 0 ;
                    qiter = 0 ;
                    break ;
                  end
                elseif iterm == 1
                  if  dxnrm <= rtol  &&  ifccnt >= 3 
                    ierr = 0 ;
                    qiter = 0 ;
                    break ;
                  end
                elseif iterm == 2
                  if  dxanrm <= rtol 
                    ierr = 0 ;
                    qiter = 0 ;
                    break ;
                  end
                end
%           
                fca = fc ;
%               ------------------------------------------------------
%               3.7 Natural monotonicity test
                if sumx > sumxa
%                 ----------------------------------------------------
%                 3.8 Output of iterate
                  if mprmon >= 3
                    if mprtim ~= 0
                      monon(5) ;
                    end
                    ncprv2(dlevf,sqrt(sumx/n),fc,niter,mprmon,lumon,qmixio,'*') ;
                    if mprtim ~= 0
                      monoff(5) ;
                    end
                  end
%                 ----------------------------------------------------
%                 3.9 Evaluation of reduced damping factor
                  th = fca-one ;
                  fcdnm = sum( ((dxq+th*dx) ./ xw) .^ 2 ) ;
                  fc = fca*fca*half*sqrt(fcnumk/fcdnm) ;
                  if qbdamp
                    fcbh = fca/fcbnd ;
                    if fc < fcbh
                      fc = fcbh ;
                      if mprmon >= 4
                        fprintf(lumon,' *** decr. rest. act. (a post) ***\n') ;
                      end
                    end
                  end
                  ncorr = ncorr+1 ;
                  nred = nred+1 ;
                  ifccnt = 0 ;
%                 ----------------------------------------------------
%                 3.10 Rank reduction, if damping factor to small
                  qredu  = fc < fcmin || (new > 0 && nred > 1) ;
                else
                  qnext = 1 ;
                end
                qred =  ~ (qnext || qredu) ;
%             UNTIL ( expression - negated above)
              end
              if ~ qredrnk || ~ qiter
                break ;
              end ;
%             End of damping-factor reduction loop
%           =======================================
            end
            if qredu
%             ------------------------------------------------------
%             3.11 Restore former values for repeting step
%                  step
              nrejr1 = nrejr1+1 ;
              x = xa ;
              f(1:m) = fa(1:m) ;
              dxq = dxqa ;
              if mprmon >= 2
                fprintf(lumon,'        %2i %40s%5.3f    %2i      %4i\n',                ...
                        niter,'not accepted damping factor             ',fc,new,irank ) ;
              end
              ifccnt = 0 ;
              fca = fck2 ;
              if niter == 0
                fc = fcmin ;
              end
              if new > 0
                qgenj = 1 ;
                qjcrfr = 1 ;
                qredu = 0 ;
                fc = fch ;
                irank = minmn ;
                rq(1:m) = rqkeep(1:m) ;
              else
%               ------------------------------------------------
%               3.12 Pseudo-rank reduction
                qrepet = 1 ;
                t3(1:m) = qu(1:m) ;
                irank = irank-1 ;
                if irank == 0
                  ierr = 3 ;
                  qiter = 0 ;
                  break ;
                end
              end
            end
            qredrnk = qredu ;
%         UNTIL ( expression - negated above)
%
          end
          if ~ qiter
            break ;
          end
%         End of pseudo-rank reduction loop
%         =================================
          if qnext
%           ------------------------------------------------------
%           4 Preparations to start the following iteration step
%           ------------------------------------------------------
%           4.1 Print values
            if mprmon >= 3
              if mprtim ~= 0
                monon(5) ;
              end
              ncprv2(dlevf,sqrt(sumx/n),fc,niter+1,mprmon,lumon,qmixio,'*') ;
              if mprtim ~= 0
                monoff(5) ;
              end
            end
%           print the natural level of the current iterate and return
%           it in one-step mode
            sumxs = sumx ;
            sumx = sumxa ;
            if mprsol >= 2 && niter ~= 0
              if mprtim ~= 0
                monon(5) ;
              end
              ncsout(n,mfout,xa,rqkeep(mcon+1:m),2,iopt,wk,mprsol,lusol) ;
              if mprtim ~= 0
                monoff(5) ;
              end
            elseif mprsol >= 1 && niter == 0
              if mprtim ~= 0
                monon(5) ;
              end
              ncsout(n,0,xa,dummy,1,iopt,wk,mprsol,lusol) ;
              if mprtim ~= 0
                monoff(5) ;
              end
            end
            niter = niter+1 ;
            wk.niter = niter ;
            xiter(1:n,niter+1) = x ;
%     exit repeat if ...
            if niter >= nitmax
              ierr = 2 ;
              qiter = 0 ;
              break ;
            end
            fca = fc ;
%           ------------------------------------------------------
%           4.2 Return, if in one-step mode
% Exit function If ...
            if mode == 1
              qsucc = 1 ;
              return ;
            end
          end
%       End Repeat
        end
%
%       End of main iteration loop
%       ==========================
%       ----------------------------------------------------------
%       9 Exits
%       ----------------------------------------------------------
%       9.1 Solution exit
        if ierr == 0
          if nonlin ~= 1
            x = x + dxq ;
            xiter(1:n,niter+2) = x ;
            if irank < minmn
              ierr = 1 ;
            end
%           print final monitor output
            if mprmon >= 2
              if mprtim ~= 0
                monon(5) ;
              end
              ncprv2(dlevf,sqrt(sumx/n),fc,niter+1,mprmon,lumon,qmixio,'*') ;
              if mprtim ~= 0
                monoff(5) ;
              end
            end
            if mprmon >= 1 && ierr == 0
              fprintf(lumon,'\n\n\n%s %3i %s\n', ...
                      ' Solution of nonlinear least squares problem obtained within ', ...
                      niter+1,' iteration steps') ;
            end
          else
            if mprmon >= 1
              fprintf(lumon,'\n\n\n%s\n\n%s\n', ...
                     ' Least squares solution of linear system of equations obtained by NLSCON', ...
                     ' No estimate available for the achieved relative accuracy') ;
            end
          end
        end
%       ----------------------------------------------------------
%       9.2 Fail exit messages
%       ----------------------------------------------------------
%       9.2.1 Termination at stationary point
        if ierr == 1 && mprerr >= 1
          fprintf(luerr,'\n Iteration terminates at stationary point\n\n') ;
        end
        wk.prec = -one ;
        wk.skap = -one ;
        prec = -one ;
        if (ierr == 0 || ierr == 1) && nonlin ~= 1 && mprmon >= 1
          if sumxk  <  tolmin
            skap = -one ;
          elseif sumxa  <  tolmin
            skap = sqrt(tolmin) ;
          else
            skap = sqrt(sumxa/sumxk) ;
          end
          if skap >= 0
            fprintf(lumon,'\n Incompatibility factor kappa %10.3e\n',skap) ;
          else 
            fprintf(lumon,'\n Incompatibility factor kappa not available\n') ;
          end
          if  ierr == 0  &&  skap < one 
            if skap >= 0
              prec=max((skap/(one-skap))*sqrt(sumxa/n),epmach) ;
            else 
              prec=epmach ;
            end
            fprintf(lumon,'\n Achieved relative accuracy %10.3e\n\n',prec) ;
          end
          wk.prec = prec ;
          wk.skap = skap ;
        end
        rtol = wk.prec ;
%       ----------------------------------------------------------
%       9.2.2 Termination after more than NITMAX iterations
        if ierr == 2 && mprerr >= 1
          fprintf(luerr,'\n\n Iteration terminates after NITMAX = %3i  Iteration steps\n',nitmax) ;
        end
%       ----------------------------------------------------------
%       9.2.3 Gauss-Newton method fails to converge
        if ierr == 3 && mprerr >= 1
          fprintf(luerr,'\n\n Gauss-Newton method fails to converge\n') ;
        end
%       ----------------------------------------------------------
%       9.2.5 Error exit due to linear solver routine NCFACT
        if ierr == 80 && mprerr >= 1
          fprintf(luerr,'\n\n Error %5i signalled by linear solver NCFACT\n',ifail) ;
        end
%       ----------------------------------------------------------
%       9.2.6 Error exit due to linear solver routine NCFIT
        if ierr == 81 && mprerr >= 1
          fprintf(luerr,'\n\n Error %5i signalled by linear solver NCFIT\n',ifail) ;
        end
%       ----------------------------------------------------------
%       9.2.7 Error exit due to fail of user function FCN
        if ierr == 82 && mprerr >= 1
          fprintf(luerr,'\n\n Error %5i signalled by user function FCN\n',ifail) ;
        end
%       ----------------------------------------------------------
%       9.2.7 Error exit due to fail of user function JAC
        if ierr == 83 && mprerr >= 1
          fprintf(luerr,'\n\n Error %5i signalled by user function JAC\n',ifail) ;
        end
        if ierr >= 80 && ierr <= 83
          wk.ifail = ifail ;
        end
        if (ierr == 82 || ierr == 83) && niter <= 1 && mprerr >= 1
          fprintf(luerr,'\n Try to find a better initial guess for the solution\n') ;
        end
%       ----------------------------------------------------------
%       9.3 Common exit
        if mprmon >= 1
          if mcon > 0
            fprintf(lumon,'\n\n   Subcondition ( 1,%4i) of constrained part %10.3e\n',irankc,condco) ;
            fprintf(lumon,'\n\n   Subcondition ( %4i,%4i) of least squares part %10.3e\n', ...
                          irankc+1,irank,cond1) ;
          else
            fprintf(lumon,'\n\n   Subcondition ( 1,%4i) of least squares part %10.3e\n',irank,cond1) ;
          end
          fprintf(lumon,'\n\n   Sensitivity ( lsq ) %10.3e  \n\n',sens1) ;
        end
        %%%%%%%%%%%%Susanna: nice print of ordered subcondition numbers
          if mprmon > 0
              fprintf(lumon,'\n\n Parameters ordered by increasing subcondition number:\n\n') ;
              fprintf(lumon,'Parameter \t subcondition\n\n') ;
              for j=1:irank
                fprintf(lumon,' %3i \t %10.2e \n',pivot(j),abs(d(1)/d(j))) ;  
              end
              fprintf(lumon,'\n') ;
          end
        %%%%%%%%%%%%%%%%
        sumxs = sumx ;
        sumx = sumxa ;
        if mprsol >= 2 && niter ~= 0
          if mprtim ~= 0
            monon(5) ;
          end
          ncsout(n,mfout,xa,rqkeep(mcon+1:m),2,iopt,wk,mprsol,lusol) ;
          if mprtim ~= 0
            monoff(5) ;
          end
        elseif mprsol >= 1 && niter == 0
          if mprtim ~= 0
            monon(5) ;
          end
          ncsout(n,0,xa,dummy,1,iopt,wk,mprsol,lusol) ;
          if mprtim ~= 0
            monoff(5) ;
          end
        end
        niter = niter+1 ;
        wk.niter = niter ;
        if mprsol >= 1
%         Print Solution or final iteration vector
          if ierr == 0
             modefi = 3 ;
          else
             modefi = 4 ;
          end
          if mprtim ~= 0
            monon(5) ;
          end
          ncsout(n,0,x,dummy,modefi,iopt,wk,mprsol,lusol) ;
          if mprtim ~= 0
            monoff(5) ;
          end
        end
%       Return the latest internal scaling to XSCAL
        xscal = xw ;
%       ----------------------------------------------------------
%       9.4 Optional computation of the statistical analysis
%           for the least squares problem solution
        iopt = iniopt(iopt,'qstat',0) ;
        iopt = iniopt(iopt,'mprsta',0) ;
        if iopt.qstat == 1  ||  iopt.mprsta == 1
          iopt.qstat = 1 ;
          if mprtim ~= 0
             monon(1) ;
          end
          [fmodel,ifail] = feval(fcn,x,'',par) ;
          if mprtim ~= 0
            monoff(1) ;
          end
          nfcn = nfcn+1 ;
%     Exit, if ...
          if ifail ~= 0
            ierr = 182 ;
            fprintf(luerr,'\n\n %s\n %s\n', ...
                          'Computation of the statistical analysis skipped', ...
                          'since FCN failed for the final Gauss-Newton iterate') ;
            return ;
          end
          f(mcon+1:mcon+mfit)=fmodel(mcon+1:mcon+mfit)-fi(1:mfit) ;
          f(1:mcon)=fmodel(1:mcon) ;
          if jacgen == 1
            if mprtim ~= 0
              monon(2) ;
            end
            [aa,ifail] = feval(fcn,x,'jacobian',par) ;
            if mprtim ~= 0
              monoff(2) ;
            end
          else
            if mprtim ~= 0
              monon(2) ;
            end
            if jacgen == 3
              [aa,eta,nfcnj,ifail] = ncjcf(fcn, par, n, m, x, fmodel, xw, eta, ...
                                           etamin, etamax, etadif, conv, nfcnj) ;
            end
            if jacgen == 2 
              [aa,nfcnj,ifail] = ncjac(fcn, par, n, m, x, fmodel, xw, ...
                                       ajdel, ajmin, nfcnj) ;
            end
            if mprtim ~= 0
              monoff(2) ;
            end
          end
          njac = njac + 1 ;
%     Exit, If ...
          if jacgen == 1  &&  ifail < 0
            ierr = 183 ;
            fprintf(luerr,'\n\n %s\n %s\n', ...
                          'Computation of the statistical analysis skipped', ...
                          'since JAC failed for the final Gauss-Newton iterate') ;
            return ;
          end
          if jacgen ~= 1  &&  ifail ~= 0
            ierr = 182 ;
            fprintf(luerr,'\n\n %s\n %s\n', ...
                          'Computation of the statistical analysis skipped', ...
                          'since FCN failed for the final Gauss-Newton iterate') ;
            return ;
          end
          a(1:m2,1:n) = aa(1:m2,1:n) ;
          wk.niwla(1) = 0 ;
          cond1 = cond ;
          mconh = mcon ;
          irank = minmn ;
          if mprtim ~= 0
            monon(3) ;
          end
          [a,d,cond1,qa,pivot,irank,wk,ifail] = ncfact(n,m1,mconh,n,a,d,pivot,cond1,irank,iopt,wk) ;
%           %%%%%%%%%%%%Susanna: nice print of ordered subcondition numbers
%           if mprmon > 0
%               fprintf(lumon,'\n\n Parameters ordered by increasing subcondition number:\n\n') ;
%               fprintf(lumon,'Parameter \t subcondition\n\n') ;
%               for j=1:irank
%                 fprintf(lumon,' %3i \t %10.2e \n',pivot(j),abs(d(1)/d(j))) ;  
%                 %fprintf(lumon,'%10.2e',vcv(i,j)) ;
%               end
%               fprintf(lumon,'\n') ;
%           end
          %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
          if mprtim ~= 0
             monoff(3) ;
          end
          if ifail ~= 0
            ierr = 180 ;
            fprintf(luerr,'\n\n %s\n %s\n', ...
                          'Computation of the statistical analysis skipped', ...
                          'since DECCON failed for the final Gauss-Newton iterate') ;
            return ;
          end
%         if mfit-(irank-mcon) > 0  !! stacon doesn't work with irank < n !!
          if mfit-(n-mcon) > 0
            irankc = wk.niwla(2) ;
            [ierr,t3,vcv,wk.rinf,xl,xr,wk.sigma2] = stacon (n,mcon,m,mfit,n,mcon,fi,x, ...
                         fmodel(mcon+1:m),a,irankc,irank,d,pivot,qa,mprmon,lumon) ;
          else
            fprintf(lumon,' \n\n Statistical analysis not available for not overdetermined system!\n\n') ;
          end
        end

function xw = ncscal(n,x,xa,xscal,iscal,qinisc)
%     ------------------------------------------------------------
%
%*    Summary :
%    
%     S C A L E : To be used in connection with NLSCON .
%       Computation of the internal scaling vector XW used for the
%       Jacobian matrix, the iterate vector and it's related
%       vectors - especially for the solution of the linear system
%       and the computations of norms to avoid numerical overflow.
%
%*    Input parameters
%     ================
%
%     N         Int     Number of unknowns
%     X(N)      Dble    Current iterate
%     XA(N)     Dble    Previous iterate
%     XSCAL(N)  Dble    User scaling passed from parameter XSCAL
%                       of interface routine NLSCON
%     ISCAL     Int     Option ISCAL passed from IOPT-field
%                       (for details see description of IOPT-fields)
%     QINISC    Logical = 1  : Initial scaling
%                       = 0 : Subsequent scaling
%
%*    Output parameters
%     =================
%
%     XW(N)     Dble   Scaling vector computed by this routine
%                      All components must be positive. The follow-
%                      ing relationship between the original vector
%                      X and the scaled vector XSCAL holds:
%                      XSCAL(I) = X(I)/XW(I) for I=1,...N
%
%*    Functions called: D1MACH
%
%*    Machine constants used
%     ======================
%
%     SMALL
%
%     ------------------------------------------------------------
%*    End Prologue
      small  = d1mach(6) ;
%*    Begin
      if iscal == 1
        xw = xscal ;
      else
        xw = zeros(n,1) ;
        for l1=1:n
          xw(l1)=max(xscal(l1),max((abs(x(l1))+abs(xa(l1)))*0.5,small)) ;
        end
      end
%     End of function NCSCAL

function [a,fw] = ncscrf(n,m,mcon,a,fw)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     S C R O W F : Row Scaling of a (M,N)-matrix in full storage
%                   mode
%
%*    Input parameters (* marks inout parameters)
%     ===========================================
%
%       N           Int    Number of columns of the matrix
%       M           Int    Number of rows of the matrix -
%                          leading dimension of A
%       MCON        Int    Number of rows to be automatically rescaled
%                          (the rows 1 to MCON will be rescaled)
%     * A(M,N)      Dble   Matrix to be scaled
%     * FW(M)       Dble   Scaling vector 
%
%*    Output parameters
%     =================
%
%       A(M,N)      Dble   The scaled matrix
%       FW(M)       Dble   Row scaling factors - FW(i) contains
%                          the factor by which the i-th row of A
%                          has been multiplied -
%                          unaltered in components MCON+1,...,M
%
%     ------------------------------------------------------------
%*    End Prologue
      one=1.0 ;
      zero=0.0 ;
%*    Begin
      for k=1:mcon
        s1 = max(abs(a(k,1:n))) ;
        if s1 > zero
          s1=one/s1 ;
          fw(k)=s1 ;
          a(k,1:n) = a(k,1:n)*s1 ;
        else
          fw(k)=one ;
        end
      end
      for k=mcon+1:m
        s1 = fw(k) ;
        a(k,1:n) = a(k,1:n)*s1 ;
      end
%     End of function N1SCRF

function [a,d,cond,ainv,pivot,irank,wk,ifail] = ncfact(n,m,mcon,ldainv,a,d,pivot,cond,irank,iopt,wk)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     F A C T : Call linear algebra subprogram for factorization of
%               a (N,N)-matrix with rank decision and casual compu-
%               tation of the rank deficient pseudo-inverse matrix
%
%*    Input parameters (* marks inout parameters)
%     ===========================================
%
%     N             Int    Number of parameters to be estimated
%     M             Int    Number of observations + equal. constraints
%     MCON          Int    Number of equality constraints
%     LDAINV        Int    Leading dimension of the matrix array AINV
%   * A(M,N)        Dble   Matrix storage.
%   * D(N)          Dble   In case of rank-reduction, D holds the output
%                          from a previous call of NCFACT
%   * PIVOT(N)      Int    In case of rank-reduction, PIVOT holds the output
%                          from a previous call of NCFACT
%   * COND          Dble   Maximum permitted subcondition for the prescribed rank
%   * IRANK         Int    Prescribed maximum rank of the matrix A
%     IOPT          Struct Option vector passed from NLSCON
%   * WK            Struct Workspace structure array
%
%*    Output parameters
%     =================
%
%   * A(M,N)         Dble   The QR-decomposition of the input matrix A.
%   * D(N)           Dble   The diagonal elements of the QR-decomposition
%   * COND           Dble   Estimated subcondition of A
%     AINV(LDAINV,N) Dble   If matrix A is rank deficient, this array
%                           holds the pseudo-inverse of A
%   * PIVOT(N)       Int    The table of column exchanges from pivoting 
%   * IRANK          Int    Actual rank of matrix A, but not greater than
%                             the input value.
%   * WK             Struct Workspace structure array
%     IFAIL          Int    Error indicator returned by this routine:
%                           = 0 matrix decomposition successfull
%                           = 2 DECCON failed due to a very bad conditioned
%                               matrix A.
%
%*    functions called:  DECCON
%
%     ------------------------------------------------------------
%*    Begin
      mprerr = iopt.mprerr ;
      luerr  = iopt.luerr ;
      irepet = -wk.niwla(1) ;
      if irepet == 0
        wk.niwla(2) = mcon ;
      end
      [a,wk.niwla(2),irank,d,pivot,cond,ainv,v,ifail] = ...
        deccon(a,m,n,mcon,m,n,wk.niwla(2),irank,cond,irepet,d,pivot) ;
      if ifail == -2  &&  mprerr > 0
        fprintf(luerr,'\n deccon failed to compute rank-deficient qr-decomposition\n') ;
      end
      if irank ~= 0
        irankc = wk.niwla(2) ;
        wk.nrwla(1) = abs(d(irankc+1)) ;
        wk.nrwla(2) = v(1) ;
      else
        cond = 1.0 ;
        wk.nrwla = 0 ;
      end

function [b,z,ifail] = ncfit(n,m,mcon,ldainv,a,ainv,b,irank,d,pivot,wk)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     F I T : Call linear algebra subprogram for (least squares) 
%             solution of the linear system A*Z = B
%
%*    Parameters
%     ==========
%
%     N,M,MCON,LDAINV,A,AINV,IRANK,IOPT,IFAIL,WK :
%                        See description for function NCFACT.          
%     B          Dble    In:  Right hand side of the linear system
%                        Out: Rhs. transformed to the upper trian-
%                             gular part of the linear system
%     Z          Dble    Out: Solution of the linear system
%
%     functions called: SOLCON
%
%     ------------------------------------------------------------
%*    Begin
      irepet = -wk.niwla(1) ;
      [z,b] = solcon(a,m,n,mcon,m,n,b,wk.niwla(2),irank,d,pivot,irepet,ainv) ;
      ifail = 0 ;

function [dxq,conv,sumx,dlevf] = nclvls(n,m,dx1,xw,f,mprmon)
%
%     ------------------------------------------------------------
%
%*    Summary :
%
%     L E V E L S : To be used in connection with NLSCON .
%     provides descaled solution, error norm and level functions
%
%*    Input parameters (* marks inout parameters)
%     ===========================================
%
%       N              Int    Number of parameters to be estimated
%       M              Int    Number of measurements plus constraints 
%       DX1(N)         Dble   array containing the scaled Newton
%                             correction
%       XW(N)          Dble   Array containing the scaling values
%       F(N)           Dble   Array containing the residuum
%       MPRMON         Int    Print information parameter (see
%                             driver routine NLSCON )
%
%*    Output parameters
%     =================
%
%       DXQ(N)         Dble   Array containing the descaled Newton
%                             correction
%       CONV           Dble   Scaled maximum norm of the Newton
%                             correction
%       SUMX           Dble   Scaled natural level function value
%       DLEVF          Dble   Standard level function value
%
%     ------------------------------------------------------------
%*    Begin
%     ------------------------------------------------------------
%     1 Descaling of solution DX1 ( stored to DXQ )
      dxq = dx1 .* xw ;
%     ------------------------------------------------------------
%     2 Evaluation of scaled natural level function SUMX and
%       scaled maximum error norm CONV
      conv = max(abs(dx1)) ;
      sumx = sum(dx1.^2) ;
%     ------------------------------------------------------------
%     3 Evaluation of (scaled) standard level function DLEVF
      dlevf = sqrt(sum(f.^2)/m) ;
%     End of function N1LVLS

function [a,nfcnnew,ifail] = ncjac(fcn, par, n, m, x, fx, yscal, ...
                                    ajdel, ajmin, nfcn)
%  ---------------------------------------------------------------------
%
%* Title
%
%  Evaluation of a dense Jacobian matrix using finite difference
%  approximation adapted for use in nonlinear systems solver NLEQ1
%
%* Environment       Matlab 6.0
%                    Sparc, Solaris
%* Latest Revision   August 2001
%
%
%* Parameter list description
%  --------------------------
%
%* External functions (to be supplied by the user)
%  -------------------------------------------------
%
%     [FOUT,FAIL] =FCN(X,FLAG,PAR) String  Name of problem-function
%       For a description refer to the driver function nleq1.
%
%* Input parameters (* marks inout parameters)
%  ----------------
%
%  FCN        String  Name of the problem-function
%  PAR        AnyType User-parameter, to be passed to the problem-function 
%  N          Int     Number of columns of the Jacobian
%  M          Int     Number of rows of the Jacobian
%  X(N)       Dble    Array containing the current scaled
%                     iterate
%  FX(N)      Dble    Array containing FCN(X)
%  YSCAL(N)   Dble    Array containing the scaling factors
%  AJDEL      Dble    Perturbation of component k: abs(Y(k))*AJDEL
%  AJMIN      Dble    Minimum perturbation is AJMIN*AJDEL
%  NFCN       Int  *  FCN - evaluation count
%
%* Output parameters (* marks inout parameters)
%  -----------------
%
%  A(M,N)     Dble    Array to contain the approximated
%                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
%  NFCNNEW    Int  *  FCN - evaluation count adjusted
%  IFAIL      Int     Return code non-zero if Jacobian could not
%                     be computed.
%  ---------------------------------------------------------------------
%
%* Begin
%
      ifail = 0 ;
      nfcnnew = nfcn ;
      for k = 1:n
         w = x(k) ;
         
         su = sign(x(k)) ;
         if su == 0
           su = 1 ;
         end
         u = max(max(abs(x(k)),ajmin),yscal(k))*ajdel*su ;
         x(k) = w + u ;
%
         [fu,ifail] = feval(fcn,x,'',par) ;
         nfcnnew = nfcnnew + 1 ;
         if ifail  ~=  0
           break ;
         end
%
         x(k) = w ;
         a(1:m,k) = (fu-fx) / u ;
      end
%* End of NCJAC

function [a,eta,nfcnnew,ifail] = ncjcf(fcn,par,n,m,x,fx,yscal, ...
                                       eta,etamin,etamax,etadif,conv,nfcn)
%
%  ---------------------------------------------------------------------
%
%* Title
%
%  Approximation of dense Jacobian matrix for nonlinear systems solver
%  NLEQ1 with feed-back control of discretization and rounding errors
%
%* Environment       Matlab 6.0
%                    Sparc, Solaris
%* Latest Revision   August 2001
%
%
%* Parameter list description
%  --------------------------
%
%* External functions (to be supplied by the user)
%  -------------------------------------------------
%
%     [FOUT,FAIL] =FCN(X,FLAG,PAR) String  Name of problem-function
%       For a description refer to the driver function nleq1.
%
%* Input parameters (* marks inout parameters)
%  ----------------
%
%  FCN        String  Name of the problem-function
%  PAR        AnyType User-parameter, to be passed to the problem-function 
%  N          Int     Number of columns of the Jacobian
%  M          Int     Number of rows of the Jacobian
%  X(N)       Dble    Array containing the current scaled
%                     iterate
%  FX(N)      Dble    Array containing FCN(X)
%  YSCAL(N)   Dble    Array containing the scaling factors
%  ETA(N)     Dble *  Array containing the scaled denominator
%                     differences
%  ETAMIN     Dble    Minimum allowed scaled denominator
%  ETAMAX     Dble    Maximum allowed scaled denominator
%  ETADIF     Dble    sqrt (1.1*EPMACH)
%                     EPMACH = machine precision
%  CONV       Dble    Maximum norm of last (unrelaxed) Newton correction
%  NFCN       Int  *  FCN - evaluation count
%
%* Output parameters (* marks inout parameters)
%  -----------------
%
%  A(M,N)     Dble    Array to contain the approximated
%                     Jacobian matrix ( dF(i)/dx(j)in A(i,j))
%  ETA(N)     Dble *  Scaled denominator differences adjusted
%  NFCNNEW    Int  *  FCN - evaluation count adjusted
%  IFAIL      Int     Return code non-zero if Jacobian could not
%                     be computed.
%
%* Constants
%  ---------
      small2 = 0.1 ;
%
%  ---------------------------------------------------------------------
%* Begin
%
      ifail = 0 ;
      nfcnnew = nfcn ;
      a = zeros(m,n) ;
      for k = 1:n
         is = 0 ;
         qfine = 0 ;
         qexit = 0 ;
         while ~ qfine
            w = x(k) ;
            su = sign(x(k)) ;
            if su == 0
              su = 1 ;
            end
            u = eta(k)*yscal(k)*su ;
            x(k) = w + u ;
            [fu,ifail] = feval(fcn,x,'',par) ;
            nfcnnew = nfcnnew + 1 ;
%           exit, if ...
            if ifail  ~=  0
              qexit = 1 ;
              break ;
            end
            x(k) = w ;
            sumd = 0.0 ;
            for i = 1:m
               hg = max (abs(fx(i)), abs(fu(i))) ;
               fhi = fu(i) - fx(i) ;
               if hg  ~=  0.0
                 sumd = sumd + (fhi/hg)^2 ;
               end
               a(i,k) = fhi / u ;
            end
            sumd = sqrt (sumd/m) ;
            qfine = 1 ;
            if sumd  ~=  0.0  &&  is  ==  0
               eta(k) = min( etamax , max(etamin , sqrt(etadif/sumd)*eta(k)) ) ;
               is = 1 ;
               qfine = conv  <  small2  ||  sumd  >=  etamin ;
            end
         end
         if qexit
           break ;
         end
      end
%* End of function NCJCF

function a = ncrnk1(n,m,xw,dx,f,fa,a,fca)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     R A N K 1 : To be used in connection with NLSCON .
%     provides Rank-1 updates of Jacobian matrix A
%
%*    Input parameters
%     ================
%
%       N          Int    Number of columns of the Jacobian
%       M          Int    Number of rows of the Jacobian
%       XW(N)      Dble   Array containing the scaling factors
%       DX(N)      Dble   Last (unrelaxed) Gauss-Newton correction
%       F(M)       Dble   FCN(x(k)),  with x(k)= current iterate
%       FA(M)      Dble   FCN(x(k-1)),  with x(k-1)= previous
%                         iterate
%     * A(M,N)     Dble   The original Jacobi matrix
%       FCA        Dble   Previous damping factor
%
%*    Output parameters:
%     ==================
%
%       A(M,N)     Dble   The rank-1 updated Jacobian matrix 
%                         ( dF(i)/dx(j)in A(i,j))
%
%     ------------------------------------------------------------
      zero=0.0 ;
      one=1.0 ;
%*    Begin
        dxj = dx ./ xw ;
        dnm = sum(dxj(1:n).^2) ;
        dxj = dxj ./ xw ;
        dnm = dnm*fca ;
        if dnm ~= zero
          s1 = fca-one ;
          dxf = f(1:m) + fa(1:m)*s1 ;
          for k=1:n
            s1 = dxj(k)/dnm ;
            a(1:m,k) = a(1:m,k)+dxf(1:m)*s1 ;
          end
        end
%       End of function NCRNK1

function del = ncprjn(n,irank,u,d,qe,pivot)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     P R J C T N :
%     To be used in connection with either DECOMP/SOLVE or 
%     DECCON/SOLCON .
%     Provides the projection to the appropriate subspace in case
%     of rank - reduction
%
%*    Input parameters (* marks inout parameters)
%     ===========================================
%
%       N              Int    Number of parameters to be estimated
%       IRANK                 Pseudo rank of decomposed Jacobian
%                             matrix
%       U(N)           Dble   Scaled Gauss-Newton correction
%       D(N)           Dble   Diagonal elements of upper
%                             triangular matrix
%       QE(N,N)        Dble   Part of pseudoinverse Jacobian
%                             matrix ( see QA of DECCON )
%       PIVOT(N)       Dble   Pivot vector resulting from matrix
%                             decomposition (DECCON)
%
%*    Output parameters
%     =================
%
%       DEL            Dble   Defekt
%
%     ------------------------------------------------------------
%*    End Prologue
      zero = 0.0 ;
%*    Begin
      v = zeros(n,1) ;
      for i=1:n
        v(i)=u(pivot(i)) ;
      end
      irk1 = irank+1 ;
      del = zero ;
      for i=irk1:n
        s =( v(i) - sum(qe(1:i-1,i).*v(1:i-1)) )/d(i) ;
        del = s*s+del ;
        v(i)=s ;
      end
%     End of function NCPRJN

function dummy = ncprv1(dlevf,dlevx,fc,niter,new,irank,mprmon,lumon,qmixio)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     N C P R V 1 : Printing of intermediate values (Type 1 routine)
%
%     Parameters
%     ==========
%
%     DLEVF, DLEVX   See descr. of internal double variables of NCINT
%     FC,NITER,NEW,IRANK,MPRMON,LUMON
%                  See parameter descr. of function NCINT
%     QMIXIO Logical  = 1 , if LUMON == LUSOL
%                     = 0 , if LUMON ~= LUSOL
%
%     ------------------------------------------------------------
%*    End Prologue
%     Print Standard - and natural level
      if qmixio
        fprintf(lumon,'  **************************************************************************\n') ;
        if mprmon >= 3
          fprintf(lumon,'        It       Normf               Normx                   New      Rank\n') ;
        end
        if mprmon == 2
          fprintf(lumon,'        It       Normf               Normx       Damp.Fct.   New      Rank\n') ;
        end
      end
      if mprmon >= 3 || niter == 0
        fprintf(lumon,'      %4i     %14.7e      %10.3e               %2i      %4i\n', ...
                      niter,dlevf,dlevx,new,irank) ;
      end
      if mprmon == 2 && niter ~= 0
        fprintf(lumon,'      %4i     %14.7e      %10.3e      %5.3f    %2i      %4i\n', ...
                      niter,dlevf,dlevx,fc,new,irank) ;
      end
      if qmixio
        fprintf(lumon,'  **************************************************************************\n') ;
      end
      dummy = 1 ;
%     End of function NCPRV1

%
function dummy = ncprv2(dlevf,dlevx,fc,niter,mprmon,lumon,qmixio,cmark)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     N C P R V 2 : Printing of intermediate values (Type 2 routine)
%
%*    Parameters
%     ==========
%
%     DLEVF,DLEVX    See descr. of internal double variables of NCINT
%     FC,NITER,MPRMON,LUMON
%                  See parameter descr. of function NCINT
%     QMIXIO Logical  = 1 , if LUMON == LUSOL
%                     = 0 , if LUMON ~= LUSOL
%     CMARK Char*1    Marker character to be printed before DLEVX
%
%     ------------------------------------------------------------
%     Print Standard - and natural level, and damping factor
      if qmixio
        fprintf(lumon,'  **************************************************************************\n') ;
        fprintf(lumon,'        It       Normf               Normx       Damp.Fct.\n') ;
      end
      fprintf(lumon,'      %4i     %14.7e    %1s %10.3e      %5.3f\n',niter,dlevf,cmark,dlevx,fc) ;
      if qmixio
        fprintf(lumon,'  **************************************************************************\n') ;
      end
      dummy = 1 ;
%     End of function NCPRV2

function dummy = ncsout(n,mfit,x,rq,mode,iopt,wk,mprint,luout)

%    global ndat % index of to be fitted data
%     ------------------------------------------------------------
%
%*    Summary :
%
%     S O L O U T : Printing of iterate (user customizable routine)
%
%*    Input parameters
%     ================
%
%     N         Int Number of equations/unknowns
%     X(N)     Dble iterate vector
%     RQ(MFIT) Dble   Linear residuum (without zero components
%                     corresponding to the equality constraints)
%     MODE          =1 This routine is called before the first
%                      Newton iteration step
%                   =2 This routine is called with an intermedi-
%                      ate iterate X(N)
%                   =3 This is the last call with the solution
%                      vector X(N)
%                   =4 This is the last call with the final, but
%                      not solution vector X(N)
%     IOPT   Struct The option structure as passed to the driver
%                   routine(certain fields may be used
%                   for user options)
%     MPRINT    Int Solution print level 
%                   (see description of IOPT-field MPRINT)
%     LUOUT     Int the solution print unit 
%                   (see description of see IOPT-field LUSOL)
%
%
%*    Workspace parameters
%     ====================
%
%     WK    see description in driver routine
%
%     ------------------------------------------------------------
%*    Begin
      qnorm = 1 ;
      if qnorm
         if mode == 1
           fprintf(luout,'%s\n%s%5i\n\n%s\n','  Start data:','  N =',n, ...
                  '  Format: iteration-number, (x(i),i=1,...N) , Normf , Normx ') ;
           fprintf(luout,'%s\n','  Initial data:') ;
         elseif mode == 3
           fprintf(luout,'%s\n','  Solution data:') ;
         elseif mode == 4
           fprintf(luout,'%s\n','  Final data:') ;
         end
         fprintf(luout,' %5i\n',wk.niter) ;
         l2 = 0 ;
         for l1=1:n
           fprintf(luout,'%18.10e ',x(l1)) ;
           l2 =l2+1 ;
           if l2 == 3
             fprintf(luout,'%1s\n',' ') ;
             l2 = 0 ;
           end
         end
         fprintf(luout,'%18.10e %18.10e \n', wk.dlevf, sqrt(wk.sumx/n) ) ;
         if mode == 2 && mprint >= 3
           if mfit ~= 0
             fprintf(luout,'\n    Residuum for current iteration :\n            '); % added Semicolon (B.R. 31.07.17)
             l2 = 0 ;
             for l1=1:mfit
               fprintf(luout,'%18.10e ',rq(l1)) ;
               l2 =l2+1 ;
               if l2 == 3
                 fprintf(luout,'\n            ') ;
                 l2 = 0 ;
               end
             end
           end
         end
         if mode == 1 && mprint >= 2
           fprintf(luout,'%s\n','  Intermediate data:') ;
         elseif mode >= 3
           fprintf(luout,'%s\n','  End data:') ;
         end
      end
      dummy = 1 ;
%     End of function NCSOUT

function norm = wnorm(n,z,xw)
%     ------------------------------------------------------------
%
%*    Summary :
%
%     W N O R M : Return the norm to be used in exit (termination)
%                 criteria
%
%*    Input parameters
%     ================
%
%     N         Int Number of equations/unknowns
%     Z(N)     Dble  The vector, of which the norm is to be computed
%     XW(N)    Dble  The scaling values of Z(N)
%
%*    Output
%     ======
%
%     WNORM(N,Z,XW)  Dble  The mean square root norm of Z(N) subject
%                          to the scaling values in XW(N):
%                          = Sqrt( Sum(1,...N)((Z(I)/XW(I))**2) / N )
%
%     ------------------------------------------------------------
%*    Begin
      norm = sqrt( sum( (z./xw).^2 ) / n ) ;
%     End of function WNORM

function [ierr,res,vcv,rinv,xl,xr,sigma2] = stacon (ndcl,mcodcl,mdcl,mfit,n,mcon,y,x, ...
                                        ymodel,a,irankc,irank,d,ipiv,ah,mprmon,lumon)
%----------------------------------------------------------------------
%
%  Statistical analysis of constrained linear least squares estimates.
%  Computation of covariance matrix, correlation coefficients  and
%  confidence intervals for final parameter estimates.
%
%  To be used in connection with DECCON and SOLCON.
%
%----------------------------------------------------------------------
%
%  Date of latest change:  July 23, 2001
%  By: U. Nowak, L. Weimann
%
%***********************************************************************
%
%   Input parameters
%   ----------------
%
%      NDCL          Int   Declared number of columns of a (see below)
%      MCODCL        Int
%      MDCL          Int   Declared number of rows of a (see below)
%      MFIT          Int   Number of measurements (observations)
%      N             Int   Number of parameters
%      MCON          Int   Number of equality constraints
%      Y(MFIT)       Dble  The vector of measurements
%      X(N)          Dble  Final (estimated) parameters
%      YMODEL(MFIT)  Dble  The model function vector evaluated for
%                          the final parameter estimates.
%      A(MDCL,NDCL)  Dble  Matrix of the model (unscaled Jacobian)
%      IRANKC        Int   Rank of the matrix constrains part
%      IRANK         Int   Rank of the matrix least squares part
%      D(N)          Dble  Diagonal elements of decomposed matrix A
%      IPIV(N)       Int   Column interchanges performed by 'DECCON'
%      AH(NDCL,NDCL) Dble  The rank deficient pseudo inverse computed
%                          by DECCON (if IRANK < N)
%      V(N)          Dble  Work array
%      MPRMON        Int   The printing level:
%                          =0 : no printing will be done
%                          =1 : Information will be printed
%      LUMON         Int   The print unit number
%
%  Output parameters
%  -----------------
%
%      IERR          Int    Error indicator
%                           =0 : no error occured
%      RES(MFIT)     Dble   The residuum YMODEL-Y
%      VCV(N,N)      Dble   Correlation matrix
%      RINV(NDCL,N)  Dble   The matrix of the correlation coefficients
%                           (lower triangle only) - only, if MPRMON > 0
%      XL(N)         Dble   The left bounds of the confidence intervals
%                           associated to the final parameter estimate
%                           vector X(N)
%      XR(N)         Dble   The right bounds of the confidence intervals
%                           associated to the final parameter estimate
%                           vector X(N)
%      SIGMA2        Dble   Estimated variance of residual
%
%***********************************************************************
%
      zero=0.0 ; one=1.0 ; hundre=100.0 ;
%
%  FISH15: Array containing upper 5% values of fisher(1,l)-distribution
%  (L=1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,
%  26,27,28,29,30,32,34,36,38,40,42,44,46,48,50,60,70,80,90,100,200,300)
%
      fish15 = [161.0,18.51,10.13,7.71,6.61,5.99,5.59,         ...
                 5.32,5.12,4.96,4.84,4.75,4.67,4.6,4.54,4.49,  ...
                 4.45,4.41,4.38,4.35,4.32,4.3,4.28,4.26,4.24,  ...
                 4.23,4.21,4.2,4.18,4.17,4.15,4.13,4.11,4.1,   ...
                 4.08,4.07,4.06,4.05,4.04,4.03,4.0,3.98,3.96,  ...
                 3.95,3.94,3.89,3.88] ;
      ierr = 0 ;
      vcv = zeros(n,n) ;  v = zeros(1,n) ;
%
%  Estimate variance and standard deviation of residual
% 
      if mprmon > 0
        fprintf(lumon,'\n\n\n          Under the assumptions of the classical linear model:\n\n\n\n') ;
      end
%
      res = ymodel(1:mfit)-y(1:mfit) ;
      sum1 = sum(res(1:mfit).^2) ;
      sum1=sum1/(mfit-(irank-mcon)) ;
%
      sigma2=sum1 ;
      if mprmon > 0
        fprintf(lumon,'\n%s\n%s\n', ...
        '   Best unbiased estimate of variance and standard deviation of residuals:' , ...
        '   -----------------------------------------------------------------------' ) ;
        fprintf(lumon,'\n   sigma2 =%10.3e     sigma =%10.3e\n',sigma2,sqrt(sigma2) ) ;
      end
%
% Computation of covariance matrix
%---------------------------------
%
      if n ~= 1
%
%  Apply pseudo-inverse to cholesky decomposition of
%  variance covariance matrix of data errors. We assume, that
%  this matrix is sigma*I
%
        m=mfit+mcon ;
        for j=1:n
           xl = zeros(n,1) ;
           xl(j)=one ;
           if j <= mcon
             xl(j)=zero ;
           end
           kred=-1 ;
           [xr,xl] = solcon(a,mdcl,ndcl,mcon,m,n,xl,irankc,irank,d,ipiv,kred,ah) ;
           rinv(1:n,j) = xr ;
        end
%
%  Product RINV*RINV^T
        for i=1:n
           for j=1:n
              vcv(i,j)=sigma2*sum(rinv(i,1:n).*rinv(j,1:n)) ;
           end
        end
      else
        vcv(1,1)=sigma2/(d(1)*d(1)) ;
      end
%
%  Pretty print out of the covariance matrix
%  -----------------------------------------
      if mprmon > 0
        fprintf(lumon,'\n\n\n   Covariance matrix of parameters\n   -----------------\n') ;
        nh2=0 ;
        while nh2 ~= n
          nh1=nh2+1 ;
          nh2=nh1+6 ;
          if nh2 > n
            nh2=n ;
          end
          if mprmon > 0
            fprintf(lumon,'\n') ;
            for l=nh1:nh2
              fprintf(lumon,'%10i',l) ;
            end
            fprintf(lumon,'\n') ;
          end
          for i=nh1:n
            ih=i ;
            if ih > nh2
              ih=nh2 ;
            end
            if mprmon > 0
              fprintf(lumon,' %3i',i) ;
              for j=nh1:ih
                fprintf(lumon,'%10.2e',vcv(i,j)) ;
              end
              fprintf(lumon,'\n') ;
            end
          end
        end
      end
%
%  Computation and pretty printout of the correlation coefficients
%  ---------------------------------------------------------------
      for i=1:n
        v(i)=sqrt(vcv(i,i)) ;
        if v(i) ~= zero ;
          rinv(i,i)=one ;
        else
          rinv(i,i)=zero ;
        end
      end
%
      if mprmon > 0
        if n ~= 1
          nm1=n-1 ;
          for j=1:nm1
            jp1=j+1 ;
            if v(j) == zero
              rinv(jp1:n,j) = zeros(n-jp1+1,1) ;
            else
              for i=jp1:n
                if v(i) ~= zero
                  rinv(i,j)=vcv(i,j)/(v(i)*v(j)) ;
                else
                  rinv(i,j)=zero ;
                end
              end
            end
          end
%
          fprintf(lumon,'\n\n\n   Correlation coefficients\n   ------------------------\n') ;
          nh2=0 ;
          while nh2 ~= n
            nh1=nh2+1 ;
            nh2=nh1+9 ;
            if nh2 > n
              nh2=n ;
            end
            fprintf(lumon,'\n   ') ;
            for l=nh1:nh2
              fprintf(lumon,'%7i',l) ;
            end
            fprintf(lumon,'\n') ;
            for i=nh1:n
              ih=i ;
              if ih > nh2
                ih=nh2 ;
              end
              fprintf(lumon,' %3i',i) ;
              for j=nh1:ih
                fprintf(lumon,'%7.2f',rinv(i,j)) ;
              end
              fprintf(lumon,'\n') ;
            end
          end
        end
      end
%
%  Standard error in parameters
%  ----------------------------
      if mprmon > 0
        fprintf(lumon,'\n\n\n   Standard deviation of parameters\n') ;
        fprintf(lumon,      '   --------------------------------\n') ;
        fprintf(lumon,'     No.  Estimate           sigma(X)\n\n') ;
        for i=1:n
          if x(i) ~= zero
            xr(i)=abs(hundre*v(i)/x(i)) ;
          else
            xr(i)=abs(hundre*v(i)) ;
          end
          fprintf(lumon,'   %4i  %10.3e   +/-  %10.3e     =%8.2f %%\n',i,x(i),v(i),xr(i)) ;
        end
      end
%
      l=mfit-(n-mcon) ;
      if l <= 30 
        lh = l ;
      elseif l <= 50
        lh = 15+fix(l/2) ;
      elseif l <= 100
        lh = 35+fix(l/10)  ;
      elseif l < 300
        lh = 45+fix(l/200) ;
      else
        lh = 47 ;
      end
%
%  Associated confidence intervals
%  -------------------------------
      if mprmon > 0
        fprintf(lumon,'\n\n\n   independent confidence intervals\n   --------------------------------\n') ;
      end
      fnma=fish15(lh) ;
      h=sqrt(fnma) ;
      if mprmon > 0
        fprintf(lumon,'      (on 95%%-probability level using F-distribution  F(alfa,1,m-n)=%6.2f)\n',fnma) ;
      end
      for i=1:n
        xr(i)=h*v(i) ;
        v(i)=hundre*xr(i)/x(i) ;
        xl(i)=x(i)-xr(i) ;
        xr(i)=x(i)+xr(i) ;
        if mprmon > 0
          fprintf(lumon,'   %4i  ( %10.3e , %10.3e )\n',i,xl(i),xr(i)) ;
        end
      end

%*    End package
%
%*    Group  Linear Solver functions (Code DECCON/SOLCON)
%

function [a,irankc,irank,d,pivot,cond,ah,v,ierr] = ...
         deccon(a,nrow,ncol,mcon,m,n,irankc,irank,cond,kred,d,pivot)
%     ------------------------------------------------------------
%
%*  Title
%
%*    Deccon - Constrained Least Squares QR-Decomposition
%
%*  Written by        P. Deuflhard, U. Nowak, L. Weimann 
%*  Purpose           Solution of least squares problems, optionally
%                     with equality constraints.
%*  Method            Constrained Least Squares QR-Decomposition
%                     (see references below)
%*  Category          D9b1. -  Singular, overdetermined or
%                              underdetermined systems of linear 
%                              equations, generalized inverses. 
%                              Constrained Least Squares solution
%*  Keywords          Linear Least Square Problems, constrained, 
%                     QR-decomposition, pseudo inverse.
%*  Version           1.3
%*  Revision          December 1993
%*  Latest Change     August 2006
%*  Library           CodeLib
%*  Code              Matlab 6.0
%*  Environment       Standard Matlab 6.0 environment on PC's,
%                     workstations and hosts.
%*  Copyright     (c) Konrad-Zuse-Zentrum fuer
%                     Informationstechnik Berlin (ZIB)
%                     Takustrasse 7, D-14195 Berlin-Dahlem
%                     phone : + 49/30/84185-0
%                     fax   : + 49/30/84185-125
%*  Contact           Lutz Weimann 
%                     ZIB, Numerical Analysis and Modelling
%                     phone: 0049+30+84185-185 ;
%                     e-mail:weimann@zib.de
%
%*    References:
%     ===========
%
%       /1/ P.Deuflhard, V.Apostolescu:
%           An underrelaxed Gauss-Newton method for equality
%           constrained nonlinear least squares problems.
%           Lecture Notes Control Inform. Sci. vol. 7, p.
%           22-32 (1978)
%       /2/ P.Deuflhard, W.Sautter:
%           On rank-deficient pseudoinverses.
%           J. Lin. Alg. Appl. vol. 29, p. 91-111 (1980)
%    
%*    Related Programs:     SOLCON
%
%  ---------------------------------------------------------------
%
%* Licence
%    You may use or modify this code for your own non commercial
%    purposes for an unlimited time. 
%    In any case you should not deliver this code without a special 
%    permission of ZIB.
%    In case you intend to use the code commercially, we oblige you
%    to sign an according licence agreement with ZIB.
%
%* Warranty 
%    This code has been tested up to a certain level. Defects and
%    weaknesses, which may be included in the code, do not establish
%    any warranties by ZIB. ZIB does not take over any liabilities
%    which may follow from acquisition or application of this code.
%
%* Software status 
%    This code is under care of ZIB and belongs to ZIB software class 1.
%
%     ------------------------------------------------------------
%
%*    Summary:
%     ========
%     Constrained QR-decomposition of (M,N)-system  with
%     computation of pseudoinverse in case of rank-defeciency .
%     First MCON rows belong to equality constraints.
%
%     ------------------------------------------------------------
%
%*    Parameters list description (* marks inout parameters)
%     ======================================================
%
%*    Input parameters
%     ================
%
%       A(NROW,NCOL) Dble   Array holding the (M,N)-Matrix to be 
%                           decomposed
%       NROW         Int    Declared number of rows of array A
%       NCOL         Int    Declared number of columns of array A and 
%                           rows and columns of array AH
%       MCON         Int    Number of equality constraints (MCON <= N)
%                           Internally reduced if equality constraints
%                           are linearly dependent
%       M            Int    Current number of rows of matrix A
%       N            Int    Current number of columns of matrix A
%     * IRANKC       Int    Prescribed maximum pseudo-rank of 
%                           constrained part of matrix A (IRANKC <= MCON)
%     * IRANK        Int    Prescribed maximum pseudo-rank of matrix A
%                           (IRANK <= N)
%     * COND         Dble   Permitted upper bound for the subcondition
%                           of the least squares part of A, .i.e.
%                           DABS(D(IRANKC+1)/D(IRANK))
%       KRED         Int    Type of operation
%                           >=0  Householder triangularization
%                                (build up pseudo-inverse,if IRANK < N)
%                           < 0  Reduction of pseudo-rank of matrix A, 
%                                skipping Householder triangularization,
%                                 build-up new pseudo-inverse
%       D(IRANK)      Dble   Diagonal elements of upper triangular matr. ;
%                            obtained from a previous call of deccon -
%                            if current call is for rank-reduction (KRED<0)
%       PIVOT(N)      Int    Index vector storing permutation of columns
%                            due to pivoting ;
%                            obtained from a previous call of deccon -
%                            if current call is for rank-reduction (KRED<0)
%
%*    Output parameters
%     =================
%
%       A(NROW,NCOL)  Dble   Array holding the (M,N)-output consisting
%                            of the transformed matrix in the upper 
%                            right triangle and the performed House-
%                            holder transf. in the lower left triangle.
%     * IRANKC        Int    New pseudo-rank of constrained part of
%                            matrix A, determined so that
%                            DABS(D(1)/D(IRANKC))<1/EPMACH
%     * IRANK         Int    New pseudo-rank of matrix A, determined
%                            so that DABS(D(IRANKC+1)/D(IRANK)) < COND
%       D(IRANK)      Dble   Diagonal elements of upper triangular matr.
%       PIVOT(N)      Int    Index vector storing permutation of columns
%                            due to pivoting
%     * COND          Dble   The sub-condition number belonging to the
%                            least squares part of A.
%                            (in case of rank reduction:
%                             sub-condition number which led to
%                             rank reduction)
%                            COND=0 indicates COND=infinity
%       AH(NCOL,NCOL) Dble   In case of rank-defect used to compute the
%                            pseudo-inverse (currently used will be an
%                            (N,N)-part of this array)
%       V(N)          Dble   V(1) holds on output the sub-condition
%                            number belonging to the constrained part
%                            of A.
%       IERR          Int    Error indicator:
%                            = 0 : DECCON computations are successfull.
%                            =-2 : Numerically negative diagonal element
%                                  encountered during computation of
%                                  pseudo inverse - due to extremely bad
%                                  conditioned Matrix A. DECCON is
%                                  unable to continue rank-reduction.
%
%*    functions called: D1MACH
%
%*    Machine constants used
%     ======================
%
%     EPMACH = relative machine precision
%
%     ------------------------------------------------------------
%*    Begin
%     --------------------------------------------------------------
%     1 Initialization

      zero=0.0 ;
      one=1.0 ;
      reduce=0.05 ;
      epmach = d1mach(3) ;
      if irank > n 
        irank = n ;
      end
      if irank > m
        irank = m ;
      end
      ah = zeros(n,n) ;
      t = 0.0 ;
%     --------------------------------------------------------------
%     1.1 Special case M=1 and N=1
      if m == 1 && n == 1
        pivot(1)=1 ;
        d(1)=a(1,1) ;
        cond = one ;
        return ;
      end
      if kred >= 0
%       ------------------------------------------------------------
%       1.1 Initialize pivot-array
        pivot(1:n)=1:n ;
%       ------------------------------------------------------------
%       2. Constrained Householder triangularization    
        jd = 1 ;
        iranc1 = irankc + 1 ;
        mh = mcon ;
        irankh = irankc ;
        idata = 0 ;
        if mh == 0
          irankh = irank ;
          mh = m ;
          idata = 1 ;
        end
        irk1 = irank ;
        for  k=1:irk1
          qloop1 = 1 ;
          while qloop1
            qloop1 = 0 ;
            level = 1 ;
            if k ~= n
              k1 = k+1 ;
              jd = 1 ;
              while jd == 1
                for j=k:n
                  d(j)=sum(a(k:mh,j).^2) ;
                end
%               ------------------------------------------------------
%               2.1 Column pivoting
                s1 = d(k) ;
                jj = k ;
                for l1=k:n
                  if d(l1) > s1
                    s1=d(l1) ;
                    jj = l1 ;
                  end
                end
                h = d(jj) ;
                if jd == 1
                  hmax = h/max(10.0,cond*reduce) ;
                end
                jd = 0 ;
                if h < hmax
                  jd = 1 ;
                end
              end
              if jj ~= k
%             ------------------------------------------------------
%             2.2 Column interchange
                i = pivot(k) ;
                pivot(k)=pivot(jj) ;
                pivot(jj)=i ;
                d(jj)=d(k) ;
                s1=a(1:m,jj) ;
                a(1:m,jj)=a(1:m,k) ;
                a(1:m,k)=s1 ;
              end
            end
            h = sum(a(k:mh,k).^2) ;
            t = sqrt(h) ;
%           ----------------------------------------------------------
%           2.3.0 A-priori test on pseudo-rank
            if  k == 1  ||  k == iranc1
              dd = t/cond ;
            end
            if t <= dd  ||  k > irankh
%             ------------------------------------------------------
%             2.3.1 Rank reduction
              irankh = k-1 ;
              if mh ~= mcon || idata == 1
                irank = irankh ;
                if irankc == irank
                  level = 4 ;
                else
                  level = 3 ;
                end
              else  
                irankc = irankh ;
                if irankc ~= mcon
                  mh = m ;
                  irankh = irank ;
                  jd = 1 ;
                  idata = 1 ;
                  qloop1 = 1 ;
                else
                  ierr = -999 ;
                  return ;
                end
              end
            end
          end
          if level == 1
%           ------------------------------------------------------
%           2.4 Householder step
            s = a(k,k) ;
            if s == 0.0
              t = -abs(t) ;
            else
              t = -abs(t)*sign(s) ;
            end
            d(k)=t ;
            a(k,k)=s-t ;
            if k ~= n
              t = one/(h-s*t) ;
              for j=k1:n
                s = sum(a(k:mh,k).*a(k:mh,j))*t ;
                s1 =-s ;
                if s ~= 0.0 
                  a(k:m,j) = a(k:m,j)+a(k:m,k)*s1 ;
                end
                d(j)=d(j)-a(k,j)^2 ;
              end
              if k == irankc
                mh = m ;
                jd = 1 ;
                irankh = irank ;
              end
              if k == irk1
                level = 3 ;
              end
            else
              level = 4 ;
            end
          end
%       Exit for-loop If ... 
          if level > 1
            break ;
          end
        end
%       enddo
      else
        k = -1 ;
        level = 3 ;
      end
%     --------------------------------------------------------------
%     3 Rank-deficient pseudo-inverse
      if level == 3
        irk1 = irank+1 ;
        for j=irk1:n
          for ii=1:irank
            i = irk1-ii ;
            s = a(i,j) ;
            if ii ~= 1
              s = s - sum(a(i,i1:irank).*v(i1:irank)) ;
            end
            i1 = i ;
            v(i)=s/d(i) ;
            ah(i,j)=v(i) ;
          end
          for i=irk1:j
            s = sum(ah(1:i-1,i).*v(1:i-1)') ;
            if i ~= j
              v(i)=-s/d(i) ;
              ah(i,j)=-v(i) ;
            end
          end
          if s > -one
            d(j)=sqrt(s+one) ;
          else 
            ierr=-2 ;
            return ;
          end
        end
      end
%    --------------------------------------------------------------
%     9 Exit
%     9.1 Subcondition of constrained part of A 
      if irankc ~= 0
        v(1) = d(irankc) ;
        if v(1) ~= zero
          v(1) = abs(d(1)/v(1)) ;
        end
      else
        v(1)=zero ;
      end
%     9.1 subcondition of least squares part of a
      if k == irank
        t = d(irank) ;
      end
      if irankc+1 <= irank  &&  t ~= zero
        cond = abs(d(irankc+1)/t) ;
      else
        cond = zero ;
      end
      ierr=0 ;

function [x,b] = solcon(a,nrow,ncol,mcon,m,n,b,irankc,irank,d,pivot,kred,ah)
%     ------------------------------------------------------------
%    
%*    Summary
%     =======
%
%     Best constrained linear least squares solution of (M,N)-
%     system . First MCON rows comprise MCON equality constraints.
%     To be used in connection with function DECCON
%     References:       See DECCON
%     Related Programs: DECCON
%    
%*    Parameters:
%     ===========
%
%*    Input parameters (* marks inout parameters)
%     ===========================================
%
%       A(M,N), NROW, NCOL, M, N, MCON, IRANKC, IRANK,
%       D(N), PIVOT(N), AH(N,N), KRED
%                           See input- respective output-parameters
%                           description of function DECCON
%     * B(M)         Dble   Right-hand side of linear system, if
%                           KRED >= 0
%                           Right-hand side of upper linear system,
%                           if KRED < 0
%
%*    Output parameters
%     =================
%
%       X(N)         Dble   Best LSQ-solution of linear system
%       B(M)         Dble   Right-hand of upper trigular system
%                           (transformed right-hand side of linear
%                            system)
%
%     ------------------------------------------------------------
      zero=0.0 ;
      v = zeros(n,1) ;
%*    Begin
%     ------------------------------------------------------------
%     1 Solution for pseudo-rank zero
      if irank == 0
        x = zeros(n,1) ;
        return ;
      end
      if (irank <= irankc) && (irank ~= n)
        iranc1 = irankc + 1 ;
        v(iranc1:n) = zeros(n-iranc1+1,1) ;
      end
      if kred >= 0 && (m ~= 1 || n ~= 1)
%       ----------------------------------------------------------
%       2 Constrained householder transformations of right-hand side
        mh = mcon ;
        if irankc == 0
          mh = m ;
        end
        for j=1:irank
          s = sum(a(j:mh,j).*b(j:mh)) / (d(j)*a(j,j)) ;
          b(j:m)=b(j:m)+a(j:m,j)*s ;
          if j == irankc
            mh = m ;
          end
        end
      end
%     ------------------------------------------------------------
%     3 Solution of upper triangular system
      irk1 = irank+1 ;
      for ii=1:irank
        i = irk1-ii ;
        i1 = i+1 ;
        s = b(i) ;
        if ii ~= 1 
          s = s - sum(a(i,i1:irank).*v(i1:irank)') ;
        end
        v(i)=s/d(i) ;
      end
      if (irank ~= n) && (irank ~= irankc)
%       ----------------------------------------------------------
%       3.2 Computation of the best constrained least squares-
%           solution
        for j=irk1:n
          s = sum(ah(1:j-1,j).*v(1:j-1)) ;
          v(j)=-s/d(j) ;
        end
        for jj=1:n
          j = n-jj+1 ;
          s = zero ;
          if jj ~= 1
            s = sum(ah(j,j1:n).*v(j1:n)') ;
          end
          if jj ~= 1 && j <= irank
            v(j)=v(j)-s ;
          else
            j1 = j ;
            v(j)=-(s+v(j))/d(j) ;
          end
        end
      end
%     ------------------------------------------------------------
%     4 Back-permutation of solution components
      x = zeros(n,1) ;
      for l1=1:n
        x(pivot(l1)) = v(l1) ;
      end

%*    Group  Time monitor package
%
%*    Begin Prologue
%     ------------------------------------------------------------
%
%*  Title
%    
%     Monitor - A package for making multiple time measurements and
%               summary statistics
%
%*  Written by        U. Nowak, L. Weimann 
%*  Version           1.0
%*  Revision          July 2001
%*  Latest Change     July 2001
%*  Library           CodeLib
%*  Code              Fortran 77, Double Precision
%*  Environment       Standard Fortran 77 environment on PC's,
%*  Copyright     (c) Konrad Zuse Zentrum fuer
%                     Informationstechnik Berlin
%                     Takustrasse 7, D-14195 Berlin-Dahlem
%                     phone : + 49/30/84185-0
%                     fax   : + 49/30/84185-125
%*  Contact           Lutz Weimann 
%                     ZIB, Numerical Software Development 
%                     phone : + 49/30/84185-185
%                     fax   : + 49/30/84185-107
%                     e-mail: weimann@zib.de
%
%  ---------------------------------------------------------------
%
%* Licence
%    You may use or modify this code for your own non commercial
%    purposes for an unlimited time. 
%    In any case you should not deliver this code without a special 
%    permission of ZIB.
%    In case you intend to use the code commercially, we oblige you
%    to sign an according licence agreement with ZIB.
%
%* Warranty 
%    This code has been tested up to a certain level. Defects and
%    weaknesses, which may be included in the code, do not establish
%    any warranties by ZIB. ZIB does not take over any liabilities
%    which may follow from aquisition or application of this code.
%
%* Software status 
%    This code is under care of ZIB and belongs to ZIB software class 1.
%
%  ---------------------------------------------------------------
%
%*    Summary:
%
%     Monitor is a package for generating time and summary statistics
%     about the execution of multiple program parts of any program.
%     Nested measurements of program parts are possible.
%     ------------------------------------------------------------
%
%*    Usage:
%
%     The usage of Monitor is naturally divided into three phases:
%     1. the initialization and setup phase before the start of
%        the program or functions package to be measured;
%     2. the run phase of the program to be measured;
%     3. the final evaluation call.
%
%     The phase 1 must start with exactly one call of the function
%     MONINI, which passes a title string and a logical unit for
%     later statistics output and possible error messages to the
%     package. This call follows a number of calls of the function
%     MONDEF, where each call associates an identification string
%     to a positive integer number, called the measurement index
%     - up to maxtab, where maxtab is a package constant. Multiple
%     measurement indices may be used for measurements of multiple
%     program parts. The index 0 must also be associated with some
%     identification string, and corresponds to all parts of the
%     measured program from the measurement start call till the final
%     evaluation call, which are not associated with specific positive
%     measurement indices. After all necessary MONDEF calls are done,
%     the measurements are started at begin of the program to be
%     measured by a parameterless call of MONSRT.
%     In phase 2, each program part to be measured must be immediately
%     preceeded by a call of the function MONON with the associated 
%     measurement index, and must be immediately followed by a call of
%     the function MONOFF with the same measurement index. Measure-
%     ments of nested program parts are possible, and nesting is allowed
%     up to the number mnest, where mnest is a package constant.
%     Calling MONOFF without a preceeding MONON call with the same 
%     measurement index, or calling one of these functions with a
%     measurement index not previously defined by a MONDEF call causes
%     an error stop of the program. 
%     Finally at the end of the program to be measured, a 
%     call of the function MONEND closes all measurements and
%     prints the summary statistics.
%     As delivered, maxtab has a value 20 and mnest a value 10, but
%     both constants may be increased, if needed, to any possible
%     integer value, by simply changing it's values where they
%     are set, below.
%
%*    functions and their parameters:
%     =================================
%
%     MONINI(CIDENT,LUMON)  : Initialize Monitor
%       CIDENT  char*20  Identification string for the total measurement
%                        ( printed in summary )
%       LUMON   int      The logical unit for printing out the summary
%
%     MONDEF(MESIND,CIDMES) : Define one measurement index
%       MESIND  int      >=1 : measurement index for a specific part
%                        = 0 : measurement index for all remaining parts
%                              (i.e. not belonging to parts with 
%                               index >=1)
%       CIDMES  char*15  Identification string for the part associated
%                        with MESIND ( printed in summary )
%
%     MONSRT (dummy)        : Start measurements
%
%     MONON(MESIND)         : Start measurement of a specific part
%       MESIND  int      >=1 : measurement index for a specific part
%
%     MONOFF(MESIND)        : Stop measurement of a specific part
%       MESIND  int      >=1 : measurement index for a specific part
%
%     MONEND (DUMMY)         : Finish measurements and print summary
%       
%
%
%*    Example:
%       Calling sequence:
%
%        monini (' example',1) ;
%        mondef (0,'solver') ;
%        mondef (1,'user function') ;
%        mondef (2,'user matrix') ;
%        monsrt (dummy) ;
%        ...
%        program to be measured (part without specific measurement index)
%        ...
%        while ~ termination      
%          ...
%          monon (2) ;
%          ...  user matrix code ...
%          monoff(2) ;
%          ...
%          program to be measured (part without specific measurement index)
%          ...
%          monon (1) ;
%          ...  user function code ...
%          monoff(1) ;
%          ...
%          program to be measured (part without specific measurement index)
%          ...
%        end
%        ...
%        monend (dummy) ;
%     ------------------------------------------------------------
%
%     initialize monitor
%
function dummy = monini(texth,iounit)
%
      global maxtab mnest ;
      global sec count asec pc1 pc2 indxo time1 time0 maxind name;
      global sec0 name0 text moni qon ioncnt indact info;
      maxtab = 20 ; mnest = 10 ;
      info = 1 ;
      moni=iounit ;
      maxind=0 ;
      text=texth ;
      sec(1:maxtab)=zeros(maxtab,1) ;
      asec(1:maxtab)=zeros(maxtab,1) ;
      count(1:maxtab)=zeros(maxtab,1) ;
      qon(1:maxtab)=zeros(maxtab,1) ;
      indact(1:mnest)=zeros(mnest,1) ;
      sec0=0.0 ;
      ioncnt=0 ;
      dummy = 1 ;
%
%     define one monitor entry
%
function dummy = mondef(indx,nameh)
      global maxtab  mnest ;
      global sec count asec pc1 pc2 indxo time1 time0 maxind name;
      global sec0 name0 text moni qon ioncnt indact ;
      if indx >= 0  &&  indx <= maxtab
        if indx > maxind
          maxind=indx ;
        end
        if indx > 0
          if count(indx) > 0
            error('\n%s\n%s%4i%s\n',' error in function mondef', ...
                  '   indx = ',indx,'   already in use') ;
          end
        end
        if indx == 0
          name0 = nameh ;
        else
          name(indx,1:length(nameh)) = nameh ;
        end
      else
        error('\n%s\n%s%4i\n%s%4i\n',' error in function mondef', ...
              '   indx > ',maxtab,'   indx=',indx) ;
      end
      dummy = 1 ;
%
%     start monitor measurements
% 
function dummy = monsrt(dummy2)
      global time1 ;
      time1 = cputime ;
      dummy = 1 ;
%      if igraph > 0
%        gmini(maxind,name0,name)
%      end

%
%     start one measurement
%
function dummy = monon(indx)
      global maxtab  mnest ;
      global sec count asec pc1 pc2 indxo time1 time0 maxind name;
      global sec0 name0 text moni qon ioncnt indact info;
      if indx > maxind || indx <= 0
        error('\n%s\n%s%4i\n',' error in function monon', ...
                      '   indx out of range    indx=',indx) ;
      end
      if qon(indx)
        error('\n%s\n%s%4i\n',' error in function monon', ...
              '   measurement is already running for    indx=',indx) ;
      end
      asec(indx) = cputime ;
      qon(indx)=1 ;
      if ioncnt == 0
        sec0=sec0+asec(indx)-time1 ;
      else
        indxo=indact(ioncnt) ;
        sec(indxo)=sec(indxo)+asec(indx)-asec(indxo) ;
      end
      ioncnt=ioncnt+1 ;
      indact(ioncnt)=indx ;
      if info > 1
        fprintf(moni,'%7s%15s%7.5f',' enter ',name(indx,:),asec(indx)) ;
      end
      dummy = 1 ;
%
%      if igraph > 0
%        call gmon(indx,sec0)
%      end

%
%     stop one measurement
%
function  dummy = monoff(indx)
      global maxtab  mnest ;
      global sec count asec pc1 pc2 indxo time1 time0 maxind name;
      global sec0 name0 text moni qon ioncnt indact info;
      if indx > maxind || indx <= 0
        error('\n%s\n%s%4i\n',' error in function monoff', ...
                      '   indx out of range    indx=',indx) ;
      end
      if ~ qon(indx)
        error('\n%s\n%s%4i\n',' error in function monoff', ...
              '   measurement has never been activated for     indx=',indx) ;
      end
      time2 = cputime ;
      qon(indx)=0 ;
      sec(indx)=sec(indx)+time2-asec(indx) ;
      count(indx)=count(indx)+1 ;
      ioncnt=ioncnt-1 ;
      if ioncnt == 0
        time1=time2 ;
      else
        asec(indact(ioncnt))=time2 ;
      end
      if info > 1
        fprintf(moni,'%6s%15s%7.5f',' exit ',name(indx,:),asec(indx)) ;
      end
      dummy = 1 ;
%
%      if igraph > 0
%        call gmoff(indx,sec(indx))
%      end
%

%
%     terminate monitor and print statistics
%
function dummy = monend(dummy2)
      global maxtab  mnest ;
      global sec count asec pc1 pc2 indxo time1 time0 maxind name;
      global sec0 name0 text moni qon ioncnt indact ;
      time0 = cputime ;
      sec0=sec0+time0-time1 ;
%
      sum1=1.0e-10 ;
      sum1 = sum1 + sum( sec(1:maxind)) ;
      for i=1:maxind
        if count(i) > 0
          asec(i)=sec(i)/count(i) ;
        end
      end
      sum0=sum1+sec0 ;
%
      pc1(1:maxind)=100.0*sec(1:maxind)/sum0 ;
      pc2(1:maxind)=100.0*sec(1:maxind)/sum1 ;
      pc10=100.0*sec0/sum0 ;
      pc20=100.0*sec0/sum1 ;
%
      fprintf(moni,'\n\n%1s\n',' ') ;
      fprintf(moni,'%s\n', ...
              ' ###########################################################################') ;
      fprintf(moni,'%s\n%s\n', ...
           ' #                                                                         #', ...
           ' #                                                                         #') ;
      fprintf(moni,'%s%-29s%3s\n',' #   Results from time monitor program for: ', text,'  #') ;
      fprintf(moni,'%s\n', ...
              ' #                                                                         #' ) ;
      fprintf(moni,'%16s%11.3f     %13s%11.3f                  %1s\n', ...
                   ' #   Total time:',sum0,'Sum of parts:',sum1,' #') ;
      fprintf(moni,'%s\n', ...
              ' #                                                                         #' ) ;
      fprintf(moni, ...
              ' #     name            calls       time    av-time    %% total      %% sum   #\n') ;
%
      i0=1 ;
      fprintf(moni,' #   %-15s%8i%11.3f%11.4f%11.2f%11.2f   #\n',name0,i0,sec0,sec0,pc10,pc20) ;
%
      for i=1:maxind
        fprintf(moni,' #   %-15s%8i%11.3f%11.4f%11.2f%11.2f   #\n', ...
                     char(name(i,:)),count(i),sec(i),asec(i),pc1(i),pc2(i) ) ;
      end
%
      fprintf(moni,'%s\n', ...
              ' #                                                                         #' ) ;
      fprintf(moni,'%s\n', ...
              ' ###########################################################################') ;
      fprintf(moni,'\n\n%1s\n',' ') ;
      dummy = 1 ;
%      if igraph > 0
%         gmend ;
%      end
%
%  end package monitor

%
%*    Group  Machine dependent functions

function dmach = d1mach(i)
%
%  DOUBLE-PRECISION MACHINE CONSTANTS
%
%  D1MACH( 1) = B**(EMIN-1), THE SMALLEST POSITIVE MAGNITUDE.
%
%  D1MACH( 2) = B**EMAX*(1 - B**(-T)), THE LARGEST MAGNITUDE.
%
%  D1MACH( 3) = B**(-T), THE SMALLEST RELATIVE SPACING.
%
%  D1MACH( 4) = B**(1-T), THE LARGEST RELATIVE SPACING.
%
%  D1MACH( 5) = LOG10(B)
%
%  TO ALTER THIS FUNCTION FOR A PARTICULAR ENVIRONMENT,
%  THE DESIRED SET OF DATA STATEMENTS SHOULD BE ACTIVATED BY
%  REMOVING THE C FROM COLUMN 1.
%  ON RARE MACHINES A STATIC STATEMENT MAY NEED TO BE ADDED.
%  (BUT PROBABLY MORE SYSTEMS PROHIBIT IT THAN REQUIRE IT.)
%
%  FOR IEEE-ARITHMETIC MACHINES (BINARY STANDARD), ONE OF THE FIRST
%  TWO SETS OF CONSTANTS BELOW SHOULD BE APPROPRIATE.
%
%  WHERE POSSIBLE, DECIMAL, OCTAL OR HEXADECIMAL CONSTANTS ARE USED
%  TO SPECIFY THE CONSTANTS EXACTLY.  SOMETIMES THIS REQUIRES USING
%  EQUIVALENT INTEGER ARRAYS.  IF YOUR COMPILER USES HALF-WORD
%  INTEGERS BY DEFAULT (SOMETIMES CALLED INTEGER*2), YOU MAY NEED TO
%  CHANGE INTEGER TO INTEGER*4 OR OTHERWISE INSTRUCT YOUR COMPILER
%  TO USE FULL-WORD INTEGERS IN THE NEXT 5 DECLARATIONS.
%
%     dmach(1) = realmin ;
%     dmach(2) = realmax ;
%     dmach(3) = eps ;
%     dmach(4) = 2.0*eps ;
%     dmach(5) = log10(2.0) ;
      if i  <  1   ||   i  >  6
        error(' d1mach - i out of bounds%10i',i) ;
      end
      if i  ==  1 
        dmach = realmin ;
      elseif i  ==  2
        dmach = realmax ;
      elseif i  ==  3
        dmach = eps ;
      elseif i  ==  4
        dmach = 2.0*eps ;
      elseif i  ==  5
        dmach = log10(2.0) ;
      elseif i  ==  6
%       dmach = sqrt(dmach(1)/dmach(3)) ;
        dmach = 4.94e-32 ;
      end

function [optval,opt] = getopt(opt,field,default)
      if isfield(opt,field)
        optval = getfield(opt,field) ;
      else
        optval = default ;
        opt = setfield(opt,field,default) ;
      end
      
function opt = iniopt(opt,field,default)
      if ~ isfield(opt,field)
        opt = setfield(opt,field,default) ;
      end
