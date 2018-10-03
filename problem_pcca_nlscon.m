function varargout = problem_pcca_nlscon(X,flag,par)
%
% Problem function.
%
% varargout = problem_pcca_nlscon(X,flag,par) 
%
% computes the problem function f(X), Jacobian df/dx(X),
% startvector for the Gauss-Newton iteration 
% and vector of to be fitted values.
%
% FX = problem_pcca_nlscon(X,'',par)          
%    returns the right hand side f(X) for a input column vector X.
% JAC = problem_pcca_nlscon(X,'jacobian',t) 
%    returns the Jacobian df/dx(X) for a input column vector X.
%
% Input:
%       X()         Float       Vector of unknowns (input)
%       flag        String      Operation flag - the following values must 
%                               be supported:
%                               '' (empty string) : FOUT must return the 
%                                   problem-function value f(X) ;
%                               'jacobian' : FOUT must return the 
%                                   associated Jacobian Jac(x)
%       par         AnyType     A (required!) user parameter
%
% Output:
%       FOUT()      Float       Vector of returned function values or 
%                               Jacobian matrix
%       FAIL        Int         evaluation-failure indicator. (output)
%                               On output: Indicates failure of FCN eval-
%                               uation, if having a value <= 2.
%                               If <0: NLEQ1 will be terminated with 
%                                    error code = 82, and FAIL stored
%                                    to wk.ifail.
%                               If =1: A new trial Newton iterate will
%                                    computed, with the damping factor
%                                    reduced to its half.
%                               If =2: A new trial Newton iterate will
%                                    computed, with the damping factor
%                                    reduced by a reduct. factor, which
%                                    must be output through F(1) by FCN,
%                                    and it is value must be >0 and < 1.
%                               Note, that if FAIL = 1 or 2, additional
%                               conditions concerning the damping factor,
%                               e.g. the minimum damping factor or the
%                               bounded damping strategy may also influ-
%                               ence the value of the reduced damping 
%                               factor.  
% Modified by Bernhard Reuter, Theoretical Physics II,
% University of Kassel, 2017

    switch flag
        case '' % Return y = f(X).
            [varargout{1:2}] = f_rhs(X,par) ;
        case 'jacobian' % Return Jacobian matrix df/dx.
            [varargout{1:2}] = jacobian(X,par);
        otherwise
            error(['Unknown flag ''' flag '''.']) ;
    end
end
% --------------------------------------------------------------------------
 

function [f,ifail] = f_rhs(x,par)

    ifail = -1 ;

    EVS = par.evs ;
    k = size(EVS,2) ;

    A = zeros(k,k) ;

    for i = 1:k-1
        A(i+1,2:k) = x(((i-1)*(k-1)+1):i*(k-1)) ;
    end

    A = fillA(A, EVS) ;

    F = trace(diag(1./A(1,:))*(A'*A)) ;

    f = k-F ;

    ifail = 0 ;

end
% --------------------------------------------------------------------------

function [JF,ifail] = jacobian(x2,par)

    EVS = par.evs ;
    k = size(EVS,2) ;

  
    ifail = -1 ;

    A = zeros(k,k) ;


    for i = 1:k-1
        A(i+1,2:k) = x2(((i-1)*(k-1)+1):i*(k-1)) ;
    end

    A = fillA(A,EVS) ;

    JF = zeros(k,k) ;
    for i = 1:k
        for j = 1:k
            JF(i,j) = A(i,j)/A(1,j) ;
        end
    end

    JF(1,:) = [] ;
    JF(:,1) = [] ;

    JF = JF(:) ;
    JF = JF' ;

    ifail = 0 ;
end

