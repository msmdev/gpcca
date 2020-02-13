function [ X, RR ]  = do_schur( P, sd, fileid, b )
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
% Perform Schur decomposition and Gram-Schmidt orthogonalization.
    % 
    % [ X, RR ]  = do_schur( P, sd, fileid )
    %
    % Input:
    %   P           stochastic (fine) transition (N,N)-matrix
    %   sd          "initial" distribution
    %   fileid      ID-string for the naming of output files
    %   b           if b < 0 then -b blocks will be sorted,
    %               if b > 0 then  b or b+1 eigenvalues will be sorted, 
    %               depending on the sizes of the blocks,
    %               if b = 0 then the whole Schur form will be sorted.
    % Output:
    %   X           (N,n)-matrix containing the ordered Schur-vectors 
    %               columnwise
    %   RR          (N,N) ordered Schur matrix 
    %
    % Written by Bernhard Reuter, Theoretical Physics II,
    % University of Kassel, 2017
%   -----------------------------------------------------------------------
    class_t1 = class(P) ;
    class_t2 = class(sd) ;
    if (strcmp(class_t1,class_t2))
        class_t = class_t1 ;
    else
        error('do_schur:DataTypeError', ...
            ['class(P) is not equal class(sd)! This will lead to ', ...
            'numeric precision errors!'])
    end

    function r = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), r = mp(expression) ;
            else
                if isnumeric(expression)
                    r = expression ;
                else
                    r = eval(expression) ;
                end ;
            end ;
        else
            r = class_t ;
        end ;
    end  % num_t
%   -----------------------------------------------------------------------
%       assertions
    N = size(P,1) ;
    assert( size(P,1)==size(P,2), 'do_schur:MatrixShapeError1', ...
        'P matrix isnt quadratic!' )
    assert( size(sd,1)==size(P,1), 'do_schur:MatrixShapeError2', ...
        'sd vector length doesnt match with the shape of P!' )
    assert(all(sum(P,2) > num_t('0.0')),'do_schur:ZeroError1', ...
        'Not all rows of P are >0!')
    assert(all(sd > num_t('0.0')),'do_schur:ZeroError2', ...
        'Not all elements of sd are >0!')
    
%       weight the stochastic matrix P_stoch by the initial distribution sd
    Pd=diag(sqrt(sd))*P*diag(num_t('1.0')./sqrt(sd)) ;
    assert(isa(Pd,num_t),'do_schur:Pd_DataTypeError', ...
        'Variable is type %s not %s', class(Pd),num_t)
    name=strcat(fileid,'-','Pd','.txt') ;
    save_t(name,Pd,'-ascii')

%       make a Schur decomposition of Pd
    [Q, R]=schur(Pd) ;
    assert(isa(Q,num_t),'do_schur:Q_DataTypeError', ...
        'Variable is type %s not %s',class(Q),num_t)
    assert(isa(R,num_t),'do_schur:R_DataTypeError', ...
        'Variable is type %s not %s',class(R),num_t)
    name=strcat(fileid,'-','Q','.txt') ;
    save_t(name,Q,'-ascii')
    name=strcat(fileid,'-','R','.txt') ;
    save_t(name,R,'-ascii')

%       sort the Schur matrix and vectors
    target=num_t('1.0') ;
    [ QQ, RR, ap ] = SRSchur_num_t( Q, R, target, b ) ;
    assert(isa(QQ,num_t),'do_schur:QQ_DataTypeError', ...
        'Variable is type %s not %s',class(QQ),num_t)
    assert(isa(RR,num_t),'do_schur:RR_DataTypeError', ...
        'Variable is type %s not %s',class(RR),num_t)
    assert((any(ap > num_t('1.0')) == 0),'do_schur:ap_ValueError', ...
        'ap from SRSchur exceeds one')
    name=strcat(fileid,'-','QQ','.txt') ;
    save_t(name,QQ,'-ascii')
    name=strcat(fileid,'-','RR','.txt') ;
    save_t(name,RR,'-ascii')

%       orthonormalize the sorted Schur vectors QQ
%       by modified Gram-Schmidt-Orthonormalization
    QQ_orthonorm=gram_schmidt_mod(QQ,sd) ;
    assert(isa(QQ_orthonorm,num_t),'do_schur:QQ_orthonorm_DataTypeError', ...
        'Variable is type %s not %s', class(QQ_orthonorm),num_t)
    dummy = ((QQ_orthonorm'*QQ_orthonorm - num_t(eye(size(QQ_orthonorm,2)))) ...
        > (num_t('10000')*eps(num_t))) ;
    assert(~any(dummy(:)),'do_schur:QQ_orthonorm_OrthogonalityError', ...
        'Schurvectors appear to not be orthogonal')
    clearvars dummy
    name=strcat(fileid,'-','QQ_orthonorm','.txt') ;
    save_t(name,QQ_orthonorm,'-ascii')

%       transform the orthonormalized Schur vectors of Pd back
%       -> orthonormalized Schur vectors X of P_stoch
    X=diag(num_t('1.0')./sqrt(sd))*QQ_orthonorm ;
    assert( size(X,1)==N, 'do_schur:MatrixShapeError4', ...
        'The shape (%d,%d) of X doesnt match with the shape (%d,%d) of P!', ...
        size(X,1), size(X,2), N, N)
    assert(isa(X,num_t),'do_schur:X_DataTypeError', ...
        'X is type %s not %s',class(X),num_t)
    dummy = (X'*diag(sd)*X - num_t(eye(size(X,2))) > (num_t('10000')*eps(num_t))) ;
    assert(~any(dummy(:)),'do_schur:X_OrthogonalityError', ...
        'Schurvectors appear to not be orthogonal')
    clearvars dummy
    dummy = ( abs(X(:,1) - 1) < ( num_t('1000') * eps(num_t) ) ) ;
    assert(all(dummy(:)), 'do_schur:FirstColumnError', ...
        'X(:,1) isnt equal 1!')
    name=strcat(fileid,'-','X','.txt') ;
    save_t(name,X,'-ascii')

end
