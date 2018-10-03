function [ Q, R, ap ] = SRSchur_num_t( Q, R, z, b )
    % SYNTAX: [ Q, R, ap ] = SRSchur( Q, R, z, b )
    %
    % INPUT: orthogonal real Q and quasi-triangular real R such that AQ=QR 
    % and a target z in the complex plane. The fourth parameter b 
    % determines the length of the ordering with respect to z to be produced:
    %
    % if b < 0 then -b blocks will be sorted,
    % if b > 0 then  b or b+1 eigenvalues will be sorted, depending on the 
    % sizes of the blocks,
    % if b = 0 then the whole Schur form will be sorted.
    %
    % OUTPUT: orthogonal real Q and quasi-triangular real R such that AQ=QR 
    % with the diagonal blocks ordered with respect to the target z. The 
    % number ofordered blocks/eigenvalues is determined by the parameter b.
    % A vector ap warns for inaccuracy of the solution if an entry of ap 
    % exceeds one.
    %
    % SUBFUNCTIONS: normalize.m, swaplist.m, select.m, swap.m, lu_complpiv.m
    % REFERENCES: 
    % http://m2matlabdb.ma.tum.de/download.jsp?MC_ID=3&MP_ID=119
    % http://citeseer.uark.edu:8080/citeseerx/viewdoc/summary?doi=10.1.1.24.8870
    % Brandts, J. H. Matlab Code for Sorting Real Schur Forms. Numer. 
    % Linear Algebr. with Appl. 2002, 9 (3), 249-261
    % SEE ALSO: schur.m, rsf2csf.m
    %
    % Modified by Bernhard Reuter (B.R.), Theoretical Physics II, 
    % University of Kassel, 2017
    %
    % multiprecision is now supported by num_t (B.R.)
    %
    % CAUTION: At the moment z is ignored and the eigenvalues are just 
    % sorted such that the highest eigenvalue is on the top left of R and 
    % the other eigenvalues are sorted in descending order (b still 
    % determines the length of the ordering). If you want to change this, 
    % change the select() function!
%   -----------------------------------------------------------------------
    class_t1 = class(Q) ;
    class_t2 = class(R) ;
    if (strcmpi(class_t1,class_t2))
        class_t = class_t1 ;
    else
        error( 'SRSchur_num_t:DataTypeError', ...
            ['class(Q) is not equal class(R)! This will lead to ' ...
            'numeric precision errors!'] )
    end

    function rr = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), rr = mp(expression) ;
            else
                if isnumeric(expression)
                    rr = expression ;
                else
                    rr = eval(expression) ;
                end
            end
        else
            rr = class_t;
        end
    end  % num_t
%   -----------------------------------------------------------------------
%   assertions
    assert( size(Q,1)==size(Q,2), 'SRSchur_num_t:MatrixShapeError1', ...
        'Q matrix isnt quadratic!' )
    assert( size(R,1)==size(R,2), 'SRSchur_num_t:MatrixShapeError2', ...
        'R matrix isnt quadratic!' )
    assert( size(Q,1)==size(R,1), 'SRSchur_num_t:MatchError', ['The ' ...
        'dimensions of R doesnt match to those of Q!'] )
    assert( size(Q,1)>=2, 'SRSchur_num_t:MatrixShapeError3', ['Q and R ' ...
        'must be at least 2x2!'] )
%   -----------------------------------------------------------------------

    r = find(abs(diag(R,-1)) > num_t('100')*num_t('eps'));
    s = 1:size(R,1)+1;
    s(r+1) = [];
    ap = []; %ap was not defined in cases were swaplist was empty! B.R. 7.07.17

    for k=1:length(s)-1;
        sk = s(k);
        if s(k+1)-sk == 2
            [Q,R] = normalize(Q,R,sk:s(k+1)-1);
            p(k)  = R(sk,sk)+sqrt(R(sk+1,sk)*R(sk,sk+1));
        else
            p(k)  = R(s(k),s(k));
        end
    end

    for k = swaplist(p,s,z,b);
        v      = s(k):s(k+1)-1;
        w      = s(k+1):s(k+2)-1;
        nrA    = norm(R([v,w],[v,w]),inf);
        [Q,R]  = swap(Q,R,v,w);
        s(k+1) = s(k)+s(k+2)-s(k+1);
        v      = s(k):s(k+1)-1;
        w      = s(k+1):s(k+2)-1;
        if length(v)==2
            [Q,R] = normalize(Q,R,v);
        end
        if length(w)==2
            [Q,R] = normalize(Q,R,w);
        end
        ap(k)  = norm(R(w,v),inf)/(num_t('10')*num_t('eps')*nrA);
    end


    R = R - tril(R,-2);
    for j=2:length(s)-1; R(s(j),s(j)-1)=num_t('0'); end

end


% ----------------------------------------------%

function [U,S] = normalize(U,S,v)
    n  = size(S,1);
    Q  = rot(S(v,v));
    S(:,v) = S(:,v)*Q;
    S(v,:) = Q'*S(v,:);
    U(:,v) = U(:,v)*Q;
end

% ----------------------------------------------%

function Q = rot(X)

    class_t = class(X);
    function rr = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), rr = mp(expression);
            else
                if isnumeric(expression)
                    rr = expression;
                else
                    rr = eval(expression);
                end
            end
        else
            rr = class_t;
        end
    end  % num_t

    c = num_t('1'); s = num_t('0');
    if X(1,1)~=X(2,2);
        tau   = (X(1,2)+X(2,1))/(X(1,1)-X(2,2));
        off   = sqrt(tau^2+num_t('1'));
        v     = [tau - off, tau + off];
        [~,w] = min(abs(v));
        c     = num_t('1')/sqrt(num_t('1')+v(w)^2);
        s     = v(w)*c;
    end
    Q = [c -s;s c];
end

% ----------------------------------------------%

function v = swaplist(p,s,z,b)
    n = length(p);
    k = 0; v = [];
    srtd = 0;
    q = diff(s);
    fini = 0;
    while ~fini
        k        = k+1;
        [~,j]  = select(p(k:n),z);
        p(k:n+1) = [p(j+k-1) p(k:n)];
        p(j+k)   = [];
        q(k:n+1) = [q(j+k-1) q(k:n)];
        q(j+k)   = [];
        v        = [v,j+k-2:-1:k];
        srtd     = srtd + q(k);
        fini     = (k==n-1)|(k==-b)|(srtd==b)|((srtd==b+1)&(b~=0));
    end
end

% ----------------------------------------------%

function [val,pos] = select(p,z)
    %y = real(z)+abs(imag(z))*i;
    %[val pos] = min(abs(p-y)); %Marcus: RAUS... stattdessen...
    [val, pos]=max(abs(p)); %fuer Permutationsmatrizen
    %[val, pos]=min(abs(p-y) + 1000* abs(imag(p))); %fuer reelle Schurwerte
    %[val, pos]=max(abs(p-(real(p)<0).*real(p))); %fuer logarythmierbare Matrizen
end

% -----------------------------------------------%

function [U,S] = swap(U,S,v,w)

    class_t1 = class(U);
    class_t2 = class(S);
    if (strcmpi(class_t1,class_t2))
        class_t = class_t1;
    else
        error('swap:U_S_DataTypeError', ...
            ['class(U) is not equal class(S)! This will lead to numeric '...
            'precision errors!'])
    end

    function rr = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), rr = mp(expression);
            else
                if isnumeric(expression)
                    rr = expression;
                else
                    rr = eval(expression);
                end;
            end;
        else
            rr = class_t;
        end;
    end  % num_t

    [p,q] = size(S(v,w)); Ip = eye(p,num_t); Iq = eye(q,num_t);
    r = [];
    for j=1:q
        r = [r;S(v,w(j))];
    end
    K = kron(Iq,S(v,v))-kron(S(w,w)',Ip);
    [L,H,P,Q] = lu_complpiv(K);
    e = min(abs(diag(H)));
    sigp = 1:p*q;
    for k = 1:p*q-1;
        sigp([k,P(k)]) = sigp([P(k),k]);
    end
    r = e*r(sigp);
    x = (H\(L\r));
    sigq = 1:p*q;
    for k = 1:p*q-1;
        sigq([k,Q(k)]) = sigq([Q(k),k]);
    end
    x(sigq) = x;
    X = [];
    for j=1:q
        X = [X,x((j-1)*p+1:j*p)];
    end
    [Q,R]      = qr([-X;e*Iq]);
    S(:,[v,w]) = S(:,[v,w])*Q;
    S([v,w],:) = Q'*S([v,w],:);
    U(:,[v,w]) = U(:,[v,w])*Q;
end

% -----------------------------------------------%

function [L,U,P,Q] = lu_complpiv(A)

    class_t=class(A);
    function rr = num_t(expression)
        if (nargin > 0)
            if(strcmpi(class_t,'mp')), rr = mp(expression);
            else
                if isnumeric(expression)
                    rr = expression;
                else
                    rr = eval(expression);
                end
            end
        else
            rr = class_t;
        end
    end  % num_t

    P = []; Q = []; n = size(A,1);
    for k=1:n-1;
        [a,r] = max(abs(A(k:n,k:n)));
        [~,c] = max(abs(a));
        cl  = c+k-1;
        rw  = r(c)+k-1;
        A([k,rw],:) = A([rw,k],:);
        A(:,[k,cl]) = A(:,[cl,k]);
        P(k) = rw; Q(k) = cl;
        if A(k,k) ~= num_t('0');
            rs = k+1:n;
            A(rs,k)  = A(rs,k)/A(k,k);
            A(rs,rs) = A(rs,rs)-A(rs,k)*A(k,rs);
        end
    end
    U = tril(A')'; L = tril(A,-1) + eye(n,num_t);
end