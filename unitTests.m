function tests = unitTests
% This file is part of GPCCA.
%
% Copyright (c) 2018, 2017 Bernhard Reuter
%
% If you use this code or parts of it, cite the following reference:
%
% Reuter, B., Weber, M., Fackeldey, K., R?blitz, S., & Garcia, M. E. (2018). Generalized
% Markov State Modeling Method for Nonequilibrium Biomolecular Dynamics: Exemplified on
% Amyloid beta Conformational Dynamics Driven by an Oscillating Electric Field. Journal of
% Chemical Theory and Computation, 14(7), 3579-3594. https://doi.org/10.1021/acs.jctc.8b00079
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
%   All test functions written by Bernhard Reuter, Theoretical Physics II,
%   University of Kassel, 2017

    tests = functiontests(localfunctions) ;
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% test do_schur (depends on gram_schmidt_mod, SRSchur_num_t, numeric_t) 
% (logic paths: 3) (assertions: 2+12) 
% The 12 assertions here can't be tested (its also not necessary, since they
% are tested elsewhere)

function do_Schur_MatrixShapeError1(testCase) % a1
%       test for the assertion if P matrix isnt quadratic
    d1 = numeric_t('[0,0,0]') ;
    d2 = numeric_t('[0,0,0]') ;
    d3 = numeric_t('[0,0,0]') ;
    d4 = numeric_t('[0,0,0]') ;
    P = [ d1; d2; d3; d4 ] ;
    sd = numeric_t('[ 0, 0, 0, 0 ]') ;
    fileid = 'test_do_schur_MatrixShapeError1' ;
    try
        [ ~, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'do_schur:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_do_schur_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function do_Schur_MatrixShapeError2(testCase) % a1
%       test for the assertion if sd vector length doesnt match with the 
%       shape of P
    d1 = numeric_t('[0,0,0]') ;
    d2 = numeric_t('[0,0,0]') ;
    d3 = numeric_t('[0,0,0]') ;
    P = [ d1; d2; d3 ] ;
    sd = numeric_t('[ 0, 0, 0, 0 ]') ;
    fileid = 'test_do_schur_MatrixShapeError2' ;
    try
        [ ~, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'do_schur:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_do_schur_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_do_schur(testCase) % l1
%       test normal execution
%       specify the test case data
    matrixfile = 'example_matrix_mu0.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_do_schur' ;
    [ X, RR ]  = do_schur( P, sd, fileid, 0 ) ;
    N = 9 ; 
    verifyTrue(testCase,all(size(P)==[N,N]))
    verifyTrue(testCase,all(size(X)==[N,N]))
    verifyTrue(testCase,all(size(RR)==[N,N]))
%       test, if P*X=X*RR (Schur decomposition)
    dummy = ( abs(X*RR - P*X) < (testCase.TestData.abstolfac * eps(numeric_t)) ) ;
    verifyTrue(testCase,all(dummy(:)))
%       test, if the first column of X is 1
    dummy = (abs(X(:,1) - 1) < (testCase.TestData.abstolfac * eps(numeric_t))) ;
    verifyTrue(testCase,all(dummy(:)))
end

function test_do_schur_DataTypeError(testCase) % l2
%       test the error for the case, that class(P)~=class(sd)
%       specify the test case data
    matrixfile = 'example_matrix_mu0.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_do_schur_DataTypeError' ;
    if isa(P,'mp') && isa(sd,'mp')
        sd = double(sd) ;
    elseif isa(P,'double') && isa(sd,'double')
        sd = single(sd) ;
    else
        error('test_do_schur_DataTypeError:DataTypeError', ...
        ['class(P)=',class(P),' class(sd)=',class(sd)])
    end
    try
        [ ~, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'do_schur:DataTypeError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_do_schur_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_do_schur_b_pos(testCase) % l3a
%       test normal execution with b>0
%       specify the test case data
    matrixfile = 'example_matrix_mu0.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_do_schur' ;
    [ X, RR ]  = do_schur( P, sd, fileid, 3 ) ;
    N = 9 ; 
    verifyTrue(testCase,all(size(P)==[N,N]))
    verifyTrue(testCase,all(size(X)==[N,3]))
    verifyTrue(testCase,all(size(RR)==[N,N]))
%       test, if P*X=X*RR (Schur decomposition)
    dummy = ( abs(X*RR(1:3,1:3) - P*X) < (testCase.TestData.abstolfac * eps(numeric_t)) ) ;
    verifyTrue(testCase,all(dummy(:)))
%       test, if the first column of X is 1
    dummy = (abs(X(:,1) - 1) < (testCase.TestData.abstolfac * eps(numeric_t))) ;
    verifyTrue(testCase,all(dummy(:)))
end

function test_do_schur_b_neg(testCase) % l3b
%       test normal execution with b<0
%       specify the test case data
    matrixfile = 'example_matrix_mu0.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_do_schur' ;
    [ X, RR ]  = do_schur( P, sd, fileid, -3 ) ;
    N = 9 ; 
    verifyTrue(testCase,all(size(P)==[N,N]))
    verifyTrue(testCase,all(size(X)==[N,3]))
    verifyTrue(testCase,all(size(RR)==[N,N]))
%       test, if P*X=X*RR (Schur decomposition)
    dummy = ( abs(X*RR(1:3,1:3) - P*X) < (testCase.TestData.abstolfac * eps(numeric_t)) ) ;
    verifyTrue(testCase,all(dummy(:)))
%       test, if the first column of X is 1
    dummy = (abs(X(:,1) - 1) < (testCase.TestData.abstolfac * eps(numeric_t))) ;
    verifyTrue(testCase,all(dummy(:)))
end
% -------------------------------------------------------------------------
% test fillA (depends on do_schur, initialize_A, numeric_t) (logic paths: 6)
% (assertions: 6)

function test_fillA_MatrixShapeError1(testCase) % a1
%       test assertion, if A matrix isnt quadratic
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    A = svecs ;
    try
        [ ~ ] = fillA( A, svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'fillA:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_fillA_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_fillA_MatchError(testCase) % a2
%       test assertion, if second dimension of schur vector matrix doesnt
%       match to the dimensions of the A matrix
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    A = zeros(3) ;
    try
        [ ~ ] = fillA( A, svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'fillA:MatchError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_fillA_MatrixShapeError2(testCase) % a3
%       test assertion for A smaller 2x2
    svecs = [ 1, 1 ]' ;
    A = 0 ;
    try
        [ ~ ] = fillA( A, svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'fillA:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_fillA_MatrixShapeError3(testCase) % a4
%       test assertion for (N,k)-matrix with k>N
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    A = zeros(4) ;
    try
        [ ~ ] = fillA( A, svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'fillA:MatrixShapeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end


function test_fillA_MatrixShapeError4(testCase) % a5
%       test assertion for (N,k)-matrix with k=N
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    A = zeros(3) ;
    try
        [ ~ ] = fillA( A, svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'fillA:MatrixShapeError4'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_fillA_FirstColumnError(testCase) % a6
%       test for the assertion if first column of schur vector matrix isnt
%       constantly equal 1
    d1 = [0,0,0] ;
    d2 = [0,0,0] ;
    d3 = [0,0,0] ;
    d4 = [0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    A = zeros(3) ; 
    try
        [ ~ ] = fillA( A, svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'fillA:FirstColumnError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_fillA_mu0(testCase) % l
%       specify the test case data
    matrixfile = 'example_matrix_mu0.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_fillA_mu0' ;
%       get Schur vectors
    [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    svecs = X(:,1:3) ;
%       initialize A
    [ T, A ] = evalc('initialize_A( svecs, 1 )') ;
%       just save
    name = 'test_fillA_mu0-initial-A.txt' ;
    save_t(name, A, '-ascii')
%       fill A
    [ T1, A ] = evalc('fillA( A, double(svecs))') ;
%       T, T1 is stored to a file
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s\n','test_fillA_mu0') ;
    fprintf(fileID,'%s\n',T) ;
    fprintf(fileID,'%s\n\n',T1) ;
    fclose(fileID) ; 
%       just save
    name = 'test_fillA_mu0-A.txt' ;
    save_t(name, A, '-ascii')
    name = 'test_fillA_mu0-svecs(N,1-3).txt' ;
    save_t(name, double(svecs), '-ascii')
%       actaul testing
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       check positivity
    [ m, k ] = size(A) ;
    [ N, n ] = size(svecs) ;
    verifyTrue( testCase, m == k )
    verifyTrue( testCase, n == k )
    A_exp = inf(k) ;
    for j = 1:k
        dummy = [] ;
        for l = 1:N
            dummy(l) = svecs(l,2:k)*A(2:k,j) ;
        end 
        A_exp(1,j) = -min(dummy) ;
    end
    testCase.verifyThat(A(1,:), IsEqualTo(A_exp(1,:), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       check partition of unity
    A_exp(1,1) = 1 - sum(A(1,2:k)) ;
    for i = 2:k
        A_exp(i,1) = - sum(A(i,2:k)) ;
    end
    testCase.verifyThat(A(:,1), IsEqualTo(A_exp(:,1), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

function test_fillA_mu1000(testCase) % l
%       specify the test case data
    matrixfile = 'example_matrix_mu1000.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_fillA_mu1000' ;
%       get Schur vectors
    [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    svecs = X(:,1:5) ;
%       initialize A
    [ T, A ] = evalc('initialize_A( svecs, 1 )') ;
%       just save
    name = 'test_fillA_mu1000-initial-A.txt' ;
    save_t(name, A, '-ascii')
%       fill A
    [ T1, A ] = evalc('fillA( A, double(svecs))') ;
%       T, T1 is stored to a file
    fileID = fopen('test_fillA_report.txt','a') ;
    fprintf(fileID,'%s\n','test_fillA_mu1000') ;
    fprintf(fileID,'%s\n',T) ;
    fprintf(fileID,'%s\n\n',T1) ;
    fclose(fileID) ; 
%       just save
    name = 'test_fillA_mu1000-A.txt' ;
    save_t(name, A, '-ascii')
    name = 'test_fillA_mu1000-svecs(N,1-5).txt' ;
    save_t(name, double(svecs), '-ascii')
%       actaul testing
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       check positivity
    [ m, k ] = size(A) ;
    [ N, n ] = size(svecs) ;
    verifyTrue( testCase, m == k )
    verifyTrue( testCase, n == k )
    A_exp = inf(k) ;
    for j = 1:k
        dummy = [] ;
        for l = 1:N
            dummy(l) = svecs(l,2:k)*A(2:k,j) ;
        end 
        A_exp(1,j) = -min(dummy) ;
    end
    testCase.verifyThat(A(1,:), IsEqualTo(A_exp(1,:), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       check partition of unity
    A_exp(1,1) = 1 - sum(A(1,2:k)) ;
    for i = 2:k
        A_exp(i,1) = - sum(A(i,2:k)) ;
    end
    testCase.verifyThat(A(:,1), IsEqualTo(A_exp(:,1), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

% -------------------------------------------------------------------------
% test getopt (depends on none) (logic paths: 2) (assertions: 0)

function test_getopt(testCase) % l1 + l2
    dummy.field = 111 ;
    [ field, dummy ] = getopt( dummy, 'field', 999 ) ;
    verifyEqual(testCase,field,111)
    verifyTrue(testCase,isstruct(dummy))
    verifyEqual(testCase,dummy.field,111)
    clearvars dummy
    dummy = [] ;
    [ field, dummy ] = getopt( dummy, 'field', 999 ) ;
    verifyEqual(testCase,field,999)
    verifyTrue(testCase,isstruct(dummy))
    verifyEqual(testCase,dummy.field,999)
end
% -------------------------------------------------------------------------
% test gram_schmidt_mod (depends on nearly_equal, numeric_t) (logic paths: 3)
% (assertions: 2)

function test_gram_schmidt_mod_MatrixShapeError1(testCase) % a1
%       test assertion for matrix with only one row
    sys = numeric_t('[ 3, 1 ]') ;
    sd = numeric_t('1') ;
    try
        [ ~ ] = gram_schmidt_mod( sys, sd ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'gram_schmidt_mod:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_gram_schmidt_mod_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_gram_schmidt_mod_MatrixShapeError2(testCase) % a2
%       test assertion for matrix with only one column
    sys = numeric_t('[ 3, 1 ]') ;
    sys = sys' ;
    sd = numeric_t('[ 9/sqrt(10), 1/sqrt(10) ]') ;
    try
        [ ~ ] = gram_schmidt_mod( sys, sd ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'gram_schmidt_mod:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_gram_schmidt_mod_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_gram_schmidt_mod_R2_mp(testCase) % l1a
%      test basis of R^2, example taken from:
%      https://de.wikipedia.org/wiki/Gram-Schmidtsches_Orthogonalisierungsverfahren
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
    x1 = numeric_t('[ 3, 1 ]') ;
    x2 = numeric_t('[ 2, 2 ]') ;
    sys = [ x1; x2 ] ;
    sys = sys' ;
    sd = (x1.^2)/norm(x1) ;
    Q = gram_schmidt_mod( sys, sd ) ;
    u1 = numeric_t('[ 3/sqrt(10), 1/sqrt(10) ]') ;
    u2 = numeric_t('[ -1/sqrt(10), 3/sqrt(10) ]') ;
    orthosys = [ u1; u2 ] ;
    orthosys = orthosys' ;
    abstol = testCase.TestData.abstolfac * eps(numeric_t) ;
    reltol = testCase.TestData.reltolfac * eps(numeric_t) ;
    [ ~, ~, c ] = verifyAlmostEqual(Q, orthosys, abstol, ...
        reltol, true) ;
    verifyTrue(testCase, c)
end

function test_gram_schmidt_mod_R2_double(testCase) % l1b
%       test basis of R^2, example taken from:
%       https://de.wikipedia.org/wiki/Gram-Schmidtsches_Orthogonalisierungsverfahren
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'double')
    
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
    
    x1 = numeric_t('[ 3, 1 ]') ;
    x2 = numeric_t('[ 2, 2 ]') ;
    sys = [ x1; x2 ] ;
    sys = sys' ;
    sd = numeric_t('[ 0.5, 0.5 ]') ;
    Q = gram_schmidt_mod( sys, sd ) ;
    u1 = numeric_t('[ sqrt(0.5), sqrt(0.5) ]') ;
    u2 = numeric_t('[ sqrt(0.5), -sqrt(0.5) ]') ;
    orthosys = [ u1; u2 ] ;
    orthosys = orthosys' ;
    testCase.verifyThat(Q, IsEqualTo(orthosys, ...
       'Within', AbsoluteTolerance(abstol) & RelativeTolerance(reltol)))
end

function test_gram_schmidt_mod_R4(testCase) % l1
%       test for subspace of R^4, example taken from:
%       http://math.bard.edu/~mbelk/math601/GramSchmidtExamples.pdf
    x1 = numeric_t('[ 1, 1, 1, 1 ]') ;
    x2 = numeric_t('[ -1, 4, 4, 1 ]') ;
    x3 = numeric_t('[ 4, -2, 2, 0 ]') ;
    sys = [ x1; x2; x3 ] ;
    sys = sys' ;
    sd = numeric_t('[ 0.25, 0.25, 0.25, 0.25 ]') ;
    Q = gram_schmidt_mod( sys, sd ) ;
    u1 = numeric_t('[ 1/2, 1/2, 1/2, 1/2 ]') ;
    u2 = numeric_t('[ -1/sqrt(2), sqrt(2)/3, sqrt(2)/3, -1/(3*sqrt(2)) ]') ;
    u3 = numeric_t('[ 1/(2*sqrt(3)), -5/(6*sqrt(3)), 7/(6*sqrt(3)), -5/(6*sqrt(3)) ]') ;
    orthosys = [ u1; u2; u3 ] ;
    orthosys = orthosys' ;
    abstol = testCase.TestData.abstolfac * eps(numeric_t) ;
    reltol = testCase.TestData.reltolfac * eps(numeric_t) ;
    [ ~, ~, c ] = verifyAlmostEqual(Q, orthosys, abstol, ...
        reltol, true) ;
    verifyTrue(testCase, c)
end
% -------------------------------------------------------------------------
% test indexsearch (depends on numeric_t) (logic paths: 6) (assertions: 1)

function test_indexsearch_MatrixShapeError(testCase) % a1
%       test assertion for (N,k)-matrix with k>N
    v1 = numeric_t('[1,0,0,0]') ;
    v2 = numeric_t('[0,1,0,0]') ;
    v3 = numeric_t('[0,0,1,0]') ;
    sys = [ v1; v2; v3 ] ;
    try
        [ ~ ] = indexsearch( sys ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'indexsearch:MatrixShapeError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_indexsearch_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_indexsearch(testCase) % l
%   test correct function for standard simplex in 6D first
    v0 = numeric_t('[0,0,0,0,0,0]') ;
    v1 = numeric_t('[1,0,0,0,0,0]') ;
    v2 = numeric_t('[0,1,0,0,0,0]') ;
    v3 = numeric_t('[0,0,1,0,0,0]') ;
    v4 = numeric_t('[0,0,0,1,0,0]') ;
    v5 = numeric_t('[0,0,0,0,1,0]') ;
    v6 = numeric_t('[0,0,0,0,0,1]') ;
    sys = [ v0; v1; v0; v2; v0; v3; v0; v4; v0; v5; v0; v6 ] ;
%       find indices
    index = indexsearch( sys ) ;
    index_exp = numeric_t('[2,4,6,8,10,12]') ;
    dummy = ( index == index_exp ) ;
    verifyTrue(testCase, all(dummy(:)) )
%   define simple system of points in 3D, which all are inside (or on the
%   surface) of the 3D-simplex with vertices v0=[0,0,0], v1=[1.5,0,0],
%   v2=[0,2,0], v3=[0,0,3]
    v3 = numeric_t('[0,0,3]') ;
    p1 = numeric_t('[0.75,1,0]') ;
    v1 = numeric_t('[1.5,0,0]') ; 
    v0 = numeric_t('[0,0,0]') ; 
    p3 = numeric_t('[0.375,0.5,0.75]') ; 
    v2 = numeric_t('[0,2,0]') ;
    p2 = numeric_t('[0,1.2,1.2]') ; 
    p4 = numeric_t('[0.15,0.2,0.6]') ; 
    p5 = numeric_t('[0,0.6,0.3]') ;
    sys = [ v3; p1; v1; v0; p3; v2; p2; p4; p5 ] ;
%       find indices
    index = indexsearch( sys ) ;
    index_exp = numeric_t('[1,6,3]') ;
    dummy = ( index == index_exp ) ;
    verifyTrue(testCase, all(dummy(:)) )
end
% -------------------------------------------------------------------------
% test initialize_A (depends on do_schur, indexsearch, numeric_t) (logic paths: 4)
% (assertions: 4)

function test_initialize_A_FirstColumnError(testCase) % a1
%       test for the assertion if first column of schur vector matrix isnt
%       constantly equal 1
    d1 = numeric_t('[0,0,0]') ;
    d2 = numeric_t('[0,0,0]') ;
    d3 = numeric_t('[0,0,0]') ;
    d4 = numeric_t('[0,0,0]') ;
    dummy = [ d1; d2; d3; d4 ] ;
    try
        [ ~ ] = initialize_A( dummy, 1 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'initialize_A:FirstColumnError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_initialize_A_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_initialize_A_MatrixShapeError1(testCase) % a2
%       test assertion for (N,k)-matrix with k>N
    d1 = numeric_t('[1,0,0,0]') ;
    d2 = numeric_t('[1,0,0,0]') ;
    d3 = numeric_t('[1,0,0,0]') ;
    dummy = [ d1; d2; d3 ] ;
    try
        [ ~ ] = initialize_A( dummy, 1 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'initialize_A:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_initialize_A_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_initialize_A_MatrixShapeError2(testCase) % a3
%       test assertion for (N,k)-matrix with k=N
    d1 = numeric_t('[1,0,0]') ;
    d2 = numeric_t('[1,0,0]') ;
    d3 = numeric_t('[1,0,0]') ;
    dummy = [ d1; d2; d3 ] ;
    try
        [ ~ ] = initialize_A( dummy, 1 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'initialize_A:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_initialize_A_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_initialize_A_condition(testCase) % a4
%       test assertion for too high matrix condition
    dummy = hilb(14) ;
    dummy(:,14) = [] ;
    dummy(:,1) = numeric_t('1') ;
    try
        [ ~ ] = initialize_A( dummy, 1 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'initialize_A:A_ConditionError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_initialize_A_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_initialize_A_init(testCase) % l1
%       test for the error if init isnt 0 or 1
    d1 = numeric_t('[1,0,0]') ;
    d2 = numeric_t('[1,0,0]') ;
    d3 = numeric_t('[1,0,0]') ;
    d4 = numeric_t('[1,0,0]') ;
    dummy = [ d1; d2; d3; d4 ] ;
    init = 2 ;
    try
        [ ~ ] = initialize_A( dummy, init ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'initialize_A:InputError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_initialize_A_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_initialize_A(testCase) % l2
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       specify the test case data
    matrixfile = 'example_matrix_mu0.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_fillA' ;
%       get Schur vectors
    [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    evs = X(:,1:4) ;
%       initialize A
    [ T, A ] = evalc('initialize_A( evs, 1 )') ;
    fileID = fopen('test_initialize_A_report.txt','a') ;
    fprintf(fileID,'%s\n','test_initialize_A') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ; 
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       find indices
    index = indexsearch( evs ) ;
    A_exp = double( pinv( X(index,1:4) ) ) ;
    testCase.verifyThat(A, IsEqualTo(A_exp, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

function test_initialize_A_conditionWarning(testCase) % l3
%       test warning for high matrix condition
    dummy = hilb(6) ;
    dummy(:,6) = [] ;
    dummy(:,1) = numeric_t('1') ;
    [ T, ~ ] = evalc('initialize_A( dummy, 1 )') ;
    fileID = fopen('test_initialize_A_report.txt','a') ;
    fprintf(fileID,'%s\n','test_initialize_A_conditionWarning') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    warning = ['Warning for 5 clusters: The condition ' ...
                'number of the initial guess for A is > 1e4'] ;
    C = strsplit(T, '\n') ;
    verifyTrue( testCase, strcmp(C(1), warning) )
end

function test_initialize_A_loadOldfile(testCase) % l4
%       test initialization by loading of an already existing matrix from 
%       file.
%       set global variable to indicate testrun (no display output, no
%       interactive input): 'minChiON' indicates that minChi criterion is
%       used and input is taken from testInput_minChiON.mat. 'minChiOFF'
%       indicates that minChi criterion isnt used and input is taken from
%       testInput_minChiOFF.mat.
    global testmode
    testmode = 'minChiON' ;
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    dummy = [ d1; d2; d3 ; d4 ] ;
    actVal = initialize_A( dummy, 0 ) ;
    e1 = [0,1,0] ;
    e2 = [0,1,0] ;
    e3 = [0,1,0] ;
    e4 = [0,1,0] ;
    expVal = [ e1; e2; e3 ; e4 ] ;
    verifyEqual(testCase, actVal, expVal, 'AbsTol', eps)
end
% -------------------------------------------------------------------------
% test initialize_optional (depends on none) (logic paths: 5) (assertions: 0)

% test for valid initialization of iopt in case of empty iopt
function test_initialize_optional_empty_iopt(testCase) % l1
%       set the ID-string for file-naming
    id = 'test_initialize_optional_empty_iopt' ;
    [ iopt, fileid, final_opt ] = initialize_optional( [], id ) ;
    verifyEqual(testCase, isempty(iopt), true)
    verifyEqual(testCase, final_opt, false)
    verifyEqual(testCase, fileid, '')
end

function test_initialize_optional_SolverError(testCase)
    id = 'test_initialize_optional_SolverError' ;
    iopt.solver = 'Newton' ;
    try
        [ ~, ~, ~ ] = initialize_optional( iopt, id ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'initialize_optional:SolverError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_initialize_optional_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

% test for valid initialization of iopt in case of only 
% iopt.solver='gauss-newton' defined
function test_initialize_optional_gauss_newton(testCase) % l3
%       set the optional solver and ID-string for file-naming
    iopt.solver = 'gauss-newton' ;
    id = 'test_initialize_optional_gauss_newton' ;
    [ iopt, fileid, final_opt ] = initialize_optional( iopt, id ) ;
    verifyEqual(testCase, isstruct(iopt), true)
    verifyEqual(testCase, iopt.solver, 'gauss-newton')
    verifyEqual(testCase, iopt.init, 1)
    verifyEqual(testCase, iopt.parallel, 0)
    verifyEqual(testCase, iopt.maxiter, 100)
    verifyEqual(testCase, iopt.display, 0)
    verifyEqual(testCase, iopt.xtol, 1e-4)
    verifyEqual(testCase, iopt.xscale, 1e-06)
    fileid_exp = strcat(id,'-',iopt.solver,'-','maxiter','_', ...
        num2str(iopt.maxiter),'-','xscale','_', num2str(iopt.xscale), ...
        '-','xtol','_',num2str(iopt.xtol),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
    verifyEqual(testCase, final_opt, true)
end

% test for valid initialization of iopt in case of only 
% iopt.solver='nelder-mead' defined
function test_initialize_optional_nelder_mead(testCase) % l4
%       set the optional solver and ID-string for file-naming
    iopt.solver = 'nelder-mead' ;
    id = 'test_initialize_optional_nelder_mead' ;
    [ iopt, fileid, final_opt ] = initialize_optional( iopt, id ) ;
    verifyEqual(testCase, isstruct(iopt), true)
    verifyEqual(testCase, iopt.solver, 'nelder-mead')
    verifyEqual(testCase, iopt.init, 1)
    verifyEqual(testCase, iopt.parallel, 0)
    verifyEqual(testCase, iopt.maxiter, 2000)
    verifyEqual(testCase, iopt.display, 0)
    verifyEqual(testCase, iopt.tolfun, 1e-8)
    verifyEqual(testCase, iopt.tolx, 1e-8)
    fileid_exp = strcat(id,'-',iopt.solver,'-','maxiter','_', ...
        num2str(iopt.maxiter),'-','tolfun','_', num2str(iopt.tolfun), ...
        '-','tolx','_',num2str(iopt.tolx),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
    verifyEqual(testCase, final_opt, true)
end

% test for valid initialization of iopt in case of only 
% iopt.solver='levenberg-marquardt' defined
function test_initialize_optional_levenberg_marquardt(testCase) % l5
%       set the optional solver and ID-string for file-naming
    iopt.solver = 'levenberg-marquardt' ;
    id = 'test_initialize_optional_levenberg_marquardt' ;
    [ iopt, fileid, final_opt ] = initialize_optional( iopt, id ) ;
    verifyEqual(testCase, isstruct(iopt), true)
    verifyEqual(testCase, iopt.solver, 'levenberg-marquardt')
    verifyEqual(testCase, iopt.init, 1)
    verifyEqual(testCase, iopt.parallel, 0)
    verifyEqual(testCase, iopt.maxiter, 500)
    verifyEqual(testCase, iopt.display, 0)
    verifyEqual(testCase, iopt.tolfun, 1e-8)
    verifyEqual(testCase, iopt.tolx, 1e-10)
    fileid_exp = strcat(id,'-',iopt.solver,'-','maxiter','_', ...
        num2str(iopt.maxiter),'-','tolfun','_', num2str(iopt.tolfun), ...
        '-','tolx','_',num2str(iopt.tolx),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
    verifyEqual(testCase, final_opt, true)
end
% -------------------------------------------------------------------------
% test initialize_workspace (depends on none) (logic paths: 4) (assertions: 0)

% test for valid initialization of wk in case of empty wk
function test_initialize_workspace_empty_wk(testCase)
%       set the solver, optional solver and ID-string for file-naming
    [ wk, fileid ] = initialize_workspace( [] ) ;
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'gpcca')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'gauss-newton')
    verifyEqual(testCase, wk.maxiter, 100)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.xtol, 1e-4)
    verifyEqual(testCase, wk.xscale, 1e-06)
    fileid_exp = strcat(wk.id,'-',wk.solver,'-','maxiter','_', ...
        num2str(wk.maxiter),'-','xscale','_', num2str(wk.xscale), ...
        '-','xtol','_',num2str(wk.xtol),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
end

% test for valid initialization of wk in case of only 
% wk.solver='gauss-newton' defined
function test_initialize_workspace_gauss_newton(testCase)
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'gauss-newton' ;
    wk.id = 'test_initialize_workspace_gauss_newton' ;
    [ wk, fileid ] = initialize_workspace( wk ) ;
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_initialize_workspace_gauss_newton')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'gauss-newton')
    verifyEqual(testCase, wk.maxiter, 100)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.xtol, 1e-4)
    verifyEqual(testCase, wk.xscale, 1e-06)
    fileid_exp = strcat(wk.id,'-',wk.solver,'-','maxiter','_', ...
        num2str(wk.maxiter),'-','xscale','_', num2str(wk.xscale), ...
        '-','xtol','_',num2str(wk.xtol),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
end

% test for valid initialization of wk in case of only 
% wk.solver='nelder-mead' defined
function test_initialize_workspace_nelder_mead(testCase)
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'nelder-mead' ;
    wk.id = 'test_initialize_workspace_nelder_mead' ;
    [ wk, fileid ] = initialize_workspace( wk ) ;
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_initialize_workspace_nelder_mead')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'nelder-mead')
    verifyEqual(testCase, wk.maxiter, 2000)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.tolfun, 1e-4)
    verifyEqual(testCase, wk.tolx, 1e-4)
    fileid_exp = strcat(wk.id,'-',wk.solver,'-','maxiter','_', ...
        num2str(wk.maxiter),'-','tolfun','_', num2str(wk.tolfun), ...
        '-','tolx','_',num2str(wk.tolx),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
end

% test for valid initialization of wk in case of only 
% wk.solver='levenberg-marquardt' defined
function test_initialize_workspace_levenberg_marquardt(testCase)
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'levenberg-marquardt' ;
    wk.id = 'test_initialize_workspace_levenberg_marquardt' ;
    [ wk, fileid ] = initialize_workspace( wk ) ;
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_initialize_workspace_levenberg_marquardt')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'levenberg-marquardt')
    verifyEqual(testCase, wk.maxiter, 500)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.tolfun, 1e-8)
    verifyEqual(testCase, wk.tolx, 1e-10)
    fileid_exp = strcat(wk.id,'-',wk.solver,'-','maxiter','_', ...
        num2str(wk.maxiter),'-','tolfun','_', num2str(wk.tolfun), ...
        '-','tolx','_',num2str(wk.tolx),'-',numeric_t) ;
    verifyEqual(testCase, fileid, fileid_exp)
end
% -------------------------------------------------------------------------
% test inputT (depends on none) (logic paths: 3) (assertions: 0)

function test_inputT(testCase) % l1 + l2 + l3
    global testmode
%       test silent input from file during testmode with use of minChi criterion
    testmode = 'minChiON' ;
    minChi_switch = inputT( 'dummy', 'minChi_switch') ;
    verifyEqual(testCase,minChi_switch,1)
%       test silent input from file during testmode without use of minChi criterion
    testmode = 'minChiOFF' ;
    minChi_switch = inputT( 'dummy', 'minChi_switch') ;
    verifyEqual(testCase,minChi_switch,0)
    clearvars -global testmode
%       test verbose standard input with testmode not defined
    minChi_switch = inputT( ['This is a test! Please type 1 to verify the ',...
        'correct function: '], 'minChi_switch' ) ;
    verifyEqual(testCase,minChi_switch,1)
%       test verbose standard input with testmode invalidly defined
    testmode = [] ;
    minChi_switch = inputT( ['This is a test! Please type 1 to verify the ',...
        'correct function: '], 'minChi_switch' ) ;
    verifyEqual(testCase,minChi_switch,1)
    clearvars -global testmode
end
% -------------------------------------------------------------------------
% test load_t (depends on none) (logic paths: 3) (assertions: 0)

function test_load_t_InputError(testCase) % l1
    name = 'count-sd_weights.txt' ;
    prec = 'single' ;
    try
        [ ~ ] = load_t( name, '-ascii', prec ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'load_t:InputError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_load_t_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_load_t_mp(testCase) % l2
%       set the precision and specify the test case data
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
    sd = testCase.TestData.sd ;
    name = 'count-sd.txt' ;
    r = load_t( name, '-ascii', prec ) ;
    abstol = eps(numeric_t) ;
    reltol = eps(numeric_t) ;
    [ ~, ~, c ] = verifyAlmostEqual(r, sd, abstol, reltol, true) ;
    verifyTrue(testCase, c)
end

function test_load_t_double(testCase) % l3
%       set the precision and specify the test case data
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'double')
    sd_weights = testCase.TestData.sd_weights ;
    name = 'count-sd_weights.txt' ;
    r = load_t( name, '-ascii', prec ) ;
    abstol = eps(numeric_t) ;
    reltol = eps(numeric_t) ;
    [ ~, ~, c ] = verifyAlmostEqual(r, sd_weights, abstol, reltol, true) ;
    verifyTrue(testCase, c)
end 

function test_load_t_double_InputError(testCase) % l4
%       set the precision and specify the test case data
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'double')
    name = 'count-sd_weights.txt' ;
    try
        [ ~ ] = load_t( name, '-mat', prec ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'load_t:InputError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_load_t_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end 
% -------------------------------------------------------------------------
% test main_nlscon with nlscon and problem_pcca_nlscon

% So far no tests for nlscon, main_nlscon and problem_pcca_nlscon
% except test_opt_soft_gauss_newton()

%--------------------------------------------------------------------------
% test nearly_equal (depends on numeric_t) (logic paths: 8) (assertions: 2)

function test_nearly_equal_DataTypeError(testCase) % a1 + a2
%       test the case that a isnt double or mp
    a = single(2.0) ;
    b = numeric_t('2.0') ;
    try
        [ ~, ~, ~ ] = nearly_equal( a, b ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'nearly_equal:a_DataTypeError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_nearly_equal_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
%       test the case that b isnt double or mp    
    a = numeric_t('2.0') ;
    b = single(2.0) ;
    try
        [ ~, ~, ~ ] = nearly_equal( a, b ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'nearly_equal:b_DataTypeError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_nearly_equal_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_nearly_equal_mixed(testCase) % l1
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
    a = ( numeric_t('2.0') + 10 * eps(numeric_t) ) ;
    b = double(2.0) ;
    [ T, a, b, c ] = evalc('nearly_equal(a, b)') ;
 %       T is stored to a file
    fileID = fopen('test_nearly_equal_report.txt','a') ;
    fprintf(fileID,'%s\n','test_nearly_equal_mixed') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, 'double'))
    verifyTrue(testCase, isa(b, 'double'))
end

function test_nearly_equal(testCase) % l2
    a = numeric_t('1.0') ;
    b = ( numeric_t('1.0') + eps(numeric_t) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('1.0') ;
    b = ( numeric_t('1.0') + 2*eps(numeric_t) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('1.0') ;
    b = ( numeric_t('1.0') + 3*eps(numeric_t) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('2.0') ;
    b = ( numeric_t('2.0') + 4*eps(numeric_t) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('2.0') ;
    b = ( numeric_t('2.0') + 6*eps(numeric_t) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_nearly_equal_exact(testCase) % l3a
    a = numeric_t('2.0') ;
    b = numeric_t('2.0') ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_nearly_equal_inf(testCase) % l3b
    a = inf ;
    b = inf ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyEqual(testCase, a, inf)
    verifyEqual(testCase, b, inf)
end

function test_nearly_equal_zero(testCase) % l4
    a = numeric_t('0.0') ;
    b = numeric_t('0.0') ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

% case b=0 isnt tested, since its equivalent to a=0, 
% case diff < realmin(num_t) is here only tested in combination with case
% a=0.
function test_nearly_equal_mixedZero_mp(testCase) % l4&l6 + l8
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
    a = numeric_t('0.0') ;
    b = ( numeric_t('0.0') + realmin(numeric_t) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('0.0') ;
    b = ( numeric_t('0.0') + realmin(numeric_t)/numeric_t('10') ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_nearly_equal_mixedZero_double(testCase) % l4&l6 + l7
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'double')
    a = 0.0 ;
    b = ( 0.0 + (eps * realmin) ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, prec))
    verifyTrue(testCase, isa(b, prec))
    
    a = 0.0 ;
    b = ( 0.0 + (eps * realmin)/10 ) ;
    [ a, b, c ] = nearly_equal(a, b) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, prec))
    verifyTrue(testCase, isa(b, prec))
end
% -------------------------------------------------------------------------
% test numeric_t (depends on none) (logic paths: 4) (assertions: 0)

function test_numeric_t_empty(testCase) % l1
%       test empty case (1)
    prec = testCase.TestData.precision ;
    verifyEqual(testCase, numeric_t, prec)
end

function test_numeric_t_double(testCase) % l2 + l3
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'double')
%       test numeric case (3)
    dummy = numeric_t('2.1') ;
    verifyTrue(testCase, isa(dummy, prec))
    dummy1 = double(2.1) ;
    [ ~, ~, c ] = nearly_equal(dummy1, dummy) ;
    verifyTrue(testCase, c)
%       test expression case (4)
    dummy = numeric_t('10.0/2.0') ;
    verifyTrue(testCase, isa(dummy, prec))
    dummy1 = double(10.0)/double(2.0) ;
    [ ~, ~, c ] = nearly_equal(dummy1, dummy) ;
    verifyTrue(testCase, c)
end

function test_numeric_t_mp(testCase) % l4
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
%       test numeric case (2a)
    dummy = numeric_t('2.1') ;
    verifyTrue(testCase, isa(dummy, prec))
    dummy1 = mp('2.1') ;
    [ ~, ~, c ] = nearly_equal(dummy1, dummy) ;
    verifyTrue(testCase, c)
%       test expression case (2b)
    dummy = numeric_t('10.0/2.0') ;
    verifyTrue(testCase, isa(dummy, prec))
    dummy1 = mp('10.0/2.0') ;
    [ ~, ~, c ] = nearly_equal(dummy1, dummy) ;
    verifyTrue(testCase, c)
end
% -------------------------------------------------------------------------
% test objective (depends on none) (logic paths: 3) (assertions: 6)
    
function test_objective_DataTypeError1(testCase) % a1
%       test assertion, if schur vectors aren't of type double
    d1 = single([1,0,0]) ;
    d2 = single([1,0,0]) ;
    d3 = single([1,0,0]) ;
    d4 = single([1,0,0]) ;
    svecs = [ d1; d2; d3; d4 ] ;
    alpha = zeros(1,4) ;
    optfile = '' ;
    try
        [ ~ ] = objective(alpha, svecs, optfile) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'objective:DataTypeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_objective_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_objective_DataTypeError2(testCase) % a2
%       test assertion, if alpha vector isn't of type double
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    alpha = single(zeros(1,4)) ;
    optfile = '' ;
    try
        [ ~ ] = objective(alpha, svecs, optfile) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'objective:DataTypeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_objective_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_objective_MatrixShapeError1(testCase) % a3
%       test assertion for schur vector (N,k)-matrix  with k>N
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    alpha = zeros(1,9) ;
    optfile = '' ;
    try
        [ ~ ] = objective(alpha, svecs, optfile) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'objective:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_objective_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end


function test_objective_MatrixShapeError2(testCase) % a4
%       test assertion for schur vector (N,k)-matrix with k=N
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    alpha = zeros(1,4) ;
    optfile = '' ;
    try
        [ ~ ] = objective(alpha, svecs, optfile) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'objective:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_objective_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_objective_FirstColumnError(testCase) % a5
%       test for the assertion if first column of schur vector matrix isnt
%       constantly equal 1
    d1 = [0,0,0] ;
    d2 = [0,0,0] ;
    d3 = [0,0,0] ;
    d4 = [0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    alpha = zeros(1,4) ;
    optfile = '' ; 
    try
        [ ~ ] = objective(alpha, svecs, optfile) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'objective:FirstColumnError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_objective_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_objective_MatrixShapeError3(testCase) % a6
%       test assertion for schur vector (N,k)-matrix with k=N
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    alpha = zeros(1,3) ;
    optfile = '' ;
    try
        [ ~ ] = objective(alpha, svecs, optfile) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'objective:MatrixShapeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_objective_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_objective1(testCase) % l
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
    k = 3 ;
%       calculate actual A by objective()
    svecs = testCase.TestData.svecs_mu0 ;
    A_init = testCase.TestData.A_mu0_init ;
    alpha = zeros(1,(k-1)^2) ;
        for i = 1:k-1
            for j = 1:k-1
                alpha(j + (i-1)*(k-1)) = A_init(i+1,j+1) ;
            end
        end
    name = 'test_objective1-optfile.txt' ;
    optfile = fopen(name, 'w') ;
    actVal = objective(alpha, svecs, optfile) ;
    fclose(optfile) ;
%       calculate expected A matrix by formula
    A = testCase.TestData.A_mu0 ;
    expVal = k - sum( sum(A.^2,1)./A(1,:) ) ;
    testCase.verifyThat(actVal, IsEqualTo(expVal, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end 

function test_objective2(testCase) % l
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
    k = 5 ;
%       calculate actual A by objective()
    svecs = testCase.TestData.svecs_mu1000 ;
    A_init = testCase.TestData.A_mu1000_init ;
    alpha = zeros(1,(k-1)^2) ;
        for i = 1:k-1
            for j = 1:k-1
                alpha(j + (i-1)*(k-1)) = A_init(i+1,j+1) ;
            end
        end
    name = 'test_objective2-optfile.txt' ;
    optfile = fopen(name, 'w') ;
    actVal = objective(alpha, svecs, optfile) ;
    fclose(optfile) ;
%       calculate expected A matrix by formula
    A = testCase.TestData.A_mu1000 ;
    expVal = k - sum( sum(A.^2,1)./A(1,:) ) ;
    testCase.verifyThat(actVal, IsEqualTo(expVal, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end 
% -------------------------------------------------------------------------
% test opt_soft (i.e. with an exponential function?) (logic paths: 23)
% (assertions: 10)

function test_opt_soft_DataTypeError1(testCase) % a1
%       test assertion, if A matrix isn't of type double
    A = single(zeros(3)) ;
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 3 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:DataTypeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_DataTypeError2(testCase) % a2
%       test assertion, if schur vector matrix isn't of type double
    A = zeros(3) ;
    d1 = single([1,0,0]) ;
    d2 = single([1,0,0]) ;
    d3 = single([1,0,0]) ;
    d4 = single([1,0,0]) ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 3 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:DataTypeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_DataTypeError3(testCase) % a3
%       test assertion, if fileid isn't of type char
    A = zeros(3) ;
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 3 ;
    fileid = 1 ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:DataTypeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_DataTypeError4(testCase) % a4
%       test assertion, if wk isn't of type struct
    A = zeros(3) ;
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 3 ;
    fileid = '' ;
    wk = 1 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:DataTypeError4'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_MatrixShapeError1(testCase) % a5
%       test assertion for the case that A isn't quadratic
    A = zeros(3) ;
    A(:,3) = [] ;
    d1 = [1,0,0] ;
    d2 = [1,0,0] ;
    d3 = [1,0,0] ;
    d4 = [1,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 3 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_MatchError(testCase) % a6
%       test assertion if the second dimension of the schur vector matrix
%       doesn't match with the dimensions of A
    A = zeros(3) ;
    d1 = [1,0] ;
    d2 = [1,0] ;
    d3 = [1,0] ;
    d4 = [1,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 2 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:MatchError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_MatrixShapeError2(testCase) % a7
%       test assertion for the case that A isn't at least a (2x2)-matrix
    A = 0 ;
    d1 = 1 ;
    d2 = 1 ;
    d3 = 1 ;
    d4 = 1 ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 1 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_MatrixShapeError3(testCase) % a8
%       test assertion for schur vector (N,k)-matrix  with k>N
    A = zeros(4) ;
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    counter = 4 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:MatrixShapeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_MatrixShapeError4(testCase) % a9
%       test assertion for schur vector (N,k)-matrix with k=N
    A = zeros(4) ;
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    d4 = [1,0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 4 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] = opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:MatrixShapeError4'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_FirstColumnError(testCase) % a10
%       test for the assertion if first column of schur vector matrix isnt
%       constantly equal 1
    A = zeros(3) ;
    d1 = [0,0,0] ;
    d2 = [0,0,0] ;
    d3 = [0,0,0] ;
    d4 = [0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    counter = 3 ;
    fileid = '' ;
    wk.id = 'gpcca' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xtol = 1e-4 ;
    wk.xscale = 1e-06 ;
    try
        [ ~, ~, ~ ] =opt_soft( A, svecs, counter, fileid, wk ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'opt_soft:FirstColumnError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_opt_soft_nelder_mead_mu0(testCase) % l5
% test the case, that nelder-mead is the choosen solver using the
% example matrix with mu=0
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
%       set variables
    A = testCase.TestData.A_mu0 ;
    svecs = testCase.TestData.svecs_mu0 ;
    counter = 3 ;
    fileid = 'test_opt_soft_nelder_mead_mu0' ;
    wk.id = 'test_opt_soft_nelder_mead_mu0' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'nelder-mead' ;
    wk.maxiter = -1 ;
    wk.display = 0 ;
    wk.tolfun = 1e-8 ;
    wk.tolx = 1e-8 ;
    [ T, chi, A, fopt ] = evalc('opt_soft( A, svecs, counter, fileid, wk )') ;
%       close all open (but invisible) plots
    close all
    name = 'test_opt_soft_nelder_mead_mu0-A.txt' ;
    save(name, 'A', '-ascii', '-double')
    name = 'test_opt_soft_nelder_mead_mu0-chi.txt' ;
    save(name, 'chi', '-ascii', '-double')
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s\n','test_opt_soft_nelder_mead_mu0') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    crispness = ( 3 - fopt ) / 3 ;
%       The crispness to compare with was taken from Fig.4 d) on page 168 
%       of (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
%       Application to Markov State Models and Data Classification. Adv. 
%       Data Anal. Classif. 2013, 7 (2), 147-179.
    verifyEqual(testCase, crispness, 0.973, 'AbsTol', 0.001)
end

function test_opt_soft_nelder_mead_mu1000(testCase) % l6+l7
% test the case, that nelder-mead is the choosen solver using the
% example matrix with mu=1000
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
%       set variables
    A = testCase.TestData.A_mu1000 ;
    svecs = testCase.TestData.svecs_mu1000 ;
    counter = 5 ;
    fileid = 'test_opt_soft_nelder_mead_mu1000' ;
    wk.id = 'test_opt_soft_nelder_mead_mu1000' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'nelder-mead' ;
    wk.maxiter = 5000 ;
    wk.display = 0 ;
    wk.tolfun = 1e-8 ;
    wk.tolx = 1e-8 ;
    [ T, chi, A, fopt ] = evalc('opt_soft( A, svecs, counter, fileid, wk )') ;
%       close all open (but invisible) plots
    close all
    name = 'test_opt_soft_nelder_mead_mu1000-A.txt' ;
    save(name, 'A', '-ascii', '-double')
    name = 'test_opt_soft_nelder_mead_mu1000-chi.txt' ;
    save(name, 'chi', '-ascii', '-double')
    fileID = fopen('test_opt_soft_report.txt','a') ;
    fprintf(fileID,'%s\n','test_opt_soft_nelder_mead_mu1000') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    crispness = ( 5 - fopt ) / 5 ;
%       The crispness to compare with was taken from page 168 of
%       (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
%       Application to Markov State Models and Data Classification. Adv. 
%       Data Anal. Classif. 2013, 7 (2), 147-179.
    verifyEqual(testCase, crispness, 0.804, 'AbsTol', 0.0025)
end

function test_opt_soft_nelder_mead_more(testCase) % l6+l7
% test the case, that nelder-mead is the choosen solver using the
% example matrices with mu=0,10,50,100,200,500,1000.
% The results to compare with are taken from page 168 of
% (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
% Application to Markov State Models and Data Classification. Adv. Data
% Anal. Classif. 2013, 7 (2), 147-179.
% There is one difference: In (1) the optimal cluster numbers for 
% mu=0,10,50,100,200,500,1000 are kopt=3,3,3,3,2,2,5 while here the
% correct optimal cluster numbers are assumed to be kopt=3,3,3,3,2,2,7
% since in all tests we obtained for mu=1000 an optimal cluster number of
% kopt=7.
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
%       set variables
    kmin = 2 ;
    kmax = 8 ;
    fileID = fopen('test_opt_soft_nelder_mead_more_report.txt','w') ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'nelder-mead' ;
    wk.maxiter = 2000 ;
    wk.display = 0 ;
    wk.tolfun = 1e-8 ;
    wk.tolx = 1e-8 ;
    kopt_vec = [] ;
%       iterate through all mu's
    i = 0 ;
    for mu = [ 0 10 50 100 200 500 1000 ]
        i = i + 1 ;
        matrixfilename = strcat('example_matrix_mu',num2str(mu)) ;
        matrixfile = strcat(matrixfilename,'.txt') ;
%           get the known test case data
        [ P, sd ] = get_knownInput( matrixfile ) ;
        fileid = strcat('test_opt_soft_nelder_mead_more-mu',num2str(mu)) ;
%           calculate schur vectors
        [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
        crisp_vec = [] ;
%           iterate through all k's
        for k = kmin:kmax
            svecs = X(:,1:k) ;
%               initialize A
            [ T, A ] = evalc('initialize_A( svecs, 1 )') ;
            name = strcat(fileid,'-n=',num2str(k),'-A_init.txt') ;
            save(name, 'A', '-ascii', '-double')
%               T is stored to a file
            fprintf(fileID,'%s\n\n',T) ;
            wk.id = fileid ;
            counter = k ;
%               make schur vectors double (needed by opt_soft)
            svecs = double(svecs) ;
            [ T1, chi, A, fopt ] = evalc('opt_soft(A,svecs,counter,fileid,wk)') ;
            name = strcat(fileid,'-n=',num2str(k),'-A.txt') ;
            save(name, 'A', '-ascii', '-double')
            name = strcat(fileid,'-n=',num2str(k),'-chi.txt') ;
            save(name, 'chi', '-ascii', '-double')
%               close all open (but invisible) plots
            close all
%               T1 is stored to a file
            fprintf(fileID,'%s\n',T1) ;
            crisp_vec(k-kmin+1,1) = k ;
            crisp_vec(k-kmin+1,2) = (k-fopt)/k ;
            [ ~ , dummy ] = max(flipud(crisp_vec(:,2))) ;
            kopt = ((kmax - dummy) + 1) ;
            kopt_vec(i) = kopt ;
        end
        name = strcat(fileid,'-crisp_vec.txt') ;
        save(name,'crisp_vec','-ascii','-double')
    end 
    name = 'test_opt_soft_nelder_mead_more-kopt_vec.txt' ;
        save(name,'kopt_vec','-ascii','-double')
    fclose(fileID) ;
    verifyEqual(testCase,kopt_vec,[3,3,3,3,2,2,7])
end

function test_opt_soft_nelder_mead_display_mu0(testCase) % l8+l9
% test the display option in the case, that nelder-mead is the choosen 
% solver using the example matrix with mu=0
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
%       set variables
    A = testCase.TestData.A_mu0 ;
    svecs = testCase.TestData.svecs_mu0 ;
    counter = 3 ;
    fileid = 'test_opt_soft_nelder_mead_display_mu0' ;
    wk.id = 'test_opt_soft_nelder_mead_display_mu0' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'nelder-mead' ;
    wk.maxiter = -1 ;
    wk.display = 1 ;
    wk.tolfun = 1e-8 ;
    wk.tolx = 1e-8 ;
    [ T, chi, A, fopt ] = evalc('opt_soft( A, svecs, counter, fileid, wk )') ;
%       close all open (but invisible) plots
    close all
    name = 'test_opt_soft_nelder_mead_display_mu0-A.txt' ;
    save(name, 'A', '-ascii', '-double')
    name = 'test_opt_soft_nelder_mead_display_mu0-chi.txt' ;
    save(name, 'chi', '-ascii', '-double')
    fileID = fopen('test_opt_soft_nelder_mead_display_mu0_report.txt','w') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    crispness = ( 3 - fopt ) / 3 ;
    verifyEqual(testCase, crispness, 0.973, 'AbsTol', 0.001)
    name = strcat(fileid,'-','n=',num2str(counter),'-','optimization','.fig') ;
    verifyTrue(testCase, exist(name,'file')==2)
end

function test_opt_soft_nelder_mead_display_mu1000(testCase) % l10+l11
% test the display option in the case, that nelder-mead is the choosen 
% solver using the example matrix with mu=1000
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
%       set variables
    A = testCase.TestData.A_mu1000 ;
    svecs = testCase.TestData.svecs_mu1000 ;
    counter = 5 ;
    fileid = 'test_opt_soft_nelder_mead_display_mu1000' ;
    wk.id = 'test_opt_soft_nelder_mead_display_mu1000' ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'nelder-mead' ;
    wk.maxiter = 5000 ;
    wk.display = 1 ;
    wk.tolfun = 1e-8 ;
    wk.tolx = 1e-8 ;
    [ T, chi, A, fopt ] = evalc('opt_soft( A, svecs, counter, fileid, wk )') ;
%       close all open (but invisible) plots
    close all
    name = 'test_opt_soft_nelder_mead_display_mu1000-A.txt' ;
    save(name, 'A', '-ascii', '-double')
    name = 'test_opt_soft_nelder_mead_display_mu1000-chi.txt' ;
    save(name, 'chi', '-ascii', '-double')
    fileID = fopen('test_opt_soft_nelder_mead_display_mu1000_report.txt','w') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    crispness = ( 5 - fopt ) / 5 ;
%       The crispness to compare with was taken from page 168 of
%       (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
%       Application to Markov State Models and Data Classification. Adv. 
%       Data Anal. Classif. 2013, 7 (2), 147-179.
    verifyEqual(testCase, crispness, 0.804, 'AbsTol', 0.0025)
    name = strcat(fileid,'-','n=',num2str(counter),'-','optimization','.fig') ;
    verifyTrue(testCase, exist(name,'file')==2)
end

function test_opt_soft_gauss_newton(testCase) % l20+l21
% test the case, that gauss-newton is the choosen solver using the
% example matrices with mu=0,10,50,100,200,500,1000.
% The results to compare with are taken from page 168 of
% (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
% Application to Markov State Models and Data Classification. Adv. Data
% Anal. Classif. 2013, 7 (2), 147-179.
% There is one difference: In (1) the optimal cluster numbers for 
% mu=0,10,50,100,200,500,1000 are kopt=3,3,3,3,2,2,5 while here the
% correct optimal cluster numbers are assumed to be kopt=3,3,3,3,2,2,7
% since in all tests we obtained for mu=1000 an optimal cluster number of
% kopt=7.

%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
    kmin = 2 ;
    kmax = 8 ;
    fileID = fopen('test_opt_soft_gauss_newton_report.txt','w') ;
    wk.schur = 1 ;
    wk.init = 1 ;
    wk.solver = 'gauss-newton' ;
    wk.maxiter = 100 ;
    wk.display = 0 ;
    wk.xscale = 1e-6 ;
    wk.xtol = 1e-4 ;
    kopt_vec = [] ;
    i = 0 ;
    for mu = [ 0 10 50 100 200 500 1000 ]
        i = i + 1 ;
        matrixfilename = strcat('example_matrix_mu',num2str(mu)) ;
        matrixfile = strcat(matrixfilename,'.txt') ;
%           get the known test case data
        [ P, sd ] = get_knownInput( matrixfile ) ;
        fileid = strcat('test_opt_soft_gauss_newton-mu',num2str(mu)) ;
%           calculate schur vectors
        [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
        crisp_vec = [] ;
        for k = kmin:kmax
            svecs = X(:,1:k) ;
%               initialize A
            [ T, A ] = evalc('initialize_A( svecs, 1 )') ;
            name = strcat(fileid,'-n=',num2str(k),'-A_init.txt') ;
            save(name, 'A', '-ascii', '-double')
%               T is stored to a file
            fprintf(fileID,'%s\n\n',T) ;
            wk.id = fileid ;
            counter = k ;
            svecs = double(svecs) ;
            [ T1, chi, A, fopt ] = evalc('opt_soft(A,svecs,counter,fileid,wk)') ;
            name = strcat(fileid,'-n=',num2str(k),'-A.txt') ;
            save(name, 'A', '-ascii', '-double')
            name = strcat(fileid,'-n=',num2str(k),'-chi.txt') ;
            save(name, 'chi', '-ascii', '-double')
%               close all open (but invisible) plots
            close all
%               T1 is stored to a file
            fprintf(fileID,'%s\n',T1) ;
            crisp_vec(k-kmin+1,1) = k ;
            crisp_vec(k-kmin+1,2) = (k-fopt)/k ;
            [ ~ , dummy ] = max(flipud(crisp_vec(:,2))) ;
            kopt = ((kmax - dummy) + 1) ;
            kopt_vec(i) = kopt ;
        end
        name = strcat(fileid,'-crisp_vec.txt') ;
        save(name,'crisp_vec','-ascii','-double')
    end 
    name = 'test_opt_soft_gauss_newton-kopt_vec.txt' ;
        save(name,'kopt_vec','-ascii','-double')
    fclose(fileID) ;
    verifyEqual(testCase,kopt_vec,[3,3,3,3,2,2,7])
end

% add more test cases for levenberg-marquardt, mb give optimset back to 
% check, if options were correctly set...

% -------------------------------------------------------------------------
% test optimization_loop (logic paths: 8) (assertions: 5) (depends on
% do_schur, initialize_A and initialize_workspace)
    
function test_optimization_loop_DataTypeError1(testCase) %a1
    parameters.parallel = 0 ;
    P = ones(4) ;
    sd = ones(4,1) ;
    EVS = ones(4) ;
    A_cell = eye(4,2) ;
    kmin = 2 ;
    kmax = 3 ;
    fileid = 'test_optimization_loop_DataTypeError1' ;
    try
        [ ~, ~, ~, ~, ~ ] =optimization_loop( P,sd, EVS, A_cell, kmin, ...
            kmax, parameters, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'optimization_loop:DataTypeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_optimization_loop_report.txt','w') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_DataTypeError1') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_optimization_loop_DataTypeError2(testCase) %a2
    parameters.parallel = 0 ;
    P = ones(4) ;
    sd = ones(4,1) ;
    EVS = single(ones(4)) ;
    A_cell = { eye(4,2), eye(4,3) } ;
    kmin = 2 ;
    kmax = 3 ;
    fileid = 'test_optimization_loop_DataTypeError2' ;
    try
        [ ~, ~, ~, ~, ~ ] =optimization_loop( P,sd, EVS, A_cell, kmin, ...
            kmax, parameters, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'optimization_loop:DataTypeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_DataTypeError2') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_optimization_loop_DataTypeError3(testCase) %a3
    parameters.parallel = 0 ;
    P = ones(4) ;
    sd = ones(4,1) ;
    EVS = ones(4) ;
    A_cell = { eye(4,2), eye(4,3) } ;
    kmin = 2 ;
    kmax = 3 ;
    fileid = 3 ;
    try
        [ ~, ~, ~, ~, ~ ] =optimization_loop( P,sd, EVS, A_cell, kmin, ...
            kmax, parameters, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'optimization_loop:DataTypeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_DataTypeError3') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_optimization_loop_DataTypeError4(testCase) %a4
    parameters = 0 ;
    P = ones(4) ;
    sd = ones(4,1) ;
    EVS = ones(4) ;
    A_cell = { eye(4,2), eye(4,3) } ;
    kmin = 2 ;
    kmax = 3 ;
    fileid = 'test_optimization_loop_DataTypeError4' ;
    try
        [ ~, ~, ~, ~, ~ ] =optimization_loop( P,sd, EVS, A_cell, kmin, ...
            kmax, parameters, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'optimization_loop:DataTypeError4'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_DataTypeError4') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_optimization_loop_FirstColumnError(testCase) %a5
    parameters.parallel = 0 ;
    P = ones(4) ;
    sd = ones(4,1) ;
    EVS = eye(4,4) ;
    A_cell = { eye(4,2), eye(4,3) } ;
    kmin = 2 ;
    kmax = 3 ;
    fileid = 'test_optimization_loop_FirstColumnError' ;
    try
        [ ~, ~, ~, ~, ~ ] =optimization_loop( P, sd, EVS, A_cell, kmin, ...
            kmax, parameters, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'optimization_loop:FirstColumnError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_FirstColumnError') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_optimization_loop_PoolInitializationError(testCase) %l1
    parameters.parallel = 2 ;
    P = ones(4) ;
    sd = ones(4,1) ;
    EVS = ones(4) ;
    A_cell = { eye(4,2), eye(4,3) } ;
    kmin = 2 ;
    kmax = 3 ;
    fileid = 'test_optimization_loop_PoolInitializationError' ;
    try
        [ ~, ~, ~, ~, ~ ] =optimization_loop( P,sd, EVS, A_cell, kmin, ...
            kmax, parameters, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'optimization_loop:PoolInitializationError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_PoolInitializationError') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_optimization_loop_serial_nelderMead(testCase) %l2 + l4
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'nelder-mead' ;
    wk.id = 'test_optimization_loop_serial_nelderMead' ;
    wk.maxiter = 2000 ;
    wk.tolfun = 1e-08 ;
    wk.tolx = 1e-08 ;
%       initialize wk
    [ T, wk, fileid ] = evalc('initialize_workspace( wk )') ;
%       T  is stored to a file
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_serial_nelderMead') ;
    fprintf(fileID,'%s\n\n\n',T) ;
%       verify that wk is fine
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_optimization_loop_serial_nelderMead')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'nelder-mead')
    verifyEqual(testCase, wk.parallel, 0)
    verifyEqual(testCase, wk.maxiter, 2000)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.tolfun, 1e-8)
    verifyEqual(testCase, wk.tolx, 1e-8)
%       specify the test case data
    matrixfile = 'count.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
%       check if the sd calculated by get_knownInput() matches the trusted 
%       testCase.TestData.sd
    [ ~, ~, c ] = verifyAlmostEqual(sd, testCase.TestData.sd, abstol, ...
        reltol, false) ;
    verifyTrue(testCase, c)
%           calculate schur vectors
    [ EVS, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    kmin = 2 ;
    kmax = 3 ;
%       initialize A
    A_cell = cell([1,kmax]) ;
    for kt = kmin:kmax
        evs = EVS(:,1:kt) ;
        [ T, A_cell{kt} ] = evalc('initialize_A(evs, wk.init)') ;
        fprintf(fileID,'iteration %d\n',kt) ;
        fprintf(fileID,'%s\n\n\n',T) ;
    end
%       run optimization loop
    [ T, Pc, A_cell, chi, val_vec, opt_vec ] = evalc(['optimization_loop('...
    'double(P), double(sd), double(EVS), A_cell, kmin, kmax, wk, fileid )']) ;
%       T  is stored to a file
    fprintf(fileID,'%s\n\n\n',T) ;
    fclose(fileID) ;
%       ensure that A_cell is a cell
    verifyTrue(testCase,iscell(A_cell))
%       test if actual Pc, chi, A are equal (within numerical 
%       errors) to the "good" (correct) values 
    testCase.verifyThat(Pc, IsEqualTo(testCase.TestData.Pc, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(chi, IsEqualTo(testCase.TestData.chi, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))    
    testCase.verifyThat(A_cell{end}, IsEqualTo(testCase.TestData.A, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       load the saved actual Pc, chi, A from files
    fileid = ['test_optimization_loop_serial_nelderMead-nelder-mead-'...
        'maxiter_2000-tolfun_1e-08-tolx_1e-08'] ;
    fileid = strcat(fileid,'-',testCase.TestData.precision) ;
    Pc_saved = load(strcat(fileid,'-n=3-Pc.txt'),'-ascii') ;
    chi_saved = load(strcat(fileid,'-n=3-chi.txt'),'-ascii') ;
    A_saved = load(strcat(fileid,'-n=3-A.txt'),'-ascii') ;
%       test if the saved actual Pc, chi, A are equal  
%       (within numerical errors) to the "good" (correct) quantities
    testCase.verifyThat(Pc_saved, IsEqualTo(testCase.TestData.Pc, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(chi_saved, IsEqualTo(testCase.TestData.chi, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))    
    testCase.verifyThat(A_saved, IsEqualTo(testCase.TestData.A, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       check if val_vec and opt_vec match the trusted values
%       after the optimization loop
    testCase.verifyThat(val_vec, IsEqualTo(testCase.TestData.val_vec_nm(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(opt_vec, IsEqualTo(testCase.TestData.opt_vec_nm(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

function test_optimization_loop_serial_gaussNewton(testCase) %l2 + l4
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'gauss-newton' ;
    wk.id = 'test_optimization_loop_serial_gaussNewton' ;
%       initialize wk
    [ T, wk, fileid ] = evalc('initialize_workspace( wk )') ;
%       T  is stored to a file
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_serial_gaussNewton') ;
    fprintf(fileID,'%s\n\n\n',T) ;
%       verify that wk is fine
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_optimization_loop_serial_gaussNewton')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'gauss-newton')
    verifyEqual(testCase, wk.parallel, 0)
    verifyEqual(testCase, wk.maxiter, 100)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.xscale, 1e-6)
    verifyEqual(testCase, wk.xtol, 1e-4)
%       specify the test case data
    matrixfile = 'count.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
%       check if the sd calculated by get_knownInput() matches the trusted 
%       testCase.TestData.sd
    [ ~, ~, c ] = verifyAlmostEqual(sd, testCase.TestData.sd, abstol, ...
        reltol, false) ;
    verifyTrue(testCase, c)
%           calculate schur vectors
    [ EVS, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    kmin = 2 ;
    kmax = 3 ;
%       initialize A
    A_cell = cell([1,kmax]) ;
    for kt = kmin:kmax
        evs = EVS(:,1:kt) ;
        [ T, A_cell{kt} ] = evalc('initialize_A(evs, wk.init)') ;
        fprintf(fileID,'iteration %d\n',kt) ;
        fprintf(fileID,'%s\n\n\n',T) ;
    end
%       run optimization loop
    [ T, Pc, A_cell, ~, val_vec, opt_vec ] = evalc(['optimization_loop('...
    'double(P), double(sd), double(EVS), A_cell, kmin, kmax, wk, fileid )']) ;
%       T  is stored to a file
    fprintf(fileID,'%s\n\n\n',T) ;
    fclose(fileID) ;
%       ensure that A_cell is a cell
    verifyTrue(testCase,iscell(A_cell))
%       test if actual Pc, chi, A are equal (within numerical 
%       errors) to the "good" (correct) values 
    testCase.verifyThat(Pc, IsEqualTo(testCase.TestData.Pc_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))   
    testCase.verifyThat(A_cell{end}, IsEqualTo(testCase.TestData.A_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       load the saved actual Pc, chi, A from files
    fileid = ['test_optimization_loop_serial_gaussNewton-gauss-newton-'...
        'maxiter_100-xscale_1e-06-xtol_0.0001'] ;
    fileid = strcat(fileid,'-',testCase.TestData.precision) ;
    Pc_saved = load(strcat(fileid,'-n=3-Pc.txt'),'-ascii') ;
    A_saved = load(strcat(fileid,'-n=3-A.txt'),'-ascii') ;
%       test if the saved actual Pc, chi, A are equal  
%       (within numerical errors) to the "good" (correct) quantities
    testCase.verifyThat(Pc_saved, IsEqualTo(testCase.TestData.Pc_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))  
    testCase.verifyThat(A_saved, IsEqualTo(testCase.TestData.A_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       check if val_vec and opt_vec match the trusted values
%       after the optimization loop
    testCase.verifyThat(val_vec, IsEqualTo(testCase.TestData.val_vec(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(opt_vec, IsEqualTo(testCase.TestData.opt_vec(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

function test_optimization_loop_parallel_nelderMead(testCase) %l3 + l5 + l6 + l7
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'nelder-mead' ;
    wk.parallel = 1 ;
    wk.id = 'test_optimization_loop_parallel_nelderMead' ;
    wk.maxiter = 2000 ;
    wk.tolfun = 1e-08 ;
    wk.tolx = 1e-08 ;
%       initialize wk
    [ T, wk, fileid ] = evalc('initialize_workspace( wk )') ;
%       T  is stored to a file
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_parallel_nelderMead') ;
    fprintf(fileID,'%s\n\n\n',T) ;
%       verify that wk is fine
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_optimization_loop_parallel_nelderMead')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'nelder-mead')
    verifyEqual(testCase, wk.parallel, 1)
    verifyEqual(testCase, wk.maxiter, 2000)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.tolfun, 1e-8)
    verifyEqual(testCase, wk.tolx, 1e-8)
%       specify the test case data
    matrixfile = 'count.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
%       check if the sd calculated by get_knownInput() matches the trusted 
%       testCase.TestData.sd
    [ ~, ~, c ] = verifyAlmostEqual(sd, testCase.TestData.sd, abstol, ...
        reltol, false) ;
    verifyTrue(testCase, c)
%           calculate schur vectors
    [ EVS, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    kmin = 2 ;
    kmax = 3 ;
%       initialize A
    A_cell = cell([1,kmax]) ;
    for kt = kmin:kmax
        evs = EVS(:,1:kt) ;
        [ T, A_cell{kt} ] = evalc('initialize_A(evs, wk.init)') ;
        fprintf(fileID,'iteration %d\n',kt) ;
        fprintf(fileID,'%s\n\n\n',T) ;
    end
%       run optimization loop
    [ T, Pc, A_cell, chi, val_vec, opt_vec ] = evalc(['optimization_loop('...
    'double(P), double(sd), double(EVS), A_cell, kmin, kmax, wk, fileid )']) ;
%       T  is stored to a file
    fprintf(fileID,'%s\n\n\n',T) ;
    fclose(fileID) ;
%       ensure that A_cell is a cell
    verifyTrue(testCase,iscell(A_cell))
%       test if actual Pc, chi, A are equal (within numerical 
%       errors) to the "good" (correct) values 
    testCase.verifyThat(Pc, IsEqualTo(testCase.TestData.Pc, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(chi, IsEqualTo(testCase.TestData.chi, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))    
    testCase.verifyThat(A_cell{end}, IsEqualTo(testCase.TestData.A, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       load the saved actual Pc, chi, A from files
    fileid = ['test_optimization_loop_parallel_nelderMead-nelder-mead-'...
        'maxiter_2000-tolfun_1e-08-tolx_1e-08'] ;
    fileid = strcat(fileid,'-',testCase.TestData.precision) ;
    Pc_saved = load(strcat(fileid,'-n=3-Pc.txt'),'-ascii') ;
    chi_saved = load(strcat(fileid,'-n=3-chi.txt'),'-ascii') ;
    A_saved = load(strcat(fileid,'-n=3-A.txt'),'-ascii') ;
%       test if the saved actual Pc, chi, A are equal  
%       (within numerical errors) to the "good" (correct) quantities
    testCase.verifyThat(Pc_saved, IsEqualTo(testCase.TestData.Pc, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(chi_saved, IsEqualTo(testCase.TestData.chi, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))    
    testCase.verifyThat(A_saved, IsEqualTo(testCase.TestData.A, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       check if val_vec and opt_vec match the trusted values
%       after the optimization loop
    testCase.verifyThat(val_vec, IsEqualTo(testCase.TestData.val_vec_nm(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(opt_vec, IsEqualTo(testCase.TestData.opt_vec_nm(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

function test_optimization_loop_parallel_gaussNewton(testCase) %l3 + l5 + l6 + l7
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       set absolute and relative tolerance for verification
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       set the solver, optional solver and ID-string for file-naming
    wk.solver = 'gauss-newton' ;
    wk.parallel = 1 ;
    wk.id = 'test_optimization_loop_parallel_gaussNewton' ;
%       initialize wk
    [ T, wk, fileid ] = evalc('initialize_workspace( wk )') ;
%       T  is stored to a file
    fileID = fopen('test_optimization_loop_report.txt','a') ;
    fprintf(fileID,'%s\n\n','test_optimization_loop_parallel_gaussNewton') ;
    fprintf(fileID,'%s\n\n\n',T) ;
%       verify that wk is fine
    verifyEqual(testCase, isstruct(wk), true)
    verifyEqual(testCase, wk.id, 'test_optimization_loop_parallel_gaussNewton')
    verifyEqual(testCase, wk.schur, 1)
    verifyEqual(testCase, wk.b, 0)
    verifyEqual(testCase, wk.init, 1)
    verifyEqual(testCase, wk.solver, 'gauss-newton')
    verifyEqual(testCase, wk.parallel, 1)
    verifyEqual(testCase, wk.maxiter, 100)
    verifyEqual(testCase, wk.display, 0)
    verifyEqual(testCase, wk.xscale, 1e-6)
    verifyEqual(testCase, wk.xtol, 1e-4)
%       specify the test case data
    matrixfile = 'count.txt' ;
%       get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
%       check if the sd calculated by get_knownInput() matches the trusted 
%       testCase.TestData.sd
    [ ~, ~, c ] = verifyAlmostEqual(sd, testCase.TestData.sd, abstol, ...
        reltol, false) ;
    verifyTrue(testCase, c)
%           calculate schur vectors
    [ EVS, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    kmin = 2 ;
    kmax = 3 ;
%       initialize A
    A_cell = cell([1,kmax]) ;
    for kt = kmin:kmax
        evs = EVS(:,1:kt) ;
        [ T, A_cell{kt} ] = evalc('initialize_A(evs, wk.init)') ;
        fprintf(fileID,'iteration %d\n',kt) ;
        fprintf(fileID,'%s\n\n\n',T) ;
    end
%       run optimization loop
    [ T, Pc, A_cell, ~, val_vec, opt_vec ] = evalc(['optimization_loop('...
    'double(P), double(sd), double(EVS), A_cell, kmin, kmax, wk, fileid )']) ;
%       T  is stored to a file
    fprintf(fileID,'%s\n\n\n',T) ;
    fclose(fileID) ;
%       ensure that A_cell is a cell
    verifyTrue(testCase,iscell(A_cell))
%       test if actual Pc, chi, A are equal (within numerical 
%       errors) to the "good" (correct) values 
    testCase.verifyThat(Pc, IsEqualTo(testCase.TestData.Pc_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))   
    testCase.verifyThat(A_cell{end}, IsEqualTo(testCase.TestData.A_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       load the saved actual Pc, chi, A from files
    fileid = ['test_optimization_loop_serial_gaussNewton-gauss-newton-'...
        'maxiter_100-xscale_1e-06-xtol_0.0001'] ;
    fileid = strcat(fileid,'-',testCase.TestData.precision) ;
    Pc_saved = load(strcat(fileid,'-n=3-Pc.txt'),'-ascii') ;
    A_saved = load(strcat(fileid,'-n=3-A.txt'),'-ascii') ;
%       test if the saved actual Pc, chi, A are equal  
%       (within numerical errors) to the "good" (correct) quantities
    testCase.verifyThat(Pc_saved, IsEqualTo(testCase.TestData.Pc_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))  
    testCase.verifyThat(A_saved, IsEqualTo(testCase.TestData.A_gn, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
%       check if val_vec and opt_vec match the trusted values
%       after the optimization loop
    testCase.verifyThat(val_vec, IsEqualTo(testCase.TestData.val_vec(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    testCase.verifyThat(opt_vec, IsEqualTo(testCase.TestData.opt_vec(1:2,1:2), ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
end

%function test_optimization_loop_parallel_catchME(testCase) %l8
%
%end

% -------------------------------------------------------------------------
% test postprocess (logic paths: 6) (assertions: 13)

% so far NO ASSERTION TESTS!

function test_postprocess(testCase) % l1
%       test error for the case that Q, R arent of the same class
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
%       check for count matrix and n_clusters=3
    matrixfile = 'count.txt' ;
%           get the known test case data
    [ P, sd1 ] = get_knownInput( matrixfile ) ;
    k = 3 ;
    val = 1.78e-01 ;
    chi = testCase.TestData.chi ;
    %pi = testCase.TestData.pi ;
    sd = testCase.TestData.sd ;
%       check if sd (and therefore likely also P) was correctly calculated
%       by get_knownInput()
    abstol = testCase.TestData.abstolfac * eps(numeric_t) ;
    reltol = testCase.TestData.reltolfac * eps(numeric_t) ;
    [ ~, ~, c ] = verifyAlmostEqual(sd1, sd, abstol, reltol, true) ;
    verifyTrue(testCase, c)
%       convert sd, pi and P to double (in case testing is in multiprecision 
%       mode), since postprocess should work with doubles
    sd = double(sd) ;
    P = double(P) ;
    %pi = double(pi) ;
    A = testCase.TestData.A ;
    fileid = 'test_postprocess' ;
    %[ T, Pc ] = evalc('postprocess( k, val, chi, pi, sd, A, P, fileid )') ;
    [ T, Pc ] = evalc('postprocess( k, val, chi, sd, A, P, fileid )') ;
    fileID = fopen('test_postprocess_report.txt','w') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
    testCase.verifyThat(Pc, IsEqualTo(testCase.TestData.Pc, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    name = strcat(fileid,'-n=3-chi.txt') ;
    chi_act = load(name,'-ascii') ;
    dummy_act = chi_act' * chi_act ;
    dummy_exp = testCase.TestData.chi' * testCase.TestData.chi ;
    [ ~, ~, c1 ] = verifyAlmostEqual(dummy_act, dummy_exp, ...
    abstol, reltol, false) ;
    [ ~, ~, c2 ] = verifyAlmostEqual(chi_act, testCase.TestData.chi, ...
    abstol, reltol, false) ;
    verifyTrue(testCase, (c1 | c2))
    %name = strcat(fileid,'-n=3-pi_weights.txt') ;
    %pi_weights_act = load(name,'-ascii') ;
    %testCase.verifyThat(pi_weights_act, IsEqualTo(testCase.TestData.pi_weights, ...
       %'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    name = strcat(fileid,'-n=3-sd_weights.txt') ;
    sd_weights_act = load(name,'-ascii') ;
    testCase.verifyThat(sd_weights_act, IsEqualTo(testCase.TestData.sd_weights, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    name = strcat(fileid,'-n=3-A.txt') ;
    A_act = load(name,'-ascii') ;
    testCase.verifyThat(A_act, IsEqualTo(testCase.TestData.A, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    name = strcat(fileid,'-n=3-A.mat') ;
    load(name,'A') ;
    testCase.verifyThat(A, IsEqualTo(testCase.TestData.A, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    name = strcat(fileid,'-n=3-Pc.txt') ;
    Pc_act = load(name,'-ascii') ;
    testCase.verifyThat(Pc_act, IsEqualTo(testCase.TestData.Pc, ...
       'Within', AbsoluteTolerance(abstol) | RelativeTolerance(reltol)))
    name = strcat(fileid,'-n=3-idx_vec.txt') ;
    idx_vec_act = load(name,'-ascii') ;
    testCase.verifyThat(idx_vec_act, IsEqualTo(testCase.TestData.idx_vec, ...
       'Within', AbsoluteTolerance(eps) & RelativeTolerance(eps)))
    name = strcat(fileid,'-n=3-chi-figure.pdf') ;
    verifyTrue( testCase, exist(name, 'file')==2 )
    name = strcat(fileid,'-n=3-chi-figure.fig') ;
    verifyTrue( testCase, exist(name, 'file')==2 )
    name = strcat(fileid,'-n=3-chi-ordered-figure.pdf') ;
    verifyTrue( testCase, exist(name, 'file')==2 )
    name = strcat(fileid,'-n=3-chi-ordered-figure.fig') ;
    verifyTrue( testCase, exist(name, 'file')==2 )
end

% -------------------------------------------------------------------------
% test save_t (depends on none) (logic paths: 2) (assertions: 0)

function test_save_t_double(testCase) % l1
%       test fouble case
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'double')
    dummy = double(1.0)/double(3.0) ;
    name = 'test_save_t_double.txt' ;
    save_t( name, dummy, '-ascii' )
    dummy1 = load(name, '-ascii') ;
    [ ~, ~, c ] = nearly_equal(dummy, dummy1) ;
    verifyTrue(testCase, c)
end

function test_save_t_mp(testCase) % l2
%       test mp case
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
    dummy = mp('1/3') ;
    name = 'test_save_t_mp.txt' ;
    save_t( name, dummy, '-ascii' )
    dummy1 = mp.read(name) ;
    [ ~, ~, c ] = nearly_equal(dummy, dummy1) ;
    verifyTrue(testCase, c)
end
% -------------------------------------------------------------------------
% test SRSchur_num_t (depends on numeric_t)

function test_SRSchur_num_t_MatrixShapeError1(testCase) % a1
%       test assertion for the case that Q isn't quadratic
    Q = zeros(4) ;
    Q(:,4) = [] ;
    R = zeros(4) ;
    try
        [ ~, ~, ~ ] = SRSchur_num_t( Q, R, 1, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'SRSchur_num_t:MatrixShapeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_SRSchur_num_t_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_SRSchur_num_t_MatrixShapeError2(testCase) % a2
%       test assertion for the case that R isn't quadratic
    Q = zeros(4) ;
    R = zeros(4) ;
    R(:,4) = [] ;
    try
        [ ~, ~, ~ ] = SRSchur_num_t( Q, R, 1, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'SRSchur_num_t:MatrixShapeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_SRSchur_num_t_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_SRSchur_num_t_MatchError(testCase) % a3
%       test assertion if dimensions of R doesnt match to those of Q
    Q = zeros(4) ;
    R = zeros(3) ;
    try
        [ ~, ~, ~ ] = SRSchur_num_t( Q, R, 1, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'SRSchur_num_t:MatchError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_SRSchur_num_t_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_SRSchur_num_t_MatrixShapeError3(testCase) % a4
%       test assertion for the case that Q, R arent at least (2x2)-matrices
    Q = 0 ;
    R = 0 ;
    try
        [ ~, ~, ~ ] = SRSchur_num_t( Q, R, 1, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'SRSchur_num_t:MatrixShapeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_SRSchur_num_t_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_SRSchur_num_t_DataTypeError(testCase) % l1
%       test error for the case that Q, R arent of the same class
    Q = zeros(4) ;
    R = single(zeros(4)) ;
    try
        [ ~, ~, ~ ] = SRSchur_num_t( Q, R, 1, 0 ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'SRSchur_num_t:DataTypeError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_SRSchur_num_t_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_SRSchur_num_t(testCase) % l
% test normal execution with test cases 1-4 from Table I in
% (1) Brandts, J. H. Matlab Code for Sorting Real Schur Forms. Numer. 
% Linear Algebr. with Appl. 2002, 9 (3), 249-261
% and a fifth test case designed from the first test matrix M1 by swapping 
% of the first block A11 with the second block A22.
% The formula for EA had a typo in (1), the correct formula was taken from
% (2) Bai, Z.; Demmel, J. W. On Swapping Diagonal Blocks in Real Schur
% Form. Linear Algebra Appl. 1993, 186 (C), 75-95.
% Results differ from Brandts et al. and Bai et al. since the sorting
% criterium in select() was changed.
    m1 = numeric_t('[ 2   -87 -20000   10000 ]') ;
    m2 = numeric_t('[ 5    2  -20000  -10000 ]') ;
    m3 = numeric_t('[ 0    0   1      -11 ]') ;
    m4 = numeric_t('[ 0    0   37      1 ]') ;
    M1 = [ m1; m2; m3; m4 ] ;
    m1 = numeric_t('[ 1   -3   3576    4888 ]') ;
    m2 = numeric_t('[ 1    1  -88     -1440 ]') ;
    m3 = numeric_t('[ 0    0   1.001  -3 ]') ;
    m4 = numeric_t('[ 0    0   1.001   1.001 ]') ;
    M2 = [ m1; m2; m3; m4 ] ;
    m1 = numeric_t('[ 1        -100     400    -1000 ]') ;
    m2 = numeric_t('[ 0.01      1       1200   -10 ]') ;
    m3 = numeric_t('[ 0         0       1.001  -0.01 ]') ;
    m4 = numeric_t('[ 0         0       100     1.001 ]') ;
    M3 = [ m1; m2; m3; m4 ] ;
    m1 = numeric_t('[ 1     -1e4    8812        4566 ]') ;
    m2 = numeric_t('[ 1e-4   1     -9           1200 ]') ;
    m3 = numeric_t('[ 0      0      1e-5+1     -1e-4 ]') ;
    m4 = numeric_t('[ 0      0      1e4         1e-5+1 ]') ;
    M4 = [ m1; m2; m3; m4 ] ;
    m1 = numeric_t('[ 1   -11   -20000   10000 ]') ;
    m2 = numeric_t('[ 37   1    -20000  -10000 ]') ;
    m3 = numeric_t('[ 0    0     2      -87 ]') ;
    m4 = numeric_t('[ 0    0     5      2 ]') ;
    M5 = [ m1; m2; m3; m4 ] ;
    M = { M1, M2, M3, M4, M5 } ;
    Q = numeric_t( eye(4) ) ;
    for i = 1:5
        R = M{i} ;
        [ T, QQ, RR, ap ] = evalc('SRSchur_num_t( Q, R, 1, 0 )') ;
        fileID = fopen('test_SRSchur_num_t_report.txt','a') ;
        fprintf(fileID,'%s\n','test_SRSchur_num_t') ;
        fprintf(fileID,'%s\n\n',T) ;
        fclose(fileID) ;
        if ~isempty(ap)
            verifyLessThanOrEqual(testCase, double(ap), 1)
        end
        EQ = norm( ( numeric_t(eye(4)) - QQ'*QQ ), 1 ) / eps(numeric_t) ;
        verifyEqual(testCase, double(EQ), 1.0, 'AbsTol', 5.0)
        EA = norm( R - QQ*RR*QQ', 1 ) / ( eps(numeric_t)*norm(R,1) ) ;
        verifyEqual(testCase, double(EA), 1.0, 'AbsTol', 5.0)
        opts.disp = 0 ;
        l = eigs( R, 4, 'sm', opts) ;
        ll = eigs( RR, 4, 'sm', opts) ;
        El1 = abs( l(1) - ll(1) ) / ( eps(numeric_t)*abs(l(1)) ) ;
        verifyEqual(testCase, double(El1), 1.0, 'AbsTol', 5.0)
        El2 = abs( l(2) - ll(2) ) / ( eps(numeric_t)*abs(l(2)) ) ;
        verifyEqual(testCase, double(El2), 1.0, 'AbsTol', 5.0)
        El3 = abs( l(3) - ll(3) ) / ( eps(numeric_t)*abs(l(3)) ) ;
        verifyEqual(testCase, double(El3), 1.0, 'AbsTol', 5.0)
        El4 = abs( l(4) - ll(4) ) / ( eps(numeric_t)*abs(l(4)) ) ;
        verifyEqual(testCase, double(El4), 1.0, 'AbsTol', 5.0)
    end
end

function test_SRSchur_num_t_partial_ordering(testCase) % l
    fileid = 'test_SRSchur_num_t_partial_ordering' ;
    matrixfile = 'count.txt' ;
%           get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
%       weight the stochastic matrix P_stoch by the initial distribution sd
    Pd=diag(sqrt(sd))*P*diag(numeric_t('1.0')./sqrt(sd)) ;
%       make a Schur decomposition of Pd
    [Q, R]=schur(Pd) ;
    name = strcat(fileid,'-Q.txt') ;
    save_t(name,Q,'-ascii')
    name = strcat(fileid,'-R.txt') ;
    save_t(name,R,'-ascii')
%       sort the Schur matrix and vectors
    target=numeric_t('1.0') ;
%       calculate the case of full sorting
    [QQ, RR, ap]=SRSchur_num_t(Q,R,target,0) ;
    name = strcat(fileid,'-QQ.txt') ;
    save_t(name,QQ,'-ascii')
    name = strcat(fileid,'-RR.txt') ;
    save_t(name,RR,'-ascii')
%       test, if the indicator vector ap is <1 everywere
    dummy = ap < numeric_t('1.0') ;
    verifyTrue(testCase,all(dummy(:)))
%       test, if Pd*QQ=QQ*RR (Schur decomposition)
    dummy = ( abs(QQ*RR - Pd*QQ) < (testCase.TestData.abstolfac * eps(numeric_t)) ) ;
    verifyTrue(testCase,all(dummy(:)))
%       calculate the case of partial sorting (first 10 eigenvalues)
    [QQ1, RR1, ap]=SRSchur_num_t(Q,R,target,10) ;
    name = strcat(fileid,'-QQ1.txt') ;
    save_t(name,QQ1,'-ascii')
    name = strcat(fileid,'-RR1.txt') ;
    save_t(name,RR1,'-ascii')
%       test, if the indicator vector ap is <1 everywere
    dummy = ap < numeric_t('1.0') ;
    verifyTrue(testCase,all(dummy(:)))
%       test, if Pd*QQ=QQ*RR (Schur decomposition)
    dummy = ( abs(QQ*RR - Pd*QQ) < (testCase.TestData.abstolfac * eps(numeric_t)) ) ;
    verifyTrue(testCase,all(dummy(:)))
%       compare the first ten rows/columns of the fully sorted Schur matrix 
%       with the fist ten rows/columns of the partially sorted one
    abstol = testCase.TestData.abstolfac * eps(numeric_t) ;
    reltol = testCase.TestData.reltolfac * eps(numeric_t) ;
    actVal = RR1(1:10,1:10) ;
    expVal = RR(1:10,1:10) ;
    [ ~, ~, c ] = verifyAlmostEqual(actVal, expVal, ...
    abstol, reltol, false) ;
    verifyTrue(testCase,c)
%       compare the first ten columns of the fully sorted Schur-vector
%       matrix with the fist ten colums of the partially sorted one
    actVal = QQ1(:,1:10) ;
    expVal = QQ(:,1:10) ;
    [ ~, ~, c ] = verifyAlmostEqual(actVal, expVal, ...
    abstol, reltol, false) ;
    verifyTrue(testCase,c)
end

% -------------------------------------------------------------------------
% test cluster_by_isa (depends on do_schur) (logic paths: 10) (assertions: 1)
% and use_minChi (depends on do_schur, cluster_by_isa) (logic paths: 3) 
% (assertions: 5)

function test_cluster_by_isa_FirstColumnError(testCase) % a1
%       test for the assertion if first column of schur vector matrix isnt
%       constantly equal 1
    d1 = [0,0,0] ;
    d2 = [0,0,0] ;
    d3 = [0,0,0] ;
    d4 = [0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    try
        [ ~, ~ ] = cluster_by_isa( svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'cluster_by_isa:FirstColumnError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_cluster_by_isa_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_cluster_by_isa_MatrixShapeError3(testCase) % l4
%       test assertion for schur vector (N,k)-matrix  with k>N
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    svecs = [ d1; d2; d3 ] ;
    try
        [ ~, ~ ] = cluster_by_isa( svecs ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'cluster_by_isa:MatrixShapeError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_cluster_by_isa_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_cluster_by_isa(testCase) % l
    import matlab.unittest.constraints.IsEqualTo;
    import matlab.unittest.constraints.AbsoluteTolerance;
    import matlab.unittest.constraints.RelativeTolerance;
    abstol = testCase.TestData.abstolfac * eps ;
    reltol = testCase.TestData.reltolfac * eps ;
%       check for example-matrix with mu=0 and n_clusters=3
    matrixfile = 'example_matrix_mu0.txt' ;
%           get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_cluster_by_isa_example_matrix_mu0' ;
    [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    [ T, chi, ~ ] = evalc('cluster_by_isa( X(:,1:3) )') ;
    fileID = fopen('test_cluster_by_isa_report.txt','a') ;
    fprintf(fileID,'%s\n','test_cluster_by_isa') ;
    fprintf(fileID,'%s\n\n',T) ;
    chi_exp = testCase.TestData.chi_isa_mu0_n3 ;
    dummy_act = chi'*chi ;
    dummy_exp = chi_exp'*chi_exp ;
    [ ~, ~, c1 ] = verifyAlmostEqual(dummy_act, dummy_exp, ...
    abstol, reltol, false) ;
    [ ~, ~, c2 ] = verifyAlmostEqual(chi, chi_exp, ...
    abstol, reltol, false) ;
    verifyTrue(testCase, (c1 | c2))
%       check for example-matrix with mu=100 and n_clusters=3
    matrixfile = 'example_matrix_mu100.txt' ;
%           get the known test case data
    [ P, sd ] = get_knownInput( matrixfile ) ;
    fileid = 'test_cluster_by_isa_example_matrix_mu100' ;
    [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
    [ T, chi, ~ ] = evalc('cluster_by_isa( X(:,1:3) )') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    chi_exp = testCase.TestData.chi_isa_mu100_n3 ;
    dummy_act = chi'*chi ;
    dummy_exp = chi_exp'*chi_exp ;
    [ ~, ~, c1 ] = verifyAlmostEqual(dummy_act, dummy_exp, ...
    abstol, reltol, false) ;
    [ ~, ~, c2 ] = verifyAlmostEqual(chi, chi_exp, ...
    abstol, reltol, false) ;
    verifyTrue(testCase, (c1 | c2))
    
end

function test_use_minChi_EmptyInput1(testCase) % a1
%       test assertion for empty kmin
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    d4 = [1,0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    fileid = '' ;
    try
        [ ~ ] = use_minChi( svecs, [], 4, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'use_minChi:EmptyInput1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_use_minChi_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_use_minChi_EmptyInput2(testCase) % a2
%       test assertion for empty kmax
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    d4 = [1,0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    fileid = '' ;
    try
        [ ~ ] = use_minChi( svecs, 1, [], fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'use_minChi:EmptyInput2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_use_minChi_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_use_minChi_InputError1(testCase) % a3
%       test assertion for noninteger kmin
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    d4 = [1,0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    fileid = '' ;
    try
        [ ~ ] = use_minChi( svecs, 1.1, 4, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'use_minChi:InputError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_use_minChi_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_use_minChi_InputError2(testCase) % a4
%       test assertion for noninteger kmax
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    d4 = [1,0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    fileid = '' ;
    try
        [ ~ ] = use_minChi( svecs, 1, 4.2, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'use_minChi:InputError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_use_minChi_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_use_minChi_InputError3(testCase) % a5
%       test assertion for kmin > kmax
    d1 = [1,0,0,0] ;
    d2 = [1,0,0,0] ;
    d3 = [1,0,0,0] ;
    d4 = [1,0,0,0] ;
    svecs = [ d1; d2; d3; d4 ] ;
    fileid = '' ;
    try
        [ ~ ] = use_minChi( svecs, 3, 2, fileid ) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'use_minChi:InputError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_use_minChi_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_use_minChi(testCase)
% The results to compare with are taken from page 168 of
% (1) R?blitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+: 
% Application to Markov State Models and Data Classification. Adv. Data
% Anal. Classif. 2013, 7 (2), 147?179.
% There is a difference: In (1) the optimal cluster numbers for 
% mu=0,10,50,100,200,500,1000 are kopt=3,3,3,3,2,2,5 while here the
% correct optimal cluster numbers are assumed to be kopt=3,3,3,3,3,3,7
% since in all tests we obtained for mu=1000 an optimal cluster number of
% kopt=7 and k=2 wasn't included in the search for the maximum min(chi),
% since min(chi) for k=2 is always (nearly) equal 0 in case of the Inner
% Simplex Algorithm.
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
    kmin = 1 ;
    kmax = 9 ;
    fileID = fopen('test_use_minChi_report.txt','a') ;
    fprintf(fileID,'%s\n','test_use_min_chi') ;
    i = 0 ;
    kopt_vec = [] ;
    for mu = [ 0 10 50 100 200 500 1000 ]
        i = i + 1 ;
        matrixfilename = strcat('example_matrix_mu',num2str(mu)) ;
        matrixfile = strcat(matrixfilename,'.txt') ;
%           get the known test case data
        [ P, sd ] = get_knownInput( matrixfile ) ;
        fileid = strcat('test_use_minChi_',matrixfilename) ;
        [ X, ~ ]  = do_schur( P, sd, fileid, 0 ) ;
        [ T, minChi_vec ] = evalc('use_minChi( X, kmin, kmax, fileid )') ;
%           close all open (but invisible) plots
        close all
%           flip minChi_vec, since we want the biggest k with maximum
%           min(chi)
        [ ~ , dummy ] = max(flipud(minChi_vec(3:8,2))) ;
%           since kmax=9 but only flipud(minChi_vec(3:8,2)) is used, we
%           have kopt=kmax-dummy instead of kopt=kmax-dummy+1
        kopt = kmax - dummy ;
        kopt_vec(i) = kopt ;
%           T, T1 (with all screen output of gpcca()) is stored to a file
        fprintf(fileID,'%s\n\n',T) ;
    end 
    fclose(fileID) ;
    verifyEqual(testCase,kopt_vec,[3,3,3,3,3,3,7])
end
%--------------------------------------------------------------------------
% test verifyAlmostEqual (depends on numeric_t) (logic paths: 18) 
% (assertions: 4)

function test_verifyAlmostEqual_DataTypeError(testCase) % a1 + a2
%       test the case that a isnt double or mp
    a = single(2.0) ;
    b = numeric_t('2.0') ;
    try
        [ ~, ~, ~ ] = verifyAlmostEqual(a, b, eps(numeric_t), 0) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'verifyAlmostEqual:DataTypeError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_verifyAlmostEqual_report.txt','w') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
%       test the case that b isnt double or mp    
    a = numeric_t('2.0') ;
    b = single(2.0) ;
    try
        [ ~, ~, ~ ] = verifyAlmostEqual(a, b, eps(numeric_t), 0) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'verifyAlmostEqual:DataTypeError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_verifyAlmostEqual_DataTypeError3(testCase) % a3
    a = numeric_t('2.0') ;
    b = numeric_t('2.0') ;
    try
        [ ~, ~, ~ ] = verifyAlmostEqual(a, b, eps, eps, 0) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'verifyAlmostEqual:DataTypeError3'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_verifyAlmostEqual_mixed(testCase) % l1+l14
    prec = testCase.TestData.precision ;
    assumeEqual(testCase, prec, 'mp')
    a = ( numeric_t('2.0') + eps ) ;
    b = double(2.0) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 2*eps, 0)') ;
 %       T is stored to a file
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_mixed') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, 'double'))
    verifyTrue(testCase, isa(b, 'double'))
end

function test_verifyAlmostEqual_FatalError2(testCase) % l2
    a = numeric_t('1.0') ;
    b = numeric_t('2.0') ;
    try
        [ ~, ~, ~ ] = verifyAlmostEqual(a, b, 0, 0) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'verifyAlmostEqual:FatalError2'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_verifyAlmostEqual_InputError(testCase) % l4
    a = numeric_t('2.0') ;
    b = numeric_t('2.0') ;
    try
        [ ~, ~, ~ ] = verifyAlmostEqual(a, b, eps, eps, true, 0) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'verifyAlmostEqual:InputError'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_verifyAlmostEqual_exact1(testCase) % l5
    a = numeric_t('2')*numeric_t(ones(3)) ;
    b = numeric_t('2')*numeric_t(ones(3)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 0, 0) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_exact2(testCase) % l5
    a = numeric_t('2')*numeric_t(ones(3)) ;
    b = numeric_t('2')*numeric_t(ones(3)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, eps(numeric_t), 0) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_exact3(testCase) % l5
    a = numeric_t('2')*numeric_t(ones(3)) ;
    b = numeric_t('2')*numeric_t(ones(3)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 0, eps(numeric_t)) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_inf(testCase) % l5
    a = inf ;
    b = inf ;
    [ a, b, c ] = verifyAlmostEqual(a, b, eps(numeric_t), 0) ;
    verifyTrue(testCase, c)
    verifyEqual(testCase, a, inf)
    verifyEqual(testCase, b, inf)
end

function test_verifyAlmostEqual_zero1(testCase) % l5
    a = numeric_t('0.0') ;
    b = numeric_t('0.0') ;
    [ a, b, c ] = verifyAlmostEqual(a, b, eps(numeric_t), 0) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_zero2(testCase) % l5
    a = numeric_t('0.0') ;
    b = numeric_t('0.0') ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 0, eps(numeric_t)) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_FatalError1(testCase) % l6
    a = numeric_t('1.0') ;
    b = numeric_t('2.0') ;
    try
        [ ~, ~, ~ ] = verifyAlmostEqual(a, b, 0, 0, true) ;
        actVal = 0 ;
    catch ME
        switch ME.identifier
            case 'verifyAlmostEqual:FatalError1'
                actVal = 1 ;
            otherwise
                rethrow(ME)
        end
    end
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s',ME.identifier) ;
    fprintf(fileID,'%s',': ') ;
    fprintf(fileID,'%s\n\n',ME.message) ;
    fclose(fileID) ;
    expVal = 1 ;
    verifyEqual(testCase, actVal, expVal)
end

function test_verifyAlmostEqual_abs_and_rel1(testCase) % l3+l7+l8+l9
% both ok 
    a = numeric_t(ones(3)) + eps(numeric_t) ;
    b = numeric_t(ones(3)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), true) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_and_rel2(testCase) % l3+l7+l8+l9+l10
% reltol unsatisfied
    a = eps(numeric_t) + eps(numeric_t) ;
    b = eps(numeric_t) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), true)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_abs_and_rel2') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_and_rel3(testCase)  % l3+l7+l8+l9+l10
% abstol unsatisfied 
    a = numeric_t('1e6') + numeric_t('1e6')*eps(numeric_t) ;
    b = numeric_t('1e6') ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), true)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_abs_and_rel3') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_and_rel4(testCase) % l3+l7+l8+l9+l10
% abstol and reltol unsatisfied, partially absA==absA==0 case 
    a = numeric_t(eye(4)) + 2*eps(numeric_t)*numeric_t(eye(4)) ;
    b = numeric_t(eye(4)) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, eps(numeric_t), eps(numeric_t), true)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_abs_and_rel4') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_or_rel1(testCase) % l3+l11+l12+l13
% both ok 
    a = numeric_t(ones(3)) + eps(numeric_t) ;
    b = numeric_t(ones(3)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), false) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_or_rel2(testCase)  % l3+l11+l12+l13
% abstol ok 
    a = eps(numeric_t) + eps(numeric_t) ;
    b = eps(numeric_t) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), false) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_or_rel3(testCase)  % l3+l11+l12+l13
% reltol ok 
    a = numeric_t('1e6') + numeric_t('1e6')*eps(numeric_t) ;
    b = numeric_t('1e6') ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), false) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end
 
function test_verifyAlmostEqual_abs_or_rel4(testCase) % l3+l11+l12+l13+l14
% both wrong
    a = eps(numeric_t) + 3*eps(numeric_t) ;
    b = eps(numeric_t) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), false)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_abs_or_rel4') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs_or_rel5(testCase) % l3+l11+l12+l13+l14
% both wrong, partially absA==absA==0 case
    a = numeric_t(eye(4))*eps(numeric_t) + numeric_t(eye(4))*3*eps(numeric_t) ;
    b = numeric_t(eye(4))*eps(numeric_t) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 2*eps(numeric_t), 2*eps(numeric_t), false)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_abs_or_rel5') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_abs(testCase) % l15+l16
    a = numeric_t(eye(4)) + eps(numeric_t) ;
    b = numeric_t(eye(4)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 3*eps(numeric_t), 0) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t(ones(4)) + 2*eps(numeric_t) ;
    b = numeric_t(ones(4)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 3*eps(numeric_t), 0) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('1.0') + 3*eps(numeric_t) ;
    b = numeric_t('1.0') ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 3*eps(numeric_t), 0)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_abs') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('2.0') + 4*eps(numeric_t) ;
    b = numeric_t('2.0') ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 5*eps(numeric_t), 0) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('2')*numeric_t(ones(3)) + 6*eps(numeric_t) ;
    b = numeric_t('2')*numeric_t(ones(3)) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 5*eps(numeric_t), 0)') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end

function test_verifyAlmostEqual_rel(testCase) % l17+l18
    a = numeric_t(eye(4)) + eps(numeric_t) ;
    b = numeric_t(eye(4)) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 0, 2*eps(numeric_t))') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n','test_verifyAlmostEqual_rel') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t(ones(4)) + 2*eps(numeric_t) ;
    b = numeric_t(ones(4)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 0, 3*eps(numeric_t)) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('1.0') + 5*eps(numeric_t) ;
    b = numeric_t('1.0') ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 0, 2*eps(numeric_t))') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('2.0') + 4*eps(numeric_t) ;
    b = numeric_t('2.0') ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 0, 3*eps(numeric_t)) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
    a = numeric_t('2')*numeric_t(ones(3)) + 10*eps(numeric_t) ;
    b = numeric_t('2')*numeric_t(ones(3)) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 0, 2*eps(numeric_t))') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
%       partially absA==absA==0 case
    a = numeric_t('2')*numeric_t(eye(3)) + eps(numeric_t)*numeric_t(eye(3)) ;
    b = numeric_t('2')*numeric_t(eye(3)) ;
    [ a, b, c ] = verifyAlmostEqual(a, b, 0, 2*eps(numeric_t)) ;
    verifyTrue(testCase, c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
    
%       partially absA==absA==0 case
    a = numeric_t('2')*numeric_t(eye(3)) + 10*eps(numeric_t)*numeric_t(eye(3)) ;
    b = numeric_t('2')*numeric_t(eye(3)) ;
    [ T, a, b, c ] = evalc('verifyAlmostEqual(a, b, 0, 2*eps(numeric_t))') ;
    fileID = fopen('test_verifyAlmostEqual_report.txt','a') ;
    fprintf(fileID,'%s\n\n',T) ;
    fclose(fileID) ;
    verifyTrue(testCase, ~c)
    verifyTrue(testCase, isa(a, numeric_t))
    verifyTrue(testCase, isa(b, numeric_t))
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function setupOnce(testCase)
% test environment setup
%   Written by Bernhard Reuter, Theoretical Physics II,
%   University of Kassel, 2017 
%       set global variable to make all plots invisible but still store
%       them
    global figureVisible
    figureVisible = false ;
%       request precision from user
    disp('You can choose in which precision the test is performed: ')
    disp('Either double precision (16 significant digits, Matlab ')
    disp('standard) or multi-precision with 50 significant digits. ')
    disp('The multi-precision mode is only applicable, if you have ')
    disp('the Matlab multiprecision toolbox from Advanpix installed!')
    testCase.TestData.precision = input(['Enter the precision to use ', ...
    'during the test IN QUOTES (mp or double): ']) ;
    prec = testCase.TestData.precision ;
    assert(ischar(prec),'setupOnce:InputError', ...
        'The precision you passed is not a string!')
    if ~strcmpi(prec,'mp') && ~strcmpi(prec,'double')
        error('setupOnce:InputError','The precision you passed is invalid!')
    end
%       create and change to temporary folder
    testCase.TestData.origPath = pwd ;
    testCase.TestData.tmpFolder = strcat('unitTestsFolder_',prec) ;
    tmp = testCase.TestData.tmpFolder ;
    if exist(tmp, 'dir') == 7
        rmdir(tmp, 's')
    end
    mkdir(testCase.TestData.tmpFolder)
    copyfile('cluster_by_isa.m', tmp)
    copyfile('count.txt', tmp)
    copyfile('do_schur.m', tmp)
    copyfile('fillA.m', tmp)
    copyfile('find_twoblocks.m', tmp)
    copyfile('getopt.m', tmp)
    copyfile('gram_schmidt_mod.m', tmp)
    copyfile('indexsearch.m', tmp)
    copyfile('initialize_A.m', tmp)
    copyfile('initialize_optional.m', tmp)
    copyfile('initialize_workspace.m', tmp)
    copyfile('inputT.m', tmp)
    copyfile('load_t.m', tmp)
    copyfile('main_nlscon.m', tmp)
    copyfile('nearly_equal.m', tmp)
    copyfile('nlscon.m', tmp)
    copyfile('numeric_t.m', tmp)
    copyfile('objective.m', tmp)
    copyfile('opt_soft.m', tmp)
    copyfile('optimization_loop.m', tmp)
    copyfile('gpcca.m', tmp)
    copyfile('plotmatrix.m',tmp)
    copyfile('postprocess.m', tmp)
    copyfile('problem_pcca_nlscon.m', tmp)
    copyfile('save_t.m', tmp)
    copyfile('SRSchur_num_t.m', tmp)
    copyfile('use_minChi.m', tmp)
    cd(testCase.TestData.tmpFolder)
%       factor abstolfac to multiply eps with for use as absolute tolerance 
%       abstol=abstolfac*eps in many (not all) comparisons 
    testCase.TestData.abstolfac = 1e8 ;
%       factor reltolfac to multiply eps with for use as relative tolerance 
%       reltol=reltolfac*eps in many (not all) comparisons 
    testCase.TestData.reltolfac = 1e8 ;
%       global variable to define precision to use
    global class_t ;
%       set precision of variables or expressions wrapped by numeric_t(),
%       i.e. 'double' or 'mp' for multipresicion (with multiprecision toolbox)
    if strcmpi(testCase.TestData.precision,'mp')
        mp.Digits(50) ;
        class_t = 'mp' ;
    elseif strcmpi(testCase.TestData.precision,'double')
        class_t = 'double' ;
    else
        error('knownInput:PrecisionError','Unknown Precision')
    end
%       create .mat files used in testmode from inputT() to get the input
%       values normaly passed interactively from keyboard
    kmin = 2 ;
    kmax = 5 ;
    koptmin = 3 ;
    koptmax = 3 ;
    Switch = 1 ;
    minChi_switch = 0 ;
    oldfileid = 'test_initialize_A_testmatrix' ;
    name = 'testInput_minChiOFF.mat' ;
    save(name,'kmin','kmax','koptmin','koptmax','Switch','minChi_switch','oldfileid')
    minChi_switch = 1 ;
    name = 'testInput_minChiON.mat' ;
    save(name,'kmin','kmax','koptmin','koptmax','Switch','minChi_switch','oldfileid')
    e1 = [0,1,0] ;
    e2 = [0,1,0] ;
    e3 = [0,1,0] ;
    e4 = [0,1,0] ;
    A = [ e1; e2; e3 ; e4 ] ;
    name = 'test_initialize_A_testmatrix-n=3-A.mat' ;
    save(name,'A')
    clearvars kmin kmax kopt minChi_switch oldfileid A
%   -----------------------------------------------------------------------    
%   -----------------------------------------------------------------------
%   Trusted results from Markov modeling in multi-precision (mp.Digits=50:  
%   50 digits) with program version PCCAwww_Schur_8 for the count.txt 
%   Butane transition count matrix from Marcus Weber, 
%   Zuse Institute Berlin, 2017.  
%   First optimization was performed with gauss-newton, maxiter=100,
%   xscale=1e-06, xtol=0.0001; second (final) optional optimization was
%   performed with nelder-mead, maxiter=2000, tolfun=1e-08, tolx=1e-08.
    if strcmp(testCase.TestData.precision,'mp')

%           trusted pi vector
        %pi = mp(['[' ...
            %'2.187571738471187133902010209426431107370120314283084e-03 '...    
            %'5.423036997059156482442849543364912557266841765301612e-03 '...   
            %'2.342044854329907260168349964151761987284772748464778e-02 '...   
            %'5.82857795969564539526224197339871353172023700420376e-02 '...   
            %'7.249429129236546999640676258221933776081948593780891e-02 '...    
            %'5.405254454504247218578406828743684371598311397463422e-02 '...   
            %'2.835011453613026984410550759359357413416088406624178e-02 '...   
            %'1.464966954395519941769762620611903248173040191869405e-02 '...   
            %'9.497264358463482494311324131462901182375049908834072e-03 '...   
            %'5.873112283083479772803973321317160093259194654200614e-03 '...   
            %'5.11881979629663654152847775468251262673084951647561e-03 '...   
            %'4.124483560128425064871186093010587059163505318218241e-03 '...   
            %'3.571750407440842955112976496339522747320363072257662e-03 '...   
            %'3.692062289670242853883241202001993248407862279938365e-03 '...   
            %'4.407199335201498052321600670009464455227404169163185e-03 '...   
            %'5.574449257298707584836971904564222834389840901133797e-03 '...   
            %'7.342640276733900326651045301307644057934297857129383e-03 '...   
            %'1.399329698364030539576898074002627013807852032174395e-02 '...   
            %'3.046962550162387511970164593953755562709520736344117e-02 '...   
            %'5.972086356741853855741812980868204329718853313486307e-02 '...   
            %'8.726492361387920248951755364896302567171102784491492e-02 '...   
            %'8.835621020173257643148704284922560801469083260596333e-02 '...   
            %'6.072468748794412023536952401498073278891664703008955e-02 '...   
            %'2.987120408479456813169795751542147880384980005856268e-02 '...   
            %'1.386374881950864541426139670828217022628172043490721e-02 '...   
            %'7.769216656478180839450417795830069882525660099732915e-03 '...   
            %'5.324104870327474208926669153463147380787585372678127e-03 '...   
            %'3.955418928941523649879773437730092147573557883468659e-03 '...   
            %'4.027280031456316269729737741278624605202506059657095e-03 '...   
            %'3.218799581014319343245023602445808184554096362520484e-03 '...   
            %'3.745858750020438513052979905013613620747426334531594e-03 '...   
            %'5.387966522927750746553881183087983429716988885592274e-03 '...   
            %'5.832573831546526856354722524286525993642042410187094e-03 '...   
            %'7.677029504691290542091061441459794707642367788514225e-03 '...   
            %'1.459368537031396585424615476070039944096715840641644e-02 '...   
            %'2.783266461995124875161728265907097088549624852133055e-02 '...   
            %'5.38541628666812454950482564318821030607325668855242e-02 '...   
            %'7.425019666813391557063689249535295617657414516216497e-02 '...   
            %'5.855544688937754390450142827946923895255202586768602e-02 '...   
            %'2.262922654046515828530060428980232174540805949777298e-02 '...   
            %'5.929305078819016054013734113714815933195344238828693e-03 '...   
            %'3.087264670715756079163608287931754032680618247919858e-03 ]']) ;
        %pi = pi' ;
        %testCase.TestData.pi = pi ;
%           trusted pi_weights
        %pi_weights = [ 2.8409817970251405e-01 ; ...
                       %2.8737681725227660e-01 ; ...
                       %4.2852500304520924e-01 ] ;
        %testCase.TestData.pi_weights = pi_weights ;
%           trusted sd vector
        sd = mp(['['...
        '2.119682047692846073089036644503324501324801279808036e-03 '...    
        '5.449182622606609008648702694595810628405739139129141e-03 '...    
        '2.342648602709593560965855121731740238964155376693496e-02 '...    
        '5.822126680997850322451632255161725741138829175623657e-02 '...    
        '7.247912813078038294255861620756886467029945508173788e-02 '...    
        '5.404189371594260860870869369594560815877618357246411e-02 '...
        '2.835574663800429935509673548967654851772234164875261e-02 '...    
        '1.468779683047542868569714542818577213417987301904714e-02 '...    
        '9.46857971304304354346847972804079388091786232065187e-03 '...    
        '5.879118132280157976303554466829975503674448832675084e-03 '...    
        '5.129230615407688846672999050142478628205769134629806e-03 '...    
        '4.149377593360995850622406639004149377593360995850611e-03 '...    
        '3.609458581212818077288406738989151627255911613258018e-03 '...   
        '3.709443583462480627905814127880817877318402239664052e-03 '...   
        '4.409338599210118482227665850122481627755836624506334e-03 '...   
        '5.579163125531170324451332300154976753486976953456981e-03 '...   
        '7.31890216467529870519422086686996950457431385292207e-03  '...  
        '1.396790481427785832125181222816577513372994050892362e-02 '...   
        '3.044543318502224666300054991751237314402839574063884e-02 '...   
        '5.97510373443983402489626556016597510373443983402488e-02 '...   
        '8.726690996350547417887316902464630305454181872719118e-02 '...   
        '8.83067539869019647052942058691196320551917212418139e-02 '...   
        '6.067089936509523571464280357946308053791931210318447e-02 '...   
        '2.982552617107433884917262410638404239364095385692149e-02 '...   
        '1.385792131180322951557266410038494225866120081987701e-02 '...   
        '7.748837674348847672849072639104134379843023546468013e-03 '...   
        '5.379193121031845223216517522371644253361995700644913e-03 '...   
        '3.959406089086637004449332600109983502474628805679159e-03 '...   
        '4.009398590211468279758036294555816627505874118882155e-03 '...   
        '3.239514072889066640003999400089986502024696295555666e-03 '...   
        '3.769434584812278158276258561215817627355896615507664e-03 '...   
        '5.369194620806878968154776783482477628355746638004297e-03 '...   
        '5.809128630705394190871369294605809128630705394190856e-03 '...   
        '7.678848172774083887416887466879968004799280107983785e-03 '...   
        '1.461780732890066490026496025596160575913612958056291e-02 '...   
        '2.78658201269809528570714392841073838924161375793631e-02 '...   
        '5.38919162125681147827825826126081087836824476328552e-02 '...   
        '7.429885517172424136379543068539719042143678448232798e-02 '...   
        '5.8661200819877018447232915062740588911663250512423e-02 '...   
        '2.265660150977353396990451432285157226416037594360853e-02 '...   
        '5.929110633404989251612258161275808628705694145878122e-03 '... 
        '2.989551567264910263460480927860820876868469729540559e-03]']) ;
        sd = sd' ;
        testCase.TestData.sd = sd ;
        name = strcat('count-sd.txt') ;
        mp.write(sd,name)
%           trusted A (nelder-mead, 3 clusters)
        A = [ ...
            2.8405568184913027e-01   2.8755576853720499e-01 ...
            4.2838854961366463e-01 ; ...
           -3.7905560917703573e-01   3.3892409425378806e-01 ...
            4.0131514923247673e-02 ; ...
           -1.9965336957536015e-01  -2.7103148561325063e-01 ...
            4.7068485518861081e-01 ] ;
        testCase.TestData.A = A ;
        name = strcat('count-A.txt') ;
        save(name,'A','-ascii','-double')
%           trusted chi (nelder-mead, 3 clusters)
        chi = [ ... 
            8.5043523873880167e-01   4.3889650578342830e-11 ...
            1.4956476121730866e-01 ; ...
            8.8090007542005921e-01   2.5589115375343415e-02 ...
            9.3510809204597284e-02 ; ...
            9.6495990190978398e-01   1.1558207813298174e-02 ...
            2.3481890276917761e-02 ; ...
            9.7558793310460756e-01   1.5709140129768342e-02 ...
            8.7029267656239927e-03 ; ...
            9.7963327053545801e-01   1.4060862546032220e-02 ...
            6.3058669185096176e-03 ; ...
            9.8568851438933369e-01   1.4311485610666249e-02 ...
            2.8588715185016204e-18 ; ...
            9.6791270535282714e-01   1.7471491809863966e-02 ...
            1.4615802837308787e-02 ; ...
            9.4601281012257943e-01   1.7835098543345616e-02 ...
            3.6152091334074828e-02 ; ...
            9.2448952433822884e-01   1.8479564079043900e-02 ...
            5.7030911582727190e-02 ; ...
            8.8318556579259355e-01   2.4948552633467999e-02 ...
            9.1865881573938366e-02 ; ...
            7.5874311920722926e-01   2.5407453476173496e-02 ...
            2.1584942731659709e-01 ; ...
            6.8245807053420782e-01   5.4751908940179971e-02 ...
            2.6279002052561207e-01 ; ...
            5.3234840569967345e-01   6.3068044907414536e-02 ...
            4.0458354939291197e-01 ; ...
            3.1537515699639257e-01   5.2337908021708332e-02 ...
            6.3228693498189903e-01 ; ...
            2.2250786517110874e-01   7.6139414102714117e-02 ...
            7.0135272072617705e-01 ; ...
            1.2375934853700991e-01   5.4508417089852046e-02 ...
            8.2173223437313792e-01 ; ...
            4.9653667316515840e-02   2.9603583297964756e-02 ...
            9.2074274938551937e-01 ; ...
            2.5736159571560040e-02   1.2561644103436354e-02 ...
            9.6170219632500353e-01 ; ...
            1.1894344993512617e-02   6.3139005313609343e-03 ...
            9.8179175447512634e-01 ; ...
            4.8084968892437262e-03   3.6671524338924661e-03 ...
            9.9152435067686373e-01 ; ...
            3.1536353492936761e-03   7.1221730872170121e-04 ...
            9.9613414734198447e-01 ; ...
            1.7064014946817492e-17   4.7395675902214775e-17 ...
            9.9999999999999978e-01 ; ...
            8.2726996336765120e-03   3.0940472642748782e-04 ...
            9.9141789563989591e-01 ; ...
            1.2722881126783015e-02   5.9164159394738403e-03 ...
            9.8136070293374311e-01 ; ...
            2.2122528189377780e-02   1.3902592390826123e-02 ...
            9.6397487941979598e-01 ; ...
            4.0589834059942542e-02   3.7980881743778290e-02 ...
            9.2142928419627912e-01 ; ...
            8.0588740622014418e-02   1.1344669076732349e-01 ...
            8.0596456861066201e-01 ; ...
            1.0568119504467928e-01   2.0230326690849965e-01 ...
            6.9201553804682103e-01 ; ...
            8.3281552679900978e-02   3.4042640424725945e-01 ...
            5.7629204307283954e-01 ; ...
            5.8911161575871167e-02   4.8967335413071073e-01 ...
            4.5141548429341805e-01 ; ...
            4.3770917595154141e-02   6.6193931626137481e-01 ...
            2.9428976614347097e-01 ; ...
            2.9817138929275813e-02   8.0565415247849947e-01 ...
            1.6452870859222457e-01 ; ...
            1.4136303417923221e-02   8.8989778149855503e-01 ...
            9.5965915083521605e-02 ; ...
            1.3671381908781288e-02   9.4419514302533569e-01 ...
            4.2133475065882837e-02 ; ...
            7.5795321139196827e-03   9.6824758368408359e-01 ...
            2.4172884201996588e-02 ; ...
            1.0405260238826743e-02   9.8297752755905710e-01 ...
            6.6172122021160209e-03 ; ...
            4.5138347223811137e-03   9.8998687175067512e-01 ...
            5.4992935269436351e-03 ; ...
            8.4518297221101869e-03   9.9154817000192241e-01 ...
            2.7596725564344715e-10 ; ...
            6.5843351926753247e-03   9.8840580521018662e-01 ...
            5.0098595971378411e-03 ; ...
            9.2202981926287932e-03   9.7351230297269720e-01 ...
            1.7267398834673802e-02 ; ...
            1.6299801690322059e-02   9.2119142752842509e-01 ...
            6.2508770781252748e-02 ; ...
            4.5801965122805573e-10   9.1035357717143850e-01 ...
            8.9646422370541698e-02] ;
        testCase.TestData.chi = chi ;
        name = strcat('count-chi.txt') ;
        save(name,'chi','-ascii','-double')
%           trusted sd_weights (nelder-mead, 3 clusters)
        sd_weights = [ 2.8405568184913021e-01 ; ... 
                       2.8755576853720499e-01 ; ...
                       4.2838854961366463e-01 ] ;
        testCase.TestData.sd_weights = sd_weights ;
        name = strcat('count-sd_weights.txt') ;
        save(name,'sd_weights','-ascii','-double')
%           trusted Pc (nelder-mead, 3 clusters)
        Pc = [ 9.5104510750059468e-01   1.8354056625839634e-02 ...
               3.0600835873566212e-02 ; ...
               1.8594617368171678e-02   9.5127058926326458e-01 ...
               3.0134793368563506e-02 ; ...
               1.9985610707576054e-02   2.0510702582310622e-02 ...
               9.5950368671011355e-01 ] ;
        testCase.TestData.Pc = Pc ;
        name = strcat('count-Pc.txt') ;
        save(name,'Pc','-ascii','-double')
%           trusted index vector (nelder-mead, 3 clusters)
        idx_vec = [ 1.0000000000000000e+00 ; ...
                    2.0000000000000000e+00 ; ...
                    3.0000000000000000e+00 ; ...
                    4.0000000000000000e+00 ; ...
                    5.0000000000000000e+00 ; ...
                    6.0000000000000000e+00 ; ...
                    7.0000000000000000e+00 ; ...
                    8.0000000000000000e+00 ; ...
                    9.0000000000000000e+00 ; ...
                    1.0000000000000000e+01 ; ...
                    1.1000000000000000e+01 ; ...
                    1.2000000000000000e+01 ; ...
                    1.3000000000000000e+01 ; ...
                    3.0000000000000000e+01 ; ...
                    3.1000000000000000e+01 ; ...
                    3.2000000000000000e+01 ; ...
                    3.3000000000000000e+01 ; ...
                    3.4000000000000000e+01 ; ...
                    3.5000000000000000e+01 ; ...
                    3.6000000000000000e+01 ; ...
                    3.7000000000000000e+01 ; ...
                    3.8000000000000000e+01 ; ...
                    3.9000000000000000e+01 ; ...
                    4.0000000000000000e+01 ; ...
                    4.1000000000000000e+01 ; ...
                    4.2000000000000000e+01 ; ...
                    1.4000000000000000e+01 ; ...
                    1.5000000000000000e+01 ; ...
                    1.6000000000000000e+01 ; ...
                    1.7000000000000000e+01 ; ...
                    1.8000000000000000e+01 ; ...
                    1.9000000000000000e+01 ; ...
                    2.0000000000000000e+01 ; ...
                    2.1000000000000000e+01 ; ...
                    2.2000000000000000e+01 ; ...
                    2.3000000000000000e+01 ; ...
                    2.4000000000000000e+01 ; ...
                    2.5000000000000000e+01 ; ...
                    2.6000000000000000e+01 ; ...
                    2.7000000000000000e+01 ; ...
                    2.8000000000000000e+01 ; ...
                    2.9000000000000000e+01 ] ;
        testCase.TestData.idx_vec = idx_vec ;
%           trusted val_vec (gauss-newton, 2-5 clusters)
        val_vec = [ 2.0000000000000000e+00 ,  4.5068631407622628e-01 ; ...
                    3.0000000000000000e+00 ,  2.0412629223762790e-01 ; ...
                    4.0000000000000000e+00 ,  2.0497574271481929e+00 ; ...
                    5.0000000000000000e+00 ,  2.7047258063177479e+00 ] ;
        testCase.TestData.val_vec = val_vec ;
%           trusted opt_vec (gauss-newton, 2-5 clusters)
        opt_vec = [ 2.0000000000000000e+00 ,  7.7465684296188686e-01 ; ...
                    3.0000000000000000e+00 ,  9.3195790258745737e-01 ; ...
                    4.0000000000000000e+00 ,  4.8756064321295178e-01 ; ...
                    5.0000000000000000e+00 ,  4.5905483873645042e-01 ] ;
        testCase.TestData.opt_vec = opt_vec ;      
%           trusted val_vec (nelder-mead, 2-4 clusters)
        val_vec = [ 2.0000000000000000e+00 ,  4.5068631407622628e-01 ;
                    3.0000000000000000e+00 ,  1.7800054504655316e-01 ;
                    4.0000000000000000e+00 ,  1.1742290998095775e+00 ] ;
        testCase.TestData.val_vec_nm = val_vec ;
%           trusted opt_vec (nelder-mead, 2-4 clusters)
        opt_vec = [ 2.0000000000000000e+00 ,  7.7465684296188686e-01 ;
                    3.0000000000000000e+00 ,  9.4066648498448224e-01 ;
                    4.0000000000000000e+00 ,  7.0644272504760564e-01 ] ;
        testCase.TestData.opt_vec_nm = opt_vec ;
%           trusted A (gauss-newton, 3 clusters)
        A = [ 2.8286204679932703e-01   2.9098466034772902e-01 ...
              4.2615329285294407e-01 ;
             -3.7746276387968841e-01   3.3990446545287428e-01 ...
              3.7558298426814150e-02 ;
             -1.9881438864461831e-01  -2.6648090341850061e-01 ...
              4.6529529206311887e-01 ] ;
        testCase.TestData.A_gn = A ;
%           trusted Pc (gauss-newton, 3 clusters)
        Pc = [ 9.5095884505792783e-01   1.8589883410165424e-02 ...
               3.0451271531903908e-02 ;
               1.8516549628659350e-02   9.5151400755663185e-01 ...
               2.9969442814710048e-02 ;
               1.9914240579993654e-02   2.0739228560594104e-02 ...
               9.5934653085941246e-01 ] ;
        testCase.TestData.Pc_gn = Pc ;
%   -----------------------------------------------------------------------
%   Trusted results from Markov modeling in double precision with program 
%   version PCCAwww_Schur_8 for the count.txt Butane transition matrix 
%   from Marcus Weber, Zuse Institute Berlin, 2017.
%   First optimization was performed with gauss-newton, maxiter=100,
%   xscale=1e-06, xtol=0.0001; second (final) optional optimization was
%   performed with nelder-mead, maxiter=2000, tolfun=1e-08, tolx=1e-08.
    elseif strcmp(testCase.TestData.precision,'double')

%           trusted pi vector
%         pi = [ 2.1875717384711871e-03 ; ...
%                5.4230369970591575e-03 ; ...
%                2.3420448543299083e-02 ; ...
%                5.8285779596956473e-02 ; ...
%                7.2494291292365468e-02 ; ...
%                5.4052544545042477e-02 ; ...
%                2.8350114536130275e-02 ; ...
%                1.4649669543955196e-02 ; ...
%                9.4972643584634829e-03 ; ...
%                5.8731122830834818e-03 ; ...
%                5.1188197962966381e-03 ; ...
%                4.1244835601284245e-03 ; ...
%                3.5717504074408435e-03 ; ...
%                3.6920622896702438e-03 ; ...
%                4.4071993352014999e-03 ; ...
%                5.5744492572987097e-03 ; ...
%                7.3426402767339024e-03 ; ...
%                1.3993296983640316e-02 ; ...
%                3.0469625501623895e-02 ; ...
%                5.9720863567418590e-02 ; ...
%                8.7264923613879344e-02 ; ...
%                8.8356210201732663e-02 ; ...
%                6.0724687487944184e-02 ; ...
%                2.9871204084794590e-02 ; ...
%                1.3863748819508656e-02 ; ...
%                7.7692166564781867e-03 ; ...
%                5.3241048703274755e-03 ; ...
%                3.9554189289415257e-03 ; ...
%                4.0272800314563153e-03 ; ...
%                3.2187995810143183e-03 ; ...
%                3.7458587500204367e-03 ; ...
%                5.3879665229277446e-03 ; ...
%                5.8325738315465198e-03 ; ...
%                7.6770295046912742e-03 ; ...
%                1.4593685370313933e-02 ; ...
%                2.7832664619951177e-02 ; ...
%                5.3854162866681099e-02 ; ...
%                7.4250196668133692e-02 ; ...
%                5.8555446889377381e-02 ; ...
%                2.2629226540465116e-02 ; ...
%                5.9293050788190088e-03 ; ...
%                3.0872646707157544e-03 ] ;
%         testCase.TestData.pi = pi ;
%           trusted pi_weights
%         pi_weights = [ 2.8409817970251633e-01 ; ...
%                        2.8737681725227632e-01 ; ...
%                        4.2852500304520702e-01 ] ;
%         testCase.TestData.pi_weights = pi_weights ;
%           trusted sd vector
        sd = [ 2.1196820476928461e-03 ; ...
               5.4491826226066090e-03 ; ...
               2.3426486027095936e-02 ; ...
               5.8221266809978502e-02 ; ...
               7.2479128130780376e-02 ; ...
               5.4041893715942611e-02 ; ...
               2.8355746638004300e-02 ; ...
               1.4687796830475429e-02 ; ...
               9.4685797130430443e-03 ; ...
               5.8791181322801582e-03 ; ...
               5.1292306154076886e-03 ; ...
               4.1493775933609959e-03 ; ...
               3.6094585812128182e-03 ; ...
               3.7094435834624808e-03 ; ...
               4.4093385992101186e-03 ; ...
               5.5791631255311704e-03 ; ...
               7.3189021646752990e-03 ; ...
               1.3967904814277858e-02 ; ...
               3.0445433185022245e-02 ; ...
               5.9751037344398343e-02 ; ...
               8.7266909963505473e-02 ; ...
               8.8306753986901965e-02 ; ...
               6.0670899365095239e-02 ; ...
               2.9825526171074340e-02 ; ...
               1.3857921311803230e-02 ; ...
               7.7488376743488473e-03 ; ...
               5.3791931210318451e-03 ; ...
               3.9594060890866369e-03 ; ...
               4.0093985902114682e-03 ; ...
               3.2395140728890665e-03 ; ...
               3.7694345848122784e-03 ; ...
               5.3691946208068789e-03 ; ...
               5.8091286307053944e-03 ; ...
               7.6788481727740835e-03 ; ...
               1.4617807328900665e-02 ; ...
               2.7865820126980953e-02 ; ...
               5.3891916212568114e-02 ; ...
               7.4298855171724243e-02 ; ...
               5.8661200819877017e-02 ; ...
               2.2656601509773534e-02 ; ...
               5.9291106334049895e-03 ; ...
               2.9895515672649104e-03 ] ;
        testCase.TestData.sd = sd ;
        name = strcat('count-sd.txt') ;
        save(name,'sd','-ascii','-double')
%           trusted A (nelder-mead, 3 clusters)
        A = [ 2.8405568184913244e-01   2.8755576853720538e-01 ...
              4.2838854961366207e-01 ; ...
              3.7905560917700676e-01  -3.3892409425382630e-01 ...
             -4.0131514923180477e-02 ; ...
             -1.9965336957541363e-01  -2.7103148561320506e-01 ...
              4.7068485518861869e-01 ] ;
        testCase.TestData.A = A ;
        name = strcat('count-A.txt') ;
        save(name,'A','-ascii','-double')
%           trusted chi (nelder-mead, 3 clusters)
        chi = [ 8.5043523873879956e-01   4.3886683668813712e-11 ...
                1.4956476121731360e-01 ; ...
                8.8090007542006055e-01   2.5589115375342489e-02 ...
                9.3510809204596937e-02 ; ...
                9.6495990190978509e-01   1.1558207813296294e-02 ...
                2.3481890276918493e-02 ; ...
                9.7558793310460945e-01   1.5709140129765899e-02 ...
                8.7029267656243241e-03 ; ...
                9.7963327053545934e-01   1.4060862546031219e-02 ...
                6.3058669185093305e-03 ; ...
                9.8568851438933425e-01   1.4311485610665552e-02 ...
                2.8071248236694572e-17 ; ...
                9.6791270535282792e-01   1.7471491809863252e-02 ...
                1.4615802837308636e-02 ; ...
                9.4601281012257998e-01   1.7835098543344384e-02 ...
                3.6152091334075327e-02 ; ...
                9.2448952433822940e-01   1.8479564079042093e-02 ...
                5.7030911582728390e-02 ; ...
                8.8318556579259333e-01   2.4948552633467080e-02 ...
                9.1865881573939379e-02 ; ...
                7.5874311920722926e-01   2.5407453476171762e-02 ...
                2.1584942731659873e-01 ; ...
                6.8245807053420782e-01   5.4751908940178097e-02 ...
                2.6279002052561395e-01 ; ...
                5.3234840569967190e-01   6.3068044907413703e-02 ...
                4.0458354939291430e-01 ; ...
                3.1537515699639085e-01   5.2337908021706729e-02 ...
                6.3228693498190214e-01 ; ...
                2.2250786517110696e-01   7.6139414102712424e-02 ...
                7.0135272072618049e-01 ; ...
                1.2375934853700912e-01   5.4508417089850582e-02 ...
                8.2173223437314014e-01 ; ...
                4.9653667316514778e-02   2.9603583297963372e-02 ...
                9.2074274938552181e-01 ; ...
                2.5736159571559877e-02   1.2561644103436256e-02 ...
                9.6170219632500376e-01 ; ...
                1.1894344993512570e-02   6.3139005313615033e-03 ...
                9.8179175447512568e-01 ; ...
                4.8084968892439985e-03   3.6671524338928122e-03 ...
                9.9152435067686306e-01 ; ...
                3.1536353492936597e-03   7.1221730872186764e-04 ...
                9.9613414734198436e-01 ; ...
               -3.0507380398593992e-18  -6.5549521997357121e-18 ...
                9.9999999999999989e-01 ; ...
                8.2726996336765744e-03   3.0940472642753200e-04 ...
                9.9141789563989580e-01 ; ...
                1.2722881126783353e-02   5.9164159394744483e-03 ...
                9.8136070293374211e-01 ; ...
                2.2122528189377291e-02   1.3902592390825667e-02 ...
                9.6397487941979698e-01 ; ...
                4.0589834059942057e-02   3.7980881743777700e-02 ...
                9.2142928419628023e-01 ; ...
                8.0588740622014210e-02   1.1344669076732210e-01 ...
                8.0596456861066357e-01 ; ...
                1.0568119504467953e-01   2.0230326690849842e-01 ...
                6.9201553804682192e-01 ; ...
                8.3281552679901921e-02   3.4042640424725878e-01 ...
                5.7629204307283921e-01 ; ...
                5.8911161575873422e-02   4.8967335413071028e-01 ...
                4.5141548429341599e-01 ; ...
                4.3770917595157881e-02   6.6193931626137592e-01 ...
                2.9428976614346619e-01 ; ...
                2.9817138929281523e-02   8.0565415247850192e-01 ...
                1.6452870859221638e-01 ; ...
                1.4136303417928717e-02   8.8989778149855780e-01 ...
                9.5965915083513417e-02 ; ...
                1.3671381908787954e-02   9.4419514302533858e-01 ...
                4.2133475065873442e-02 ; ...
                7.5795321139255963e-03   9.6824758368408681e-01 ...
                2.4172884201987460e-02 ; ...
                1.0405260238833199e-02   9.8297752755905998e-01 ...
                6.6172122021064947e-03 ; ...
                4.5138347223880404e-03   9.8998687175067734e-01 ...
                5.4992935269345851e-03 ; ...
                8.4518297221162914e-03   9.9154817000192508e-01 ...
                2.7595845559331095e-10 ; ...
                6.5843351926822497e-03   9.8840580521018928e-01 ...
                5.0098595971283938e-03 ; ...
                9.2202981926351631e-03   9.7351230297269997e-01 ...
                1.7267398834664754e-02 ; ...
                1.6299801690329144e-02   9.2119142752842731e-01 ...
                6.2508770781243561e-02 ; ...
                4.5803069549102952e-10   9.1035357717144483e-01 ...
                8.9646422370524254e-02 ] ;
        testCase.TestData.chi = chi ;
        name = strcat('count-chi.txt') ;
        save(name,'chi','-ascii','-double')
%           trusted sd_weights (nelder-mead, 3 clusters)
        sd_weights = [ 2.8405568184913244e-01 ; ...
                       2.8755576853720533e-01 ; ...
                       4.2838854961366213e-01 ] ;
        testCase.TestData.sd_weights = sd_weights ;
        name = strcat('count-sd_weights.txt') ;
        save(name,'sd_weights','-ascii','-double')
%           trusted Pc (nelder-mead, 3 clusters)
        Pc = [ 9.5104510750059479e-01   1.8354056625839690e-02 ...
               3.0600835873566056e-02 ; ...
               1.8594617368171685e-02   9.5127058926326480e-01 ...
               3.0134793368563478e-02 ; ...
               1.9985610707576221e-02   2.0510702582310643e-02 ...
               9.5950368671011388e-01 ] ;
        testCase.TestData.Pc = Pc ;
        name = strcat('count-Pc.txt') ;
        save(name,'Pc','-ascii','-double')
%           trusted index vector (nelder-mead, 3 clusters)
        idx_vec = [ 1.0000000000000000e+00 ; ...
                    2.0000000000000000e+00 ; ...
                    3.0000000000000000e+00 ; ...
                    4.0000000000000000e+00 ; ...
                    5.0000000000000000e+00 ; ...
                    6.0000000000000000e+00 ; ...
                    7.0000000000000000e+00 ; ...
                    8.0000000000000000e+00 ; ...
                    9.0000000000000000e+00 ; ...
                    1.0000000000000000e+01 ; ...
                    1.1000000000000000e+01 ; ...
                    1.2000000000000000e+01 ; ...
                    1.3000000000000000e+01 ; ...
                    3.0000000000000000e+01 ; ...
                    3.1000000000000000e+01 ; ...
                    3.2000000000000000e+01 ; ...
                    3.3000000000000000e+01 ; ...
                    3.4000000000000000e+01 ; ...
                    3.5000000000000000e+01 ; ...
                    3.6000000000000000e+01 ; ...
                    3.7000000000000000e+01 ; ...
                    3.8000000000000000e+01 ; ...
                    3.9000000000000000e+01 ; ...
                    4.0000000000000000e+01 ; ...
                    4.1000000000000000e+01 ; ...
                    4.2000000000000000e+01 ; ...
                    1.4000000000000000e+01 ; ...
                    1.5000000000000000e+01 ; ...
                    1.6000000000000000e+01 ; ...
                    1.7000000000000000e+01 ; ...
                    1.8000000000000000e+01 ; ...
                    1.9000000000000000e+01 ; ...
                    2.0000000000000000e+01 ; ...
                    2.1000000000000000e+01 ; ...
                    2.2000000000000000e+01 ; ...
                    2.3000000000000000e+01 ; ...
                    2.4000000000000000e+01 ; ...
                    2.5000000000000000e+01 ; ...
                    2.6000000000000000e+01 ; ...
                    2.7000000000000000e+01 ; ...
                    2.8000000000000000e+01 ; ...
                    2.9000000000000000e+01 ] ;
        testCase.TestData.idx_vec = idx_vec ;
%           trusted val_vec (gauss-newton, 2-5 clusters)
        val_vec = [ 2.0000000000000000e+00 ,  4.5068631407624360e-01 ; ...
                    3.0000000000000000e+00 ,  2.0412629225068146e-01 ; ...
                    4.0000000000000000e+00 ,  2.0497563172193218e+00 ; ...
                    5.0000000000000000e+00 ,  2.7047259529717986e+00 ] ;
        testCase.TestData.val_vec = val_vec ;
%           trusted opt_vec (gauss-newton, 2-5 clusters)
        opt_vec = [ 2.0000000000000000e+00 ,  7.7465684296187820e-01 ; ...
                    3.0000000000000000e+00 ,  9.3195790258310618e-01 ; ...
                    4.0000000000000000e+00 ,  4.8756092069516954e-01 ; ...
                    5.0000000000000000e+00 ,  4.5905480940564025e-01 ] ;
        testCase.TestData.opt_vec = opt_vec ; 
 %           trusted val_vec (nelder-mead, 2-4 clusters)
        val_vec = [ 2.0000000000000000e+00 ,  4.5068631407624360e-01 ;
                    3.0000000000000000e+00 ,  1.7800054504654916e-01 ;
                    4.0000000000000000e+00 ,  1.1742290998137550e+00 ] ;
        testCase.TestData.val_vec_nm = val_vec ;
%           trusted opt_vec (nelder-mead, 2-4 clusters)
        opt_vec = [ 2.0000000000000000e+00 ,  7.7465684296187820e-01 ;
                    3.0000000000000000e+00 ,  9.4066648498448358e-01 ;
                    4.0000000000000000e+00 ,  7.0644272504656125e-01 ] ;
        testCase.TestData.opt_vec_nm = opt_vec ;  
%           trusted A (gauss-newton, 3 clusters)
        A = [ 2.8286204679941895e-01   2.9098466035046661e-01 ...
              4.2615329285011433e-01 ;
              3.7746276388008371e-01  -3.3990446545367009e-01 ...
             -3.7558298426413637e-02 ;
             -1.9881438864515949e-01  -2.6648090341476616e-01 ...
              4.6529529205992559e-01 ] ;
        testCase.TestData.A_gn = A ;
%           trusted Pc (gauss-newton, 3 clusters)
        Pc = [ 9.5095884505793293e-01   1.8589883410363803e-02 ...
               3.0451271531701438e-02 ;
               1.8516549628664006e-02   9.5151400755682702e-01 ...
               2.9969442814510940e-02 ;
               1.9914240580011240e-02   2.0739228560775341e-02 ...
               9.5934653085921318e-01 ] ;
        testCase.TestData.Pc_gn = Pc ;
    end
%   -----------------------------------------------------------------------
%   example transition matrix from Susanna Roeblitz,
%   Zuse Institute Berlin, 2017.
    for i = [ 0, 10, 50, 100, 200, 500, 1000 ]
        mu = i ;
        W =[1000 100 100 10 0 0 0 0 0;...
            100 1000 100 0 0 0 0 0 0;...
            100 100 1000 0 mu 0 0 0 0;...
            10 0 0 1000 100 100 10 0 0;...
            0 0 mu 100 1000 100 0 0 0;...
            0 0 0 100 100 1000 0 mu 0;...
            0 0 0 10 0 0 1000 100 100;...
            0 0 0 0 0 mu 100 1000 100;...
            0 0 0 0 0 0 100 100 1000];
        name = strcat('example_matrix_mu',num2str(i),'.txt') ;
        save(name,'-ascii','W','-double')
    end

    if strcmp(testCase.TestData.precision,'mp')
%       Trusted results from testing with test_fillA(testCase) in   
%       multi-precision (mp.Digits=50: 50 digits) with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=0 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%           trusted initial A (3,3)-matrix
        testCase.TestData.A_mu0_init = [ 3.3390904470272770e-01 ...
            3.3304547764863618e-01   3.3304547764863618e-01 ; ...
           -7.6900386296779886e-49   4.0366140439774489e-01 ...
           -4.0366140439774489e-01 ; ...
            4.5730890701570093e-01  -2.2865445350785046e-01 ...
           -2.2865445350785046e-01 ] ;
%           trusted A (3,3)-matrix after using fillA() on A_init
        testCase.TestData.A_mu0 = [ 3.3390904470272770e-01 ...
            3.3304547764863618e-01   3.3304547764863618e-01 ; ...
           -0.0000000000000000e+00   4.0366140439774489e-01 ...
           -4.0366140439774489e-01 ; ...
            4.5730890701570093e-01  -2.2865445350785046e-01 ...
           -2.2865445350785046e-01 ] ;
%           trusted (cut to double precision by double()) schur vector (9,3)-matrix
        v1 = [ 1.0000000000000000e+00   1.1983873824613118e+00 ...
              -6.6212084649124303e-01 ] ;
        v2 = [ 1.0000000000000000e+00   1.2386618947283068e+00 ...
              -7.3016081598267124e-01 ] ;
        v3 = [ 1.0000000000000000e+00   1.2386618947283068e+00 ...
              -7.3016081598267124e-01 ] ;
        v4 = [ 1.0000000000000000e+00  -2.2437686940366006e-48 ...
               1.3208168935772102e+00 ] ;
        v5 = [ 1.0000000000000000e+00  -2.4972075324166187e-48 ...
               1.4565448979422639e+00 ] ;
        v6 = [ 1.0000000000000000e+00  -2.4731316796753597e-48 ...
               1.4565448979422639e+00 ] ;
        v7 = [ 1.0000000000000000e+00  -1.1983873824613118e+00 ...
              -6.6212084649124303e-01 ] ;
        v8 = [ 1.0000000000000000e+00  -1.2386618947283068e+00 ...
              -7.3016081598267124e-01 ] ;
        v9 = [ 1.0000000000000000e+00  -1.2386618947283068e+00 ...
              -7.3016081598267124e-01 ] ;
        testCase.TestData.svecs_mu0 = [ v1; v2; v3; v4; v5; v6; v7; v8; v9 ] ; 
%       Trusted results from testing with test_fillA(testCase) in   
%       multi-precision (mp.Digits=50: 50 digits) with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=1000 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%           trusted initial A (5,5)-matrix
        testCase.TestData.A_mu1000_init = [ ...
            9.1894070188323768e-02   1.7530349447685256e-01 ...
            1.7530349447685256e-01   2.7874947042898557e-01 ...
            2.7874947042898557e-01 ; ...
            9.5617462375120710e-51   2.3004217991775647e-01 ...
           -2.3004217991775647e-01   2.4680580481257491e-01 ...
           -2.4680580481257491e-01 ; ...
           -8.2063944467677397e-02   2.1685904674531994e-01 ...
            2.1685904674531994e-01  -1.7582707451148125e-01 ...
           -1.7582707451148125e-01 ; ...
            1.2649128975068150e-51  -1.5813354667522092e-01 ...
            1.5813354667522092e-01   3.0646480920422281e-01 ...
           -3.0646480920422281e-01 ; ...
           -2.6659946136620422e-01  -4.6750705493354562e-02 ...
           -4.6750705493354562e-02   1.8005043617645666e-01 ...
            1.8005043617645666e-01 ] ;
%           trusted A (5,5)-matrix after using fillA() on A_init
        testCase.TestData.A_mu1000 = [ ...
            1.1624352928638547e-01   1.4701250216352996e-01 ...
            1.4701250216352996e-01   2.9486573319327730e-01 ...
            2.9486573319327730e-01 ; ...
           -0.0000000000000000e+00   1.9291729793400017e-01 ...
           -1.9291729793400017e-01   2.0697555985554739e-01 ...
           -2.0697555985554739e-01 ; ...
           -6.8820224317862874e-02   1.8186169747481573e-01 ...
            1.8186169747481573e-01  -1.4745158531588429e-01 ...
           -1.4745158531588429e-01 ; ...
           -0.0000000000000000e+00  -1.3261349091810171e-01 ...
            1.3261349091810171e-01   2.5700661906731492e-01 ...
           -2.5700661906731492e-01 ; ...
           -2.2357485803610264e-01  -3.9205939465147714e-02 ...
           -3.9205939465147714e-02   1.5099336848319903e-01 ...
            1.5099336848319903e-01 ] ;
%           trusted (cut to double precision by double()) schur vector (9,5)-matrix
        testCase.TestData.svecs_mu1000 = [ 1.0000000000000000e+00 ...
           -1.3488263927262389e+00   1.3104358647426715e+00 ...
            1.0433110064103206e+00  -6.0005035322490063e-01 ; ...
            1.0000000000000000e+00  -1.3990235434236398e+00 ...
            1.4737789968291537e+00   1.1266779128181024e+00 ...
           -1.0896513975018773e-01 ; ...
            1.0000000000000000e+00  -8.8508490246081217e-01 ...
           -4.3510374799711482e-01  -9.8936037088281537e-01 ...
            6.5386362361676453e-01 ; ...
            1.0000000000000000e+00   2.4862023997663969e-49 ...
           -1.4466972854287909e+00  -2.2942203823682567e-49 ...
           -2.9609371304557914e+00 ; ...
            1.0000000000000000e+00  -7.2188567221853384e-01 ...
           -6.8838572855836444e-01  -1.0501513257627346e+00 ...
            5.5658671487432088e-01 ; ...
            1.0000000000000000e+00   7.2188567221853384e-01 ...
           -6.8838572855836444e-01   1.0501513257627346e+00 ...
            5.5658671487432088e-01 ; ...
            1.0000000000000000e+00   1.3488263927262389e+00 ...
            1.3104358647426715e+00  -1.0433110064103206e+00 ...
           -6.0005035322490063e-01 ; ...
            1.0000000000000000e+00   8.8508490246081217e-01 ...
           -4.3510374799711482e-01   9.8936037088281537e-01 ...
            6.5386362361676453e-01 ; ...
            1.0000000000000000e+00   1.3990235434236398e+00 ...
            1.4737789968291537e+00  -1.1266779128181024e+00 ...
           -1.0896513975018773e-01 ] ; 
%       Trusted results from testing with test_use_minChi(testCase) in   
%       multi-precision (mp.Digits=50: 50 digits) with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=0 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%       The results are trusted since they are (visually) similar to the
%       results (for the same test matrix but obtained by - more precise -
%       Nelder-Mead (instead of unoptimized calculation by the Inner 
%       Simplex Algorithm as in use_minChi() and cluster_by_isa()) 
%       optimization) published in Fig. 4 a) from
%       (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+:
%       Application to Markov State Models and Data Classification. Adv.
%       Data Anal. Classif. 2013, 7 (2), 147?179.
%           trusted chi by ISA for n_cluster=3
        testCase.TestData.chi_isa_mu0_n3 = [ ...
            3.1115284081506680e-02   9.6818509177611722e-01   6.9962414237611918e-04 ; ...
            1.8709352970645370e-50   1.0000000000000000e+00  -9.6219529563319043e-50 ; ...
            3.2073176521106348e-50   1.0000000000000000e+00  -9.7555911918365141e-50 ; ...
            9.3793037467239504e-01   3.1034812663802466e-02   3.1034812663802466e-02 ; ...
            1.0000000000000000e+00  -1.8976629441654589e-49   2.0313011796700687e-49 ; ...
            1.0000000000000000e+00  -1.8442076499636150e-49   1.9243905912663809e-49 ; ...
            3.1115284081506680e-02   6.9962414237611918e-04   9.6818509177611722e-01 ; ...
            1.3363823550460978e-50   8.9537617788088554e-50   1.0000000000000000e+00 ; ...
            1.3363823550460978e-50   8.9537617788088554e-50   1.0000000000000000e+00  ] ;
%       Trusted results from testing with test_use_minChi(testCase) in   
%       multi-precision (mp.Digits=50: 50 digits) with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=100 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%       The results are trusted since they are (visually) similar to the
%       results (for the same test matrix but obtained by - more precise -
%       Nelder-Mead (instead of unoptimized calculation by the Inner 
%       Simplex Algorithm as in use_minChi() and cluster_by_isa()) 
%       optimization) published in Fig. 4 b) from
%       (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+:
%       Application to Markov State Models and Data Classification. Adv.
%       Data Anal. Classif. 2013, 7 (2), 147?179.
%           trusted chi by ISA for n_cluster=3
        testCase.TestData.chi_isa_mu100_n3 = [ ...
           -4.2580228676557226e-03   4.2801027432807288e-02   9.6145699543484842e-01 ; ...
           -9.3546764853226848e-51  -1.0691058840368783e-50   1.0000000000000000e+00 ; ...
           -2.7309451028470270e-02   2.5981853655796699e-01   7.6749091447050322e-01 ; ...
           -5.3455294201843913e-51   1.0000000000000000e+00  -2.4054882390829761e-50 ; ...
           -2.3685358779865925e-02   8.3488897155844621e-01   1.8879638722141970e-01 ; ...
            1.8879638722141970e-01   8.3488897155844621e-01  -2.3685358779865925e-02 ; ...
            9.6145699543484842e-01   4.2801027432807288e-02  -4.2580228676557226e-03 ; ...
            7.6749091447050322e-01   2.5981853655796699e-01  -2.7309451028470270e-02 ; ...
            1.0000000000000000e+00  -8.0182941302765869e-51  -4.0091470651382935e-51 ] ;
    elseif strcmp(testCase.TestData.precision,'double')
%       Trusted results from testing with test_fillA(testCase) in   
%       double precision with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=0 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%           trusted initial A (3,3)-matrix
        testCase.TestData.A_mu0_init = [ ...
            3.3390904470272748e-01   3.3304547764863601e-01 ...
            3.3304547764863596e-01 ; ...
            1.1824450993798861e-14  -4.0366140439775017e-01 ...
            4.0366140439773823e-01 ; ...
           -4.5730890701570082e-01   2.2865445350784039e-01 ...
            2.2865445350786040e-01 ] ;
%           trusted A (3,3)-matrix after using fillA() on A_init
        testCase.TestData.A_mu0 = [ ...
            3.3390904470272720e-01   3.3304547764863651e-01 ...
            3.3304547764863612e-01 ; ...
            1.1934897514720417e-14  -4.0366140439774961e-01 ...
            4.0366140439773768e-01 ; ...
           -4.5730890701570021e-01   2.2865445350784008e-01 ...
            2.2865445350786009e-01 ] ;
%           trusted (cut to double precision by double()) schur vector (9,3)-matrix
        testCase.TestData.svecs_mu0 = [ ...
            1.0000000000000002e+00   1.1983873824612949e+00 ...
            6.6212084649127290e-01 ; ...
            1.0000000000000002e+00   1.2386618947282897e+00 ...
            7.3016081598270310e-01 ; ...
            1.0000000000000002e+00   1.2386618947282877e+00 ...
            7.3016081598270255e-01 ; ...
            1.0000000000000002e+00   2.9934085326177830e-14 ...
           -1.3208168935772115e+00 ; ...
            1.0000000000000002e+00   3.8344472783066848e-14 ...
           -1.4565448979422635e+00 ; ...
            1.0000000000000002e+00   3.8841087123289861e-14 ...
           -1.4565448979422631e+00 ; ...
            1.0000000000000002e+00  -1.1983873824613247e+00 ...
            6.6212084649121505e-01 ; ...
            1.0000000000000002e+00  -1.2386618947283270e+00 ...
            7.3016081598263805e-01 ; ...
            1.0000000000000002e+00  -1.2386618947283281e+00 ...
            7.3016081598263893e-01 ] ; 
%       Trusted results from testing with test_fillA(testCase) in   
%       double precision with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=1000 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%           trusted initial A (5,5)-matrix
        testCase.TestData.A_mu1000_init = [ ...
            9.1894070188323573e-02   1.7530349447685195e-01 ...
            1.7530349447685301e-01   2.7874947042898623e-01 ...
            2.7874947042898501e-01 ; ...
           -1.9968968873256783e-16   2.3004217991775602e-01 ...
           -2.3004217991775622e-01   2.4680580481257552e-01 ...
           -2.4680580481257508e-01 ; ...
            8.2063944467677424e-02  -2.1685904674531978e-01 ...
           -2.1685904674531978e-01   1.7582707451147850e-01 ...
            1.7582707451148361e-01 ; ...
           -2.6690050013010803e-16  -1.5813354667522189e-01 ...
            1.5813354667521951e-01   3.0646480920422386e-01 ...
           -3.0646480920422126e-01 ; ...
           -2.6659946136620399e-01  -4.6750705493354333e-02 ...
           -4.6750705493354631e-02   1.8005043617645519e-01 ...
            1.8005043617645786e-01 ] ;
%           trusted A (5,5)-matrix after using fillA() on A_init
        testCase.TestData.A_mu1000 = [ ...
            1.1624352928638476e-01   1.4701250216352987e-01 ...
            1.4701250216353140e-01   2.9486573319327480e-01 ...
            2.9486573319327908e-01 ; ...
           -2.0948669444255076e-16   1.9291729793399881e-01 ...
           -1.9291729793399898e-01   2.0697555985554686e-01 ...
           -2.0697555985554650e-01 ; ...
            6.8820224317862555e-02  -1.8186169747481468e-01 ...
           -1.8186169747481468e-01   1.4745158531588126e-01 ...
            1.4745158531588554e-01 ; ...
           -1.8621039506004511e-16  -1.3261349091810185e-01 ...
            1.3261349091809985e-01   2.5700661906731453e-01 ...
           -2.5700661906731237e-01 ; ...
           -2.2357485803610139e-01  -3.9205939465147326e-02 ...
           -3.9205939465147575e-02   1.5099336848319703e-01 ...
            1.5099336848319928e-01 ] ;
%           trusted (cut to double precision by double()) schur vector (9,5)-matrix
        testCase.TestData.svecs_mu1000 = [ ...
            1.0000000000000002e+00   1.3488263927262381e+00 ...
           -1.3104358647426639e+00  -1.0433110064103177e+00 ...
           -6.0005035322489464e-01 ; ...
            1.0000000000000002e+00   1.3990235434236424e+00 ...
           -1.4737789968291495e+00  -1.1266779128181197e+00 ...
           -1.0896513975018769e-01 ; ...
            1.0000000000000002e+00   8.8508490246081195e-01 ...
            4.3510374799711093e-01   9.8936037088282069e-01 ...
            6.5386362361676154e-01 ; ...
            1.0000000000000002e+00  -1.5699540956179371e-15 ...
            1.4466972854287901e+00  -2.8905778341536754e-15 ...
           -2.9609371304557919e+00 ; ...
            1.0000000000000002e+00   7.2188567221853306e-01 ...
            6.8838572855835944e-01   1.0501513257627397e+00 ...
            5.5658671487431721e-01 ; ...
            1.0000000000000002e+00  -7.2188567221853406e-01 ...
            6.8838572855836977e-01  -1.0501513257627297e+00 ...
            5.5658671487432432e-01 ; ...
            1.0000000000000002e+00  -1.3488263927262389e+00 ...
           -1.3104358647426779e+00   1.0433110064103133e+00 ...
           -6.0005035322490197e-01 ; ...
            1.0000000000000002e+00  -8.8508490246081095e-01 ...
            4.3510374799711915e-01  -9.8936037088281159e-01 ...
            6.5386362361676742e-01 ; ...
            1.0000000000000002e+00  -1.3990235434236402e+00 ...
           -1.4737789968291592e+00   1.1266779128180928e+00 ...
           -1.0896513975019061e-01 ] ; 
%       Trusted results from testing with test_use_minChi(testCase) in   
%       double precision with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=0 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017. 
%       The results are trusted since they are (visually) similar to the
%       results (for the same test matrix but obtained by - more precise -
%       Nelder-Mead (instead of unoptimized calculation by the Inner 
%       Simplex Algorithm as in use_minChi() and cluster_by_isa()) 
%       optimization) published in Fig. 4 a) from
%       (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+:
%       Application to Markov State Models and Data Classification. Adv.
%       Data Anal. Classif. 2013, 7 (2), 147?179.
%           trusted chi by ISA for n_cluster=3
        testCase.TestData.chi_isa_mu0_n3 = [ ...
            3.1115284081507055e-02   6.9962414237666952e-04   9.6818509177611567e-01 ; ...
           -2.2306101257299255e-17   2.7238701370796497e-16   9.9999999999999922e-01 ; ...
            2.3155133771019864e-16   9.5036998707745997e-16   9.9999999999999833e-01 ; ...
            9.3793037467239537e-01   3.1034812663803282e-02   3.1034812663800999e-02 ; ...
            9.9999999999999956e-01  -8.3121127279257571e-16   9.4045844324124387e-16 ; ...
            9.9999999999999933e-01  -8.9620175089935441e-16   1.2085348725220209e-15 ; ...
            3.1115284081505195e-02   9.6818509177611567e-01   6.9962414237894190e-04 ; ...
            4.1989789562932371e-16   9.9999999999999956e-01  -3.2491549094043100e-16 ; ...
            1.3725993281327047e-17   1.0000000000000002e+00  -5.6591874961648636e-16 ] ;
%       Trusted results from testing with test_use_minChi(testCase) in   
%       double precision with program version   
%       PCCAwww_Schur_9 for the example transition matrix with mu=100 from 
%       Susanna Roeblitz, Zuse Institute Berlin, 2017.
%       The results are trusted since they are (visually) similar to the
%       results (for the same test matrix but obtained by - more precise -
%       Nelder-Mead (instead of unoptimized calculation by the Inner 
%       Simplex Algorithm as in use_minChi() and cluster_by_isa()) 
%       optimization) published in Fig. 4 b) from
%       (1) Roeblitz, S.; Weber, M. Fuzzy Spectral Clustering by PCCA+:
%       Application to Markov State Models and Data Classification. Adv.
%       Data Anal. Classif. 2013, 7 (2), 147?179.
%           trusted chi by ISA for n_cluster=3
        testCase.TestData.chi_isa_mu100_n3 = [ ... 
           -4.2580228676599614e-03   4.2801027432808406e-02   9.6145699543485141e-01 ; ...
           -1.4772846838530374e-15   5.0378037267330542e-17   1.0000000000000013e+00 ; ...
           -2.7309451028471876e-02   2.5981853655796916e-01   7.6749091447050255e-01 ; ...
            2.5476043241765107e-15   1.0000000000000000e+00  -2.7488885789634078e-15 ; ...
           -2.3685358779863556e-02   8.3488897155844455e-01   1.8879638722141870e-01 ; ...
            1.8879638722141809e-01   8.3488897155844921e-01  -2.3685358779867570e-02 ; ...
            9.6145699543484642e-01   4.2801027432807524e-02  -4.2580228676540859e-03 ; ...
            7.6749091447050422e-01   2.5981853655796533e-01  -2.7309451028469895e-02 ; ...
            9.9999999999999867e-01   1.8729778932510503e-17   1.1929305952124212e-15 ] ;
    end
%   -----------------------------------------------------------------------
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function teardownOnce(testCase)
% test environment tear down
%   Written by Bernhard Reuter, Theoretical Physics II,
%   University of Kassel, 2017 
    cd(testCase.TestData.origPath)
%       clear all variables, globals and functions from workspace
    clear all
    %rmdir(testCase.TestData.tmpFolder)
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [ P, sd ] = get_knownInput( matrixfile )
% Load the count transition matrix for Butane from Marcus Weber and
% calculate the stochastic transition matrix and the "initial distribution"
% from it. 
%   The loading and calculations are performed in double or
%   multiprecision, depending on the user input from keyboard (see
%   setupOnce() function).
%   Written by Bernhard Reuter, Theoretical Physics II,
%   University of Kassel, 2017
%   ------------------------------------------------------------
%   read the count matrix from file, calculate the stochastic matrix,
    Tc = load_t(matrixfile,'-ascii',numeric_t) ;
    dummy = (mod(Tc,1) ~= 0) ;
    assert(~any(dummy(:)), ...
        'main:Tc_DataError','Tc doesnt seem to be a count matrix')
    clearvars dummy
    assert(isa(Tc,numeric_t),'main:Tc_DataTypeError', ...
        'Variable is type %s not %s',class(Tc),numeric_t)
    assert(size(Tc,1)==size(Tc,2),'main:Tc_MatrixShapeError', ...
        'Matrix is not quadratic but %d x %d',size(Tc,1),size(Tc,2))
%       find rows of the count matrix Tc with nonzero rowsum
    [row, ~] = find(sum(Tc,2) > numeric_t('0.0001')) ;
%       calculate stochastic matrix P from the count matrix Tc
%       while accounting for rows in Tc with zero rowsum
    P = Tc ;
    Ts = Tc(row, :) ;
    P(row, :) = diag(numeric_t('1.0')./sum(Ts,2))*Ts ;
    assert(isa(P,numeric_t),'main:P_DataTypeError', ...
        'Variable is type %s not %s',class(P),numeric_t)
    %   calculate "initial distribution"
    sd = sum(Tc,2) ; sd=sd/sum(sd) ;
    assert(isa(sd,numeric_t),'main:sd_DataTypeError', ...
        'Variable is type %s not %s',class(sd),numeric_t)
end
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
function [ actVal, expVal, c ] = verifyAlmostEqual(actVal, expVal, ...
    abstol, reltol, varargin)
% Checks if an actual value actVal and an expected value expVal are almost 
%   equal with respect to a given absolute or relative tolerance or both. 
%   The precision can bei either mp or double.
%   Written by Bernhard Reuter, Theoretical Physics II, University of
%   Kassel, 2017
%
%   Input:
%       actVal      the actual value
%       expVal      the expected value
%       abstol      the absolute tolerance, will be used only if abstol>0
%       reltol      the relative tolerance, will be used only if reltol>0
%       varargin:
%       combi       logical; if true, both abstol and reltol have to be
%                   fullfilled; if false, either abstol or reltol have to
%                   be fullfilled
%   -----------------------------------------------------------------------   
%       catch inappropriate use
    assert( (isa(actVal,'mp') || isa(actVal,'double')), ...
        'verifyAlmostEqual:DataTypeError1', ...
        'actVal is neither mp nor double')
    assert( (isa(expVal,'mp') || isa(expVal,'double')), ...
        'verifyAlmostEqual:DataTypeError2', ...
        'expVal is neither mp nor double')
    assert(all(size(actVal)==size(expVal)), ...
        'verifyAlmostEqual:MatchError', ...
        'expVal and actVal do not have the same shape')
%   -----------------------------------------------------------------------
    class_t1 = class(actVal);
    class_t2 = class(expVal);
    if (strcmp(class_t1,class_t2)) % l2
        class_t = class_t1;
    else % l1
        disp(['class(actVal)=' class_t1 ' is not equal class(expVal)=' class_t2 '!'])
        disp('will convert to datatype double before comparison ')
        actVal = double(actVal) ;
        expVal = double(expVal) ;
        class_t = 'double' ;
    end
%   -----------------------------------------------------------------------
    function r = num_t(expression)
        if (nargin > 0)
            if(strcmp(class_t,'mp')), r = mp(expression);
            else
                if isnumeric(expression)
                    r = expression;
                else
                    r = eval(expression);
                end;
            end;
        else
            r = class_t;
        end;
    end  % num_t
%   -----------------------------------------------------------------------
%   check for correct use of varargin
    if nargin == 5 % l3
        assert(isa(varargin{1}, 'logical'), 'verifyAlmostEqual:DataTypeError3', ...
        'combi is not of class logical')
        combi = varargin{1} ;
    elseif nargin > 5 % l4
        error('verifyAlmostEqual:InputError','Too many input arguments!')
    end
%   -----------------------------------------------------------------------
    absA = abs(actVal) ;
    absB = abs(expVal) ;
	diff = abs(actVal - expVal) ;
    idx = (actVal == expVal) ; % to handle elementwise infinities/equalities
    diff(idx) = num_t('0') ; % to handle elementwise infinities/equalities
    

%       shortcut, handles infinities/equalities, if ALL elements of actVal
%       and expVal are infinite or equal
    if (actVal == expVal) % l5
        c = true ;
        return
    elseif nargin == 5 % l6
%           use absolute and relative error, almost equal if both are
%           fullfilled
        if ( (combi == true) && (abstol > num_t('0'))  && (reltol > num_t('0')) ) % l7+l8+l9
            c1 = ( diff < abstol ) ;
            c1 = all(c1(:)) ;
            relerr = diff ./ min(absA + absB, realmax(num_t)) ;
            idx = (actVal == expVal) ; % to handle elementwise infinities/equalities
            relerr(idx) = num_t('0') ; % to handle elementwise infinities/equalities
            c2 = ( relerr < reltol ) ;
            c2 = all(c2(:)) ;
            c = ( c1 == true && c2 == true ) ;
            if c == false % l10
                disp('Expected value:')
                disp(expVal)
                disp('Actual value:')
                disp(actVal)
                disp('abstol:')
                disp(abstol)
                disp('reltol:')
                disp(reltol)
                disp('absolute error:')
                disp(diff)
                disp('relative error:')
                disp(relerr)
            end 
            return
        % use absolute and relative error, almost equal if at least one is
        % fullfilled
        elseif ( (combi == false) && (abstol > num_t('0'))  && (reltol > num_t('0')) ) % l11+l12+l13
            c1 = ( diff < abstol ) ;
            c1 = all(c1(:)) ;
            relerr = diff ./ min(absA + absB, realmax(num_t)) ;
            idx = (actVal == expVal) ; % to handle elementwise infinities/equalities
            relerr(idx) = num_t('0') ; % to handle elementwise infinities/equalities
            c2 = ( relerr < reltol ) ;
            c2 = all(c2(:)) ;
            c = ( c1 == true || c2 == true ) ;
            if c == false % l14
                disp('Expected value:')
                disp(expVal)
                disp('Actual value:')
                disp(actVal)
                disp('abstol:')
                disp(abstol)
                disp('reltol:')
                disp(reltol)
                disp('absolute error:')
                disp(diff)
                disp('relative error:')
                disp(relerr)
            end
            return
        else % l6
            error('verifyAlmostEqual:FatalError1','Nothing applied (1)!')
        end
%       use absolute error 
    elseif abstol > num_t('0') % l15
        c = ( diff < abstol ) ;
        c = all(c(:)) ;
        if c == false % l16
            disp('Expected value:')
            disp(expVal)
            disp('Actual value:')
            disp(actVal)
            disp('abstol:')
            disp(abstol)
            disp('absolute error:')
            disp(diff)
        end
        return
%       use relative error 
    elseif reltol > num_t('0') % l17
        relerr = diff ./ min(absA + absB, realmax(num_t)) ;
        idx = (actVal == expVal) ; % to handle elementwise infinities/equalities
        relerr(idx) = num_t('0') ; % to handle elementwise infinities/equalities
        c = ( relerr < reltol ) ;
        c = all(c(:)) ;
        if c == false % l18
            disp('Expected value:')
            disp(expVal)
            disp('Actual value:')
            disp(actVal)
            disp('reltol:')
            disp(reltol)
            disp('relative error:')
            disp(relerr)
        end
        return
    else % l2
        error('verifyAlmostEqual:FatalError2','Nothing applied (2)!')
    end
end
% end of function verifyAlmostEqual()
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
