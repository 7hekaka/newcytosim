function [sys, rhs, con] = cytosim_dump(path)

% This is used to explore Cytosim's linear system in matlab
% - load the matrices and vector from Cytosim's dump
% - plot convergence pattern of BICGstab, with and without preconditionning
%
% F. Nedelec, 16.10.2014, 03.2018, 06.2018, 26.01.2019, 30.06.2019,
% 11.08.2019, 17.08.2019, 7.01.2020, 3.06.2020

if nargin < 1
    path = '.';
end

abstol = 0.0005;

%% Loading

if isfolder(path)

    cwd = pwd;
    cd(path);
    
    dim = load('ord.txt');
    time_step = load('stp.txt');
    
    obj = fread(fopen('obj.bin'), dim, 'double');    
    drg = fread(fopen('drg.bin'), dim, 'double');
    sys = fread(fopen('sys.bin'), [dim, dim], 'double');
    ela = fread(fopen('ela.bin'), [dim, dim], 'double');  %elasticity matrix
    mob = fread(fopen('mob.bin'), [dim, dim], 'double');  %projection matrix
    con = fread(fopen('con.bin'), [dim, dim], 'double');  %preconditionner
    %pts = fread(fopen('pts.bin'), dim, 'double');
    rhs = fread(fopen('rhs.bin'), dim, 'double');
    sol = fread(fopen('sol.bin'), dim, 'double');
    
    cd(cwd);
else
    error(['cannot find dump directory ',path]);
end

fprintf(1, 'loaded cytosim system of size %i with time_step %f\n', dim, time_step);

%% Check matrix

%figure('name', 'System matrix'); imshow(abs(sys)); 
%imshow(abs(mob)); set(gcf, 'name','Projection matrix');
%imshow(abs(ela)); set(gcf, 'name','Elasticity matrix');

if ( 1 )
    
    mat = eye(dim) - time_step * mob * ela;
    err0 = max(max(abs(mat-sys)));
    
    nbo = 0;
    nbv = 0;
    PRC = zeros(dim);
    for o = 0:max(max(obj))
        i = find(obj==o);
        if ~isempty(i)
            nbo = nbo + 1;
            nbv = nbv + length(i)^2;
            PRC(i,i) = inv(mat(i,i));
        end
    end
    
    fprintf(1, '%i mecables with %i block scalars\n', nbo, nbv);
    
    fprintf(1, '    norm8(matrix - reconstituted_matrix) : %e\n', err0);
    if ( err0 > 1e-8 )
        imshow(abs(mat));
        set(gcf, 'name', 'Reconstituted matrix');
    end

    fprintf(1, '    norm8(preconditionner - reconstituted_preconditionner) : %e\n', max(max(abs(PRC-con))));
    fprintf(1, '    norm(rhs) = %f\n', norm(rhs));
    fprintf(1, '    norm(sol) = %f\n', norm(sol));
    fprintf(1, '    norm(con*rhs) = %f\n', norm(con*rhs));

    %figure('name', 'Mecable footprint'); imshow(blk);

end

if ( 0 )
    mag = time_step * max(max(abs(mob * ela)));
    fprintf(1, 'norm8(system matrix - eye) : %e\n', mag);
end
if ( 0 )
    figure('Position', [50 50 1000 1000]);
    plot(reshape(mat,1,dim*dim), reshape(sys,1,dim*dim), '.')
    xl = xlim;
    ylim(xl);
    xlabel('Reconstituted matrix');
    ylabel('cytosim matrix');
end
if ( 0 )
    figure; hold on;
    plot(abs(sys), '^b');
    plot(abs(mat), 'vr');
end
if ( 0 )
    figure('name', 'System matrix structure');
    spy(mat)
    drawnow;
end

%% TESTING DIFFERENT ITERATIVE SOLVERS

tol = abstol / norm(rhs);

maxit = dim;
sss = sparse(sys);

solution = bicgstab(sss, rhs, abstol*0.001, maxit);
     
% BCGS
[vec,~,~,itr,rv0] = bicgstab(sss, rhs, tol, maxit);

report('BCGS', 2*itr, vec, rv0(end));

if 0
    figure;
    plot(vec, sol, 'k.');
    xlabel('matlab solution');
    ylabel('cytosim solution');
    xl = xlim;
    ylim(xl);
    drawnow;
end

figure('Name', 'Convergence');
convergence_axes = [];
convergence_plot(rv0/rv0(1), 'k');
xlabel('Mat*vec operations');
ylabel('Relative residual');
title('Solver Convergence');
hold on;

if 0 %% without preconditionning
 
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = tol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        [vec,rv0,itr] = cgstab(sss, rhs, [], [], OPT);
        convergence_plot(rv0(:,2), rv0(:,1)/rv0(1,1),'m');
        report(sprintf('BiCGS(%i)', OPT.ell), itr, vec, rv0(end,1));
    end
if 0
    % GMRES
    for i = 2:6
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit);
        convergence_plot(rv0/rv0(1),'g');
        report(sprintf('GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
end
    % IDRS-STAB
    for i = 0:4
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = IDRstabg5(sss, rhs, RS, tol/2, maxit, [], [], []);
        convergence_plot(rv0/rv0(1),'m');
        report(sprintf('IDRSTAB %03i', RS), itr, vec, rv0(end));
    end

    % IDRS
    OPT.smoothing = 0;
    for i = 0:4
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, [], [], [], OPT);
        convergence_plot(rv0/rv0(1),'r');
        report(sprintf('IDRS %03i', RS), itr, vec, rv0(end));
    end
end

    fprintf(2, 'Elasticity          has %i non-zero elements\n', nnz(ela));
    fprintf(2, 'Mobility/Projection has %i non-zero elements\n', nnz(mob));
    fprintf(2, 'Block conditionner  has %i non-zero elements\n', nnz(con));

if 1 %% WITH TRADITIONAL PRECONDITIONNING

    [vec,~,~,itr,rv0] = bicgstab(sss, rhs, tol, maxit, @precondition);
    report('P BCGS', 2*itr, vec, rv0(end));
    convergence_plot(rv0/rv0(1),'k');

end
if 0
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = tol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        OPT.TypePrecond = 'right';
        [vec,rv0,itr] = cgstab(sss, rhs, OPT, @precondition);
        convergence_plot(rv0(:,2), rv0(:,1)/rv0(1, 1),'m');
        report(sprintf('P BiCGS(%i)', OPT.ell), itr, vec, rv0(end));
    end
    % GMRES
    for i = 1:6
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = gmres(sss, rhs, RS, tol/3, maxit, [], @precondition);
        convergence_plot(rv0/rv0(1),'g');
        report(sprintf('P GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    % IDRS-STAB
    for i = 0:5
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = IDRstabg5(sss, rhs, RS, tol/2, maxit, @precondition, [], []);
        convergence_plot(rv0/rv0(1),'m');
        report(sprintf('P IDRSTAB %03i', RS), itr, vec, rv0(end));
    end
    % IDRS 
    OPT.smoothing = 1;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs_f(sss, rhs, RS, tol, maxit, @precondition, [], [], OPT);
        convergence_plot(rv0/rv0(1),'r-');
        report(sprintf('P IDRS %03i', RS), itr, vec, rv0(end));
    end
end


DRY = eye(dim) - diag(time_step./drg) * ela;
sDRY = sparse(DRY);

fprintf(2, 'DRY   has %i non-zero elements\n', nnz(sDRY));

if 0
    e = eig(DRY);
    figure('name', 'eigenvalues');
    plot(e, 'x');
    figure;
end
if ( 1 )
    figure('name', 'System matrix structure');
    spy(sDRY)
    drawnow;
    figure;
end

%% check ILU of matrix without projection:
 
if 1

    OPT.type = 'nofill';
    OPT.milu = 'off';
    [L, U] = ilu(sDRY, OPT);
    fprintf(2, 'incomplete LU      has %i + %i non-zero elements\n', nnz(L), nnz(U));
    
    [vec,~,~,itr,rv0] = bicgstab(sss, rhs, tol, maxit, L, U);
    convergence_plot(rv0/rv0(1),'b--');
    report('i BCGS', 2*itr, vec, rv0(end));
    if 0
        % IDRS-STAB
        for i = 0:5
            RS = min(2^i, dim);
            [vec,~,~,itr,rv0] = IDRstabg5(sss, rhs, RS, tol/2, maxit, L, U, []);
            convergence_plot(rv0/rv0(1),'m');
            report(sprintf('i IDRSTAB %03i', RS), itr, vec, rv0(end));
        end
    end
    OPT.smoothing = 1;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, L, U, [], OPT);
        convergence_plot(rv0/rv0(1), 'r--');
        report(sprintf('i IDRS %03i', RS), itr, vec, rv0(end));
    end
end

if 1

    S = ( DRY - DRY' );
    imshow(S);
    
    fprintf(2, 'Symmetrized LU has 2 * %i non-zero elements\n', nnz(S));
    [vec,~,~,itr,rv0] = bicgstab(sss, rhs, tol, maxit, S, S');
    convergence_plot(rv0/rv0(1),'g--');
    report('s BCGS', 2*itr, vec, rv0(end));
    
end

 %% incomplete Cholesky
if 1 

    clear OPT;
    OPT.michol = 'off';
    OPT.type = 'nofill';
    %OPT.diagcomp = 1;

    L = ichol(sDRY, OPT);
    U = L';
    fprintf(2, 'incomplete Cholesky has %i + %i non-zero elements\n', nnz(L), nnz(U));
    
    [vec,~,~,itr,rv0] = bicgstab(sss, rhs, tol, maxit, L, U);
    convergence_plot(rv0/rv0(1),'b--');
    report('i BCGS', 2*itr, vec, rv0(end));
    if 0
    % IDRS-STAB
    for i = 0:5
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = IDRstabg5(sss, rhs, RS, tol/2, maxit, L, U, []);
        convergence_plot(rv0/rv0(1),'m');
        report(sprintf('i IDRSTAB %03i', RS), itr, vec, rv0(end));
    end
    end
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, L, U, [], OPT);
        convergence_plot(rv0/rv0(1),'r--');
        report(sprintf('i IDRS %03i', RS), itr, vec, rv0(end));
    end
end



if 0 %% check different preconditionners

    % preconditionner = incomplete LU factorization
    [L, U] = ilu(sss, OPT);
    fprintf(2, 'incomplete LU      has %i + %i non-zero elements\n', nnz(L), nnz(U));

    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit, L, U);
    convergence_plot(rv0/rv0(1),'b--');
    report('i BCGS', 2*itr, vec, rv0(end));

    for i = 1:5
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, L, U);
        convergence_plot(rv0/rv0(1),'k--');
        report(sprintf('i GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    
    for i = 1:3
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, L, U, [], OPT);
        convergence_plot(rv0/rv0(1),'r--');
        report(sprintf('i IDRS %03i', RS), itr, vec, rv0(end));
    end

    
    if ( 0 )
        % preconditionner = Symmetric successive over-relaxation
        L = tril(sss);
        U = triu(sss);
        V = diag(sss);
        D = diag(V);
        M = (D+L)*diag(1./V)*(D+U);
        [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit, M);
        convergence_plot(rv0/rv0(1),'b--');
        report('s BCGS', 2*itr, vec);

        for i = 2:7
            RS = min(2^i, dim);
            [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, M);
            convergence_plot(rv0/rv0(1),'k--');
            report(sprintf('s GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
        end
    end
  
end

 %%
 
    function y = precondition(x)
        y = con * x;
    end

    function y = precondition2(x, mode)
        if ( strcmp(mode, 'notransp') )
            y = con * x;
        else
            y = con' * x;
        end
    end

    function convergence_plot(data, col)
        if isempty(convergence_axes)
            semilogy(data, col, 'Linewidth', 2);
            convergence_axes = gca;
        else
            semilogy(convergence_axes, data, col, 'Linewidth', 2);
        end
    end

    function report(s, i, v, r)
        tr = norm(sss*v-rhs);
        fprintf(1, '%-14s     converged after %4i vecmuls residual %f %f error %e\n', s, i, tr, r, norm(v-solution));
    end

end
