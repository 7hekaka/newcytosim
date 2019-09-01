function [sys, rhs, con] = cytosim_dump(path)

% This is used to explore Cytosim's linear system in matlab
% - load the matrices and vector from Cytosim's dump
% - plot convergence pattern of BICGstab, with and without preconditionning
%
% F. Nedelec, 16.10.2014, 03.2018, 06.2018, 26.01.2019, 30.06.2019,
% 11.08.2019, 17.08.2019

if nargin < 1
    path = '.';
end

abstol = 0.01;

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
    pts = fread(fopen('pts.bin'), dim, 'double');
    rhs = fread(fopen('rhs.bin'), dim, 'double');
    sol = fread(fopen('sol.bin'), dim, 'double');
    
    cd(cwd);
else
    error(['cannot find dump directory ',path]);
end

    fprintf(1, 'loaded system of size %i with time_step %f\n', dim, time_step);

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
    BLK = zeros(dim);
    for o = 0:max(max(obj))
        i = find(obj==o);
        if ~isempty(i)
            nbo = nbo + 1;
            nbv = nbv + length(i)^2;
            PRC(i,i) = mat(i,i);
            PRC(i,i) = inv(mat(i,i));
        end
    end
    
    fprintf(1, '%i mecables with %i block scalars\n', nbo, nbv);
    
    fprintf(1, '    norm8(system matrix - reconstituted matrix) : %e\n', err0);
    if ( err0 > 1e-8 )
        imshow(abs(mat));
        set(gcf, 'name', 'Reconstituted matrix');
    end

    fprintf(1, '    norm8(preconditionner - reconstituted preconditionner) : %e\n', max(max(abs(PRC-con))));
    fprintf(1, '    norm(rhs) = %f\n', norm(rhs));
    fprintf(1, '    norm(sol) = %f\n', norm(sol));
    fprintf(1, '    norm(con*rhs) = %f\n', norm(con*rhs));

    %figure('name', 'Mecable footprints'); imshow(blk);

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
   
tol = abstol / norm(rhs);

maxit = dim;
sss = sparse(sys);

if 1 %% without preconditionning
     
    % BCGS
    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit);
    
    if 0
        figure;
        plot(vec, sol, 'k.');
        xlabel('matlab solution');
        ylabel('cytosim solution');
        xl = xlim;
        ylim(xl);
        drawnow;
    end
    report('BCGS', 2*itr, vec, rv0(end));

    figure;
    semilogy(rv0/rv0(1),'k', 'Linewidth', 2);
    xlabel('Mat*vec operations');
    ylabel('Relative residual');
    title('Solver Convergence');
    hold on;
        
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = tol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        [vec,rv0,itr] = cgstab(sss, rhs, [], [], OPT);
        semilogy(rv0(:,2), rv0(:,1)/rv0(1,1),'m', 'Linewidth', 1);
        report(sprintf('BiCGS(%i)', OPT.ell), itr, vec, rv0(end,1));
    end

    % GMRES
    for i = 2:6
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit);
        semilogy(rv0/rv0(1),'g', 'Linewidth', 2);
        report(sprintf('GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    
    % IDRS
    OPT.smoothing = 0;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, [], [], [], OPT);
        semilogy(rv0/rv0(1),'r', 'Linewidth', 2);
        report(sprintf('IDRS %03i', RS), itr, vec, rv0(end));
    end

    % IDRS-STAB
    for i = 0:3
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = IDRstabg5(sss, rhs, RS, tol, maxit, [], [], []);
        semilogy(rv0/rv0(1),'m', 'Linewidth', 2);
        report(sprintf('IDRSTAB %03i', RS), itr, vec, rv0(end));
    end
end

    
    fprintf(2, 'Elasticity         has %i non-zero elements\n', nnz(ela));
    fprintf(2, 'Mobility           has %i non-zero elements\n', nnz(mob));
    fprintf(2, 'Block conditionner has %i non-zero elements\n', nnz(con));

    function y = solve(x)
        y = con * x;
    end

    function y = mfun2(x, mode)
        if ( strcmp(mode, 'notransp') )
            y = con * x;
        else
            y = con' * x;
        end
    end

if 0 %% WITH PRECONDITIONNING

    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit, @solve);
    report('P BCGS', 2*itr, vec, rv0(end));

    figure;
    semilogy(rv0/rv0(1),'k', 'Linewidth', 2);
    xlabel('Mat*vec operations');
    ylabel('Relative residual');
    title('Convergence with Preconditionning');
    hold on;
    
    if ( 0 )
        % BiCGStab(L)
        for i = 1:4
            OPT.Tol = tol;
            OPT.ell = min(2^i, dim);
            OPT.MaxIt = maxit;
            OPT.TypePrecond = 'right';
            [vec,rv0,itr] = cgstab(sss, rhs, OPT, @solve);
            semilogy(rv0(:,2), rv0(:,1)/rv0(1, 1),'m', 'Linewidth', 1);
            report(sprintf('P BiCGS(%i)', OPT.ell), itr, vec, rv0(end));
        end
    end
    
    for i = 1:6
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol/3, maxit, [], @solve);
        semilogy(rv0/rv0(1),'g', 'Linewidth', 2);
        report(sprintf('P GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    
    % IDRS 
    OPT.smoothing = 1;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, @solve, [], [], OPT);
        semilogy(rv0/rv0(1),'r-', 'Linewidth', 2);
        report(sprintf('P IDRS %03i', RS), itr, vec, rv0(end));
    end
  
    % IDRS-STAB
    for i = 0:3
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = IDRstabg5(sss, rhs, RS, tol/3, maxit, @solve, [], []);
        semilogy(rv0/rv0(1),'m', 'Linewidth', 2);
        report(sprintf('P IDRSTAB %03i', RS), itr, vec, rv0(end));
    end

end
    
if 0 %% check different preconditionners

    % preconditionner = incomplete LU factorization
    [L, U] = ilu(sss);
    fprintf(2, 'incomplete LU      has %i + %i non-zero elements\n', nnz(L), nnz(U));

    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit, L, U);
    semilogy(rv0/rv0(1),'b--', 'Linewidth', 2);
    report('i BCGS', 2*itr, vec, rv0(end));

    for i = 1:5
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, L, U);
        semilogy(rv0/rv0(1),'k--', 'Linewidth', 2);
        report(sprintf('i GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    
    for i = 1:3
        RS = min(2^i, dim);
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, L, U, [], OPT);
        semilogy(rv0/rv0(1),'r--', 'Linewidth', 2);
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
        semilogy(rv0/rv0(1),'b--', 'Linewidth', 2);
        report('s BCGS', 2*itr, vec);

        for i = 2:7
            RS = min(2^i, dim);
            [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, M);
            semilogy(rv0/rv0(1),'k--', 'Linewidth', 2);
            report(sprintf('s GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
        end
    end
  
 end

    function report(s, i, v, r)
        tr = norm(sss*v-rhs);
        fprintf(1, '%-14s     converged after %4i vecmuls residual %f %f error %e\n', s, i, tr, r, norm(v-sol));
    end

end
