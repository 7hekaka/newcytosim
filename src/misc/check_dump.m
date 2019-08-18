function check_dump(path)

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
    
    ord = load('ord.txt');
    time_step = load('stp.txt');
    
    obj = fread(fopen('obj.bin'), ord, 'double');    
    drg = fread(fopen('drg.bin'), ord, 'double');
    sys = fread(fopen('sys.bin'), [ord, ord], 'double');
    ela = fread(fopen('ela.bin'), [ord, ord], 'double');  %elasticity matrix
    mob = fread(fopen('mob.bin'), [ord, ord], 'double');  %projection matrix
    con = fread(fopen('con.bin'), [ord, ord], 'double');  %preconditionner
    pts = fread(fopen('pts.bin'), ord, 'double');
    rhs = fread(fopen('rhs.bin'), ord, 'double');
    sol = fread(fopen('sol.bin'), ord, 'double');
    
    cd(cwd);
else
    error(['cannot find dump directory ',path]);
end

    fprintf(1, 'loaded system of size %i with time_step %f\n', ord, time_step);

%% Check matrix

%figure('name', 'System matrix'); imshow(abs(sys)); 
%imshow(abs(mob)); set(gcf, 'name','Projection matrix');
%imshow(abs(ela)); set(gcf, 'name','Elasticity matrix');

if ( 1 )
    
    mat = eye(ord) - time_step * mob * ela;
    err0 = max(max(abs(mat-sys)));
       
    nbo = 0;
    nbv = 0;
    PRC = zeros(ord);
    BLK = zeros(ord);
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
    plot(reshape(mat,1,ord*ord), reshape(sys,1,ord*ord), '.')
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

if 1 %% without preconditionning
    
    tol = abstol / norm(rhs);

    maxit = ord;
    sss = sparse(sys);
    
    % BCGS
    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit);
    
    if 1
        figure;
        plot(vec, sol, 'k.');
        xlabel('matlab solution');
        ylabel('cytosim solution');
        xl = xlim;
        ylim(xl);
        drawnow;
    end
    report('BCGS', 2*itr, vec);
    
    figure;
    semilogy(rv0/rv0(1),'k', 'Linewidth', 2);
    xlabel('Mat*vec operations');
    ylabel('Relative residual');
    title('Solver Convergence');
    hold on;

    if ( 0 )
        % QMR
        [vec,flg,res,itr,rv0] = qmr(sss, rhs, tol, maxit);
        semilogy(rv0/rv0(1),'m:', 'Linewidth', 2);
        report('QMR', 2*itr, vec);
    end
    
    % GMRES
    for i = 4:6
        RS = 2^i;
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit);
        semilogy(rv0/rv0(1),'g', 'Linewidth', 2);
        report(sprintf('GMRES %03i', RS), itr(1)*RS+itr(2), vec);
    end
    
    % IDRS
    OPT.smoothing = 0;
    for i = 0:3
        RS = 2^i;
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, [], [], [], OPT);
        semilogy(rv0/rv0(1),'r', 'Linewidth', 2);
        report(sprintf('IDRS %03i', RS), itr, vec);
    end

    % IDRS-STAB
    for i = 0:3
        RS = 2^i;
        [vec,flg,res,itr] = IDRstabg4(sss, rhs, RS, tol, maxit, [], [], []);
        semilogy(rv0/rv0(1),'m', 'Linewidth', 2);
        report(sprintf('IDRSTAB %03i', RS), itr, vec);
    end
end

    
    fprintf(2, 'Elasticity         has %i non-zero elements\n', nnz(ela));
    fprintf(2, 'Mobility           has %i non-zero elements\n', nnz(mob));
    fprintf(2, 'Block conditionner has %i non-zero elements\n', nnz(con));

    function y = mfun1(x)
        y = con * x;
    end

    function y = mfun2(x, mode)
        if ( strcmp(mode, 'notransp') )
            y = con * x;
        else
            y = con' * x;
        end
    end

if 1 %% WITH PRECONDITIONNING

    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit, @mfun1);
    report('P BCGS', 2*itr, vec);

    figure;
    semilogy(rv0/rv0(1),'k', 'Linewidth', 2);
    xlabel('Mat*vec operations');
    ylabel('Relative residual');
    title('Convergence with Preconditionning');
    hold on;

    for i = 2:6
        RS = 2^i;
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol/5, maxit, [], @mfun1);
        semilogy(rv0/rv0(1),'g', 'Linewidth', 2);
        report(sprintf('P GMRES %03i', RS), itr(1)*RS+itr(2), vec);
    end
    
    
    if ( 0 )
        % QMR method
        [vec,flg,res,itr,rv0] = qmr(sss, rhs, tol, maxit, @mfun2);
        semilogy(rv0/rv0(1),'m-', 'Linewidth', 2);
        report('P QMR', 2*itr, vec);
    end
    
    
    % IDRS 
    OPT.smoothing = 1;

    for i = 0:3
        RS = 2^i;
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, @mfun1, [], [], OPT);
        semilogy(rv0/rv0(1),'r-', 'Linewidth', 2);
        report(sprintf('P IDRS %03i', RS), itr, vec);
    end    
  
    % IDRS-STAB
    for i = 0:3
        RS = 2^i;
        [vec,flg,res,itr] = IDRstabg4(sss, rhs, RS, tol/5, maxit, @mfun1, [], []);
        semilogy(rv0/rv0(1),'m', 'Linewidth', 2);
        report(sprintf('P IDRSTAB %03i', RS), itr, vec);
    end

end
    
if 0 %% check different preconditionners

    % preconditionner = incomplete LU factorization
    [L, U] = ilu(sss);
    fprintf(2, 'incomplete LU      has %i + %i non-zero elements\n', nnz(L), nnz(U));

    [vec,flg,res,itr,rv0] = bicgstab(sss, rhs, tol, maxit, L, U);
    semilogy(rv0/rv0(1),'b--', 'Linewidth', 2);
    report('i BCGS', 2*itr, vec);

    for i = 2:6
        RS = 2^i;
        [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, L, U);
        semilogy(rv0/rv0(1),'k--', 'Linewidth', 2);
        report(sprintf('i GMRES %03i', RS), itr(1)*RS+itr(2), vec);
    end
    
    for i = 1:3
        RS = 2^i;
        [vec,flg,res,itr,rv0] = idrs(sss, rhs, RS, tol, maxit, L, U, [], OPT);
        semilogy(rv0/rv0(1),'r--', 'Linewidth', 2);
        report(sprintf('i IDRS %03i', RS), itr, vec);
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
        fprintf(1, 'BCGS-SSOR     converged after %4i vecmuls %f\n', 2*itr, res);
        
        for i = 2:7
            RS = 2^i;
            [vec,flg,res,itr,rv0] = gmres(sss, rhs, RS, tol, maxit, M);
            semilogy(rv0/rv0(1),'k--', 'Linewidth', 2);
            fprintf(1, 'GMRES-SSOR %03i converged after %4i vecmuls %f\n', RS, (itr(1)-1)*RS+itr(2), res);
        end
    end
  
 end

    function report(s, i, v)
        fprintf(1, '%-14s     converged after %4i vecmuls residual %f error %e\n', s, i, norm(sss*v-rhs), norm(v-sol));
    end

end
