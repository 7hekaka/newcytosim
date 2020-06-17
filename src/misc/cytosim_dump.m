function [sys, rhs, con] = cytosim_dump(path)

% This is used to explore Cytosim's linear system in matlab
% - load the matrices and vector from Cytosim's dump
% - plot convergence pattern of BICGstab, with and without preconditionning
%
% F. Nedelec, 16.10.2014, 03.2018, 06.2018, 26.01.2019, 30.06.2019,
% 11.08.2019, 17.08.2019, 7.01.2020, 3.06.2020, 15.06.2002, 17.06.2020

if nargin < 1
    path = '.';
end

mulcnt = 0;
abstol = 0.005;

%% Loading

if isfolder(path)

    cwd = pwd;
    cd(path);
    
    ord = load('ord.txt');
    dim = ord(1);
    time_step = load('stp.txt');
    
    obj = fread(fopen('obj.bin'), dim, 'double');    
    drg = fread(fopen('drg.bin'), dim, 'double');
    sys = fread(fopen('sys.bin'), [dim, dim], 'double');
    ela = fread(fopen('ela.bin'), [dim, dim], 'double');  % elasticity matrix
    mob = fread(fopen('mob.bin'), [dim, dim], 'double');  % projection matrix
    con = fread(fopen('con.bin'), [dim, dim], 'double');  % preconditionner
    %pts = fread(fopen('pts.bin'), dim, 'double');
    rhs = fread(fopen('rhs.bin'), dim, 'double');
    sol = fread(fopen('sol.bin'), dim, 'double');
    
    cd(cwd);
else
    error(['cannot find dump directory ',path]);
end

fprintf(1, '----------------------- loaded system of size %i with time_step %f -----------------------\n', dim, time_step);

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
    fprintf(2, '    norm(ela-transpose(ela)) = %f\n', norm(ela-ela',1));

    fprintf(2, '    norm8(matrix - reconstituted_matrix) : %e\n', err0);
    if ( err0 > 1e-8 )
        imshow(abs(mat));
        set(gcf, 'name', 'Reconstituted matrix');
    end

    fprintf(2, '    norm8(preconditionner - reconstituted_preconditionner) : %e\n', norm(PRC-con,1));
    fprintf(2, '    norm(rhs) = %f', norm(rhs));
    fprintf(2, '    norm(sol) = %f', norm(sol));
    fprintf(2, '    norm(con*rhs) = %f\n', norm(con*rhs));

    %figure('name', 'Mecable footprint'); imshow(PRC);

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
    figure('name', 'System matrix values');
    plot(abs(sys), '^b'); hold on;
    plot(abs(mat), 'vr');
end
if ( 0 )
    figure('name', 'System matrix structure');
    spy(mat)
    drawnow;
end

%% RECALCULATE SOLUTION

reltol = abstol / norm(rhs);

maxit = dim;
system = sparse(sys);

solution = bicgstab(@multiply, rhs, reltol*0.001, maxit);

fprintf(2, '    norm(sol - matlab_sol) = %f\n', norm(sol-solution));
    
if 0
    figure('Name', 'Validation of solution');;
    plot(solution, sol, 'k.');
    xlabel('matlab solution');
    ylabel('cytosim solution');
    xl = xlim;
    ylim(xl);
    drawnow;
end

%% BCGS

mulcnt = 0;
[vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit);

report('bicgstab', mulcnt, vec, rv0(end));

figure('Name', 'Convergence');
convergence_axes = [];
convergence_plot(rv0/rv0(1), 'k');
xlabel('Mat*vec operations');
ylabel('Relative residual');
title('Solver Convergence');
hold on;

% Matlab BiCGStab(L)
mulcnt = 0;
[vec,~,~,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit);
convergence_plot(rv0/rv0(1), 'b');
report('bicgstab(1)', mulcnt, vec, rv0(end));

if 0 %% without preconditionning
 
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = reltol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        [vec,rv0,itr] = cgstab(@multiply, rhs, [], [], OPT);
        convergence_plot(rv0(:,2), rv0(:,1)/rv0(1,1),'m');
        report(sprintf('BiCGS(%i)', OPT.ell), itr, vec, rv0(end,1));
    end
if 0
    % GMRES
    for i = 2:6
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = gmres(@multiply, rhs, RS, reltol, maxit);
        convergence_plot(rv0/rv0(1),'g');
        report(sprintf('GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
end
    % IDRS-STAB
    for i = 0:4
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = IDRstabg5(@multiply, rhs, RS, reltol/2, maxit, [], [], []);
        convergence_plot(rv0/rv0(1),'m');
        report(sprintf('IDRSTAB %03i', RS), itr, vec, rv0(end));
    end

    % IDRS
    OPT.smoothing = 0;
    for i = 0:4
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs(@multiply, rhs, RS, reltol, maxit, [], [], [], OPT);
        convergence_plot(rv0/rv0(1),'r');
        report(sprintf('IDRS %03i', RS), itr, vec, rv0(end));
    end
end

 %% USING PRECONDITIONNER CALCULATED BY CYTOSIM
 
 fprintf(2, '    Elasticity            has %9i elements\n', nnz(ela));
 fprintf(2, '    Mobility/Projection   has %9i elements\n', nnz(mob));
 fprintf(2, '--- Block preconditionner has %9i elements\n', nnz(con));

 if 1
     mulcnt = 0;
     [vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @precondition);
     report('P bicgstab', mulcnt, vec, rv0(end));
     convergence_plot(rv0/rv0(1),'k--');
     
     mulcnt = 0;
     % Matlab BiCGStab(L)
     [vec,~,~,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, @precondition);
     convergence_plot(rv0/rv0(1), 'b--');
     report('P bicgstab(1)', mulcnt, vec, rv0(end));
     
     % checking the reconstituted block preconditionner:
     mulcnt = 0;
     [vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @preconditionPRC);
     report('R bicgstab', mulcnt, vec, rv0(end));
     convergence_plot(rv0/rv0(1),'k--');
end

if 0
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = reltol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        OPT.TypePrecond = 'right';
        [vec,rv0,itr] = cgstab(@multiply, rhs, OPT, @precondition);
        convergence_plot(rv0(:,2), rv0(:,1)/rv0(1, 1),'m');
        report(sprintf('P BiCGS(%i)', OPT.ell), itr, vec, rv0(end));
    end

    % GMRES
    for i = 1:6
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = gmres(@multiply, rhs, RS, reltol/3, maxit, [], @precondition);
        convergence_plot(rv0/rv0(1),'g');
        report(sprintf('P GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    % IDRS-STAB
    for i = 0:5
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = IDRstabg5(@multiply, rhs, RS, reltol/2, maxit, @precondition, [], []);
        convergence_plot(rv0/rv0(1),'m');
        report(sprintf('P IDRSTAB %03i', RS), itr, vec, rv0(end));
    end
    % IDRS 
    OPT.smoothing = 1;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs_f(@multiply, rhs, RS, reltol, maxit, @precondition, [], [], OPT);
        convergence_plot(rv0/rv0(1),'r-');
        report(sprintf('P IDRS %03i', RS), itr, vec, rv0(end));
    end
end

%% a possible sparse symmetric preconditionner
% we average all point drag coefficient to derive a matrix that is
% symmetric and ammenable to incomplete Cholesky factorization

if 1
    
    val = time_step / mean(1./drg);
    DRY = eye(dim) - diag(val.*ones(length(ela), 1)) * ela;
    DRY = sparse(DRY);
    
    fprintf(2, '--- DRY is symmetric with %i elements; ', nnz(DRY));
    
    clear OPT;
    OPT.michol = 'off';
    OPT.type = 'nofill';
    L = ichol(DRY, OPT);
    fprintf(2, 'incomplete Cholesky has %i elements\n', nnz(L));
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    convergence_plot(rv0/rv0(1),'k:');
    report('y bicgstab', mulcnt, vec, rv0(end));
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, L, L');
    convergence_plot(rv0/rv0(1),'b:');
    report('y bicgstab(1)', mulcnt, vec, rv0(end));
end

%% a smaller sparse symmetric preconditionner
% We project the previous preconditionner to make it isotropic in X, Y Z

if 1
    
    val = time_step / mean(1./drg);
    SML = eye(dim) - diag(val.*ones(length(ela), 1)) * ela;

    if ( ord(2) == 3 )
        SML = SML(1:3:end, 1:3:end) + SML(2:3:end, 2:3:end) + SML(3:3:end, 3:3:end);
        SML = (1/3) * sparse(SML);
    elseif ( ord(2) == 2 )
        SML = SML(1:2:end, 1:2:end) + SML(2:2:end, 2:2:end);
        SML = (1/2) * sparse(SML);
    end
    
    fprintf(2, '--- SML is symmetric with %i elements; ', nnz(SML));
    
    clear OPT;
    OPT.michol = 'off';
    OPT.type = 'nofill';
    SML = ichol(SML, OPT);
    fprintf(2, 'incomplete Cholesky has %i elements\n', nnz(SML));
    
    if ( ord(2) == 3 )
        L = sparse(dim, dim);
        L(1:3:end, 1:3:end) = SML;
        L(2:3:end, 2:3:end) = SML;
        L(3:3:end, 3:3:end) = SML;
    elseif ( ord(2) == 2 )
        L = sparse(dim, dim);
        L(1:2:end, 1:2:end) = SML;
        L(2:2:end, 2:2:end) = SML;
    end
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    convergence_plot(rv0/rv0(1),'k:');
    report('s bicgstab', mulcnt, vec, rv0(end));
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, L, L');
    convergence_plot(rv0/rv0(1),'b:');
    report('s bicgstab(1)', mulcnt, vec, rv0(end));
end


%% incomplete LU factorization

% the WET matrix is not necessarily symmetric because of the diagonal matrix
% may not commute with 'ela'. However, if all point drags are equal, then
% surely WET is symmetric and we can use incomplete Cholesky factorization

WET = eye(dim) - diag(time_step./drg) * ela;

sWET = sparse(WET);

fprintf(2, '--- WET has %i elements; ', nnz(sWET));
fprintf(2, 'norm(WET-transpose(WET)) = %.2f; ', norm(WET-WET',1));

if 0 && ( norm(WET-WET',1) < 1 )
    
    clear OPT;
    OPT.michol = 'off';
    OPT.type = 'nofill';
    L = ichol(sWET, OPT);
    fprintf(2, 'incomplete Cholesky has %i elements\n', nnz(L));
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    convergence_plot(rv0/rv0(1),'k.-');
    report('w bicgstab', mulcnt, vec, rv0(end));
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, L, L');
    convergence_plot(rv0/rv0(1),'b.-');
    report('w bicgstab(1)', mulcnt, vec, rv0(end));
    
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs(system, rhs, RS, reltol, maxit, L, L', [], OPT);
        convergence_plot(rv0/rv0(1),'r.-');
        report(sprintf('w IDRS %03i', RS), 2*itr, vec, rv0(end));
    end

else

    clear OPT;
    OPT.type = 'nofill';
    OPT.milu = 'off';
    [L, U] = ilu(sWET, OPT);
    fprintf(2, '    incomplete LU has %i + %i elements\n', nnz(L), nnz(U));
    
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, U);
    report('i bicgstab', mulcnt, vec, rv0(end));
    convergence_plot(rv0/rv0(1),'k-.');
    
    % Matlab BiCGStab(L)
    mulcnt = 0;
    [vec,~,~,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, L, U);
    convergence_plot(rv0/rv0(1), 'b-.');
    report('i bicgstab(1)', mulcnt, vec, rv0(end));
    
end

if 0
    % IDRS-STAB
    for i = 0:5
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = IDRstabg5(@multiply, rhs, RS, reltol/2, maxit, L, U, []);
        convergence_plot(rv0/rv0(1),'m');
        report(sprintf('i IDRSTAB %03i', RS), itr, vec, rv0(end));
    end
    
    OPT.smoothing = 1;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs(@multiply, rhs, RS, reltol, maxit, L, U, [], OPT);
        convergence_plot(rv0/rv0(1), 'r--');
        report(sprintf('i IDRS %03i', RS), itr, vec, rv0(end));
    end
end


if 0 %% check more preconditionners

    % preconditionner = incomplete LU factorization
    [L, U] = ilu(system, OPT);
    fprintf(2, 'incomplete LU      has %i + %i elements\n', nnz(L), nnz(U));

    [vec,~,~,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, U);
    convergence_plot(rv0/rv0(1),'b--');
    report('i bicgstab', 2*itr, vec, rv0(end));

    for i = 1:5
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = gmres(@multiply, rhs, RS, reltol, maxit, L, U);
        convergence_plot(rv0/rv0(1),'k--');
        report(sprintf('i GMRES %03i', RS), itr(1)*RS+itr(2), vec, rv0(end));
    end
    
    for i = 1:3
        RS = min(2^i, dim);
        [vec,~,~,itr,rv0] = idrs(@multiply, rhs, RS, reltol, maxit, L, U, [], OPT);
        convergence_plot(rv0/rv0(1),'r--');
        report(sprintf('i IDRS %03i', RS), itr, vec, rv0(end));
    end
  
end

 %%
 
    function y = multiply(x)
        mulcnt = mulcnt + 1;
        y = system * x;
    end

    function y = precondition(x)
        y = con * x;
    end

    function y = preconditionPRC(x)
        y = PRC * x;
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
        tr = norm(system*v-rhs);
        fprintf(1, '    %-14s     converged after %4i vecmuls residual %f %f error %e\n', s, i, tr, r, norm(v-solution));
    end

end
