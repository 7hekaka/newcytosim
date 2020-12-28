function [sys, rhs, con] = cytosim_dump(path)

% This is used to explore Cytosim's linear system in matlab
% - load the matrices and vector from Cytosim's dump
% - plot convergence pattern of BICGstab, with and without preconditionning
%
% F. Nedelec, 16.10.2014, 03.2018, 06.2018, 26.01.2019, 30.06.2019,
% 11.08.2019, 17.08.2019, 7.01.2020, 03.06.2020, 19.06.2020, 27.12.2020

if nargin < 1
    path = '.';
end

abstol = 0.001;


%% Loading

if isfolder(path)

    cwd = pwd;
    cd(path);
    
    ord = load('ord.txt');
    dim = ord(1);
    if ord(3) == 4
        precision = 'single';
    else
        precision = 'double';
    end
    stp = load('stp.txt');
    time_step = stp(1);
    if ( length(stp) > 1 )
        abstol = stp(2);
    end
    obj = fread(fopen('obj.bin'), dim, precision);    
    drg = fread(fopen('drg.bin'), dim, precision);
    sys = fread(fopen('sys.bin'), [dim, dim], precision);
    ela = fread(fopen('ela.bin'), [dim, dim], precision);  % elasticity matrix
    mob = fread(fopen('mob.bin'), [dim, dim], precision);  % projection matrix
    con = fread(fopen('con.bin'), [dim, dim], precision);  % preconditionner
    %pts = fread(fopen('pts.bin'), dim, precision);
    rhs = fread(fopen('rhs.bin'), dim, precision);
    sol = fread(fopen('sol.bin'), dim, precision);
    
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
    err0 = norm(mat-sys, 1);
    
    nbo = 0;
    nbv = 0;
    for o = 0:max(max(obj))
        i = find(obj==o);
        if ~isempty(i)
            nbo = nbo + 1;
            nbv = nbv + length(i)^2;
        end
    end
    
    fprintf(1, '%i mecables with %i block scalars (%s)\n', nbo, nbv, precision);
    fprintf(2, '    norm(ela-transpose(ela)) = %f\n', norm(ela-ela',1));

    fprintf(2, '    norm8(matrix - reconstituted_matrix) : %e\n', err0);
    if ( err0 > 1e-5 )
        figure;
        imshow(abs(mat));
        set(gcf, 'name', 'Reconstituted matrix');
    end

    norm_rhs = norm(rhs);
    fprintf(2, '    norm(rhs) = %f', norm(rhs));
    fprintf(2, '    norm(sol) = %f', norm(sol));
    fprintf(2, '    norm(con*rhs) = %f\n', norm(con*rhs));

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
    uuu = eye(dim) - time_step * ela;
    figure('name', 'System structure', 'Position', [150 300 1800 600]);
    subplot(1,3,1);
    spy(uuu);
    subplot(1,3,2);
    rcm = symrcm(uuu);
    spy(uuu(rcm,rcm));
    title('Reverse Cuthill-McKee ordering');
    subplot(1,3,3);
    amd = symamd(uuu);
    spy(uuu(amd,amd));
    title('Approximate Minimum Degree permutation');
    drawnow;
end

%% RECALCULATE SOLUTION

reltol = abstol / norm(rhs);

maxit = dim;
system = sparse(sys);

mulcnt = 0;
solution = bicgstab(@multiply, rhs, reltol*0.001, maxit);

fprintf(2, '    %i matvecs; norm(sol - matlab_sol) = %f\n', mulcnt, norm(sol-solution));
    
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
convergence_axes = [];

mulcnt = 0;
[vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit);
report('bicgstab', mulcnt, vec, res, rv0, '-.');

if 0
    % Matlab BiCGStab(L)
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstabl(@multiply, rhs, reltol, maxit);
    report('bicgstab(1)', mulcnt, vec, res, rv0, ':');
end

if 0
    % Matlab CGS
    mulcnt = 0;
    [vec,~,res,itr,rv0] = cgs(@multiply, rhs, reltol, maxit);
    report('CGS', mulcnt, vec, res, rv0, ':');
end

if 0 %% no preconditionning
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = reltol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        [vec,rv0,itr] = cgstab(system, rhs, [], [], OPT);
        report(sprintf('BiCGS(%i)', OPT.ell), 2*itr, vec, rv0(end,1), rv0(:,1), ':');
    end
end
if 1
    % GMRES
    for i = 4:6
        mulcnt = 0;
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = gmres(@multiply, rhs, RS, reltol, maxit);
        report(sprintf('GMRES %03i', RS), mulcnt, vec, res, rv0, ':');
    end
end
if 0 %% no preconditionning
    % IDRS-STAB
    for i = 0:4
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = IDRstabg5(system, rhs, RS, reltol, maxit, [], [], []);
        report(sprintf('IDRSTAB %03i', RS), 2*itr, vec, res, rv0, ':');
    end
    % IDRS
    OPT.smoothing = 0;
    for i = 0:4
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = idrs(system, rhs, RS, reltol, maxit, [], [], [], OPT);
        report(sprintf('IDRS %03i', RS), 2*itr, vec, res, rv0, ':');
    end
end
if 1
    % CORS
    mulcnt = 0;
    [vec,~,res,itr,rv0] = CORS(@multiply, rhs, reltol, maxit);
    report('CORS', mulcnt, vec, res, rv0, ':');
    % BiCOR
    mulcnt = 0;
    [vec,~,res,itr,rv0] = BiCOR(@multiply, rhs, reltol, maxit);
    report('BiCOR', mulcnt, vec, res, rv0, ':');
    % BiCGCR2
    [vec,~,res,itr,rv0] = BiCGCR2(system, rhs, reltol, maxit);
    report('BiCGCR2', 2*itr, vec, res, rv0, ':');
end

drawnow;

%% USING LOADED PRECONDITIONNER CALCULATED BY CYTOSIM
fprintf(1, '  --  --  --  --  --  --  --  -- PRECONDITIONNED --  --  --  --  --  --  --  --  --\n');

% Calculate block-diagonal preconditionner
PRC = zeros(dim);
for o = 0:max(max(obj))
    i = find(obj==o);
    if ~isempty(i)
        PRC(i,i) = inv(mat(i,i));
    end
end

%figure('name', 'Mecable footprint'); imshow(PRC);
fprintf(2, '    norm8(preconditionner - reconstituted_preconditionner) : %e\n', norm(PRC-con,1));

fprintf(2, '    Elasticity            has %9i elements\n', nnz(ela));
fprintf(2, '    Mobility/Projection   has %9i elements\n', nnz(mob));
fprintf(2, '    Given Preconditionner has %9i elements\n', nnz(con));
fprintf(2, '    Block preconditionner has %9i elements\n', nnz(PRC));

convergence_axes = [];

mulcnt = 0;
[vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @precondition);
report('P bicgstab', mulcnt, vec, res, rv0);

% checking the reconstituted block preconditionner:
mulcnt = 0;
[vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @preconditionPRC);
report('R bicgstab', mulcnt, vec, res, rv0);

if 0
    mulcnt = 0;
    % Matlab BiCGStab(L)
    [vec,~,res,~,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, @precondition);
    report('P bicgstab(1)', mulcnt, vec, res, rv0);
end
if 0
    % GMRES
    for i = 2:6
        mulcnt = 0;
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = gmres(@multiply, rhs, RS, reltol/3, maxit, [], @precondition);
        str = sprintf('P GMRES %03i', RS);
        report(str, mulcnt, vec, res, rv0);
    end
end
if 1
    % CORS
    mulcnt = 0;
    [vec,~,res,itr,rv0] = CORS(@multiply, rhs, reltol, maxit, @preconditionPRC);
    report('R CORS', mulcnt, vec, res, rv0);
    % BiCOR
    mulcnt = 0;
    [vec,~,res,itr,rv0] = BiCOR(@multiply, rhs, reltol, maxit, @preconditionPRC);
    report('R BiCOR', mulcnt, vec, res, rv0);
    % BiCGCR2
    [vec,~,res,itr,rv0] = BiCGCR2(system, rhs, reltol, maxit, @preconditionPRC);
    report('R BiCGCR2', 2*itr, vec, res, rv0);
end

if 0
    % BiCGStab(L)
    for i = 1:4
        OPT.Tol = reltol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        OPT.TypePrecond = 'right';
        [vec,rv0,itr] = cgstab(@multiply, rhs, OPT, @precondition);
        report(sprintf('P BiCGS(%i)', OPT.ell), itr, vec, rv0(end), rv0);
    end
    % IDRS-STAB
    for i = 0:5
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = IDRstabg5(@multiply, rhs, RS, reltol/2, maxit, @precondition, [], []);
        report(sprintf('P IDRSTAB %03i', RS), itr, vec, res, rv0);
    end
    % IDRS
    OPT.smoothing = 1;
    for i = 0:3
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = idrs_f(@multiply, rhs, RS, reltol, maxit, @precondition, [], [], OPT);
        report(sprintf('P IDRSF %03i', RS), itr, vec, res, rv0);
    end
end

%% a possible sparse symmetric preconditionner
% we average all point drag coefficient to derive a matrix that is
% symmetric and ammenable to incomplete Cholesky factorization

fprintf(2, 'Mecable drag coefficients : %f +/- %f\n', mean(1./drg), var(1./drg));

if 0
    
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
    [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    report('y bicgstab', mulcnt, vec, res, rv0);
    
    try
        OPT.michol = 'off';
        OPT.type = 'ict';
        OPT.droptol = 0.1;
        L = ichol(DRY, OPT);
        fprintf(2, 'dropped incomplete Cholesky has %i elements\n', nnz(L));
        
        mulcnt = 0;
        [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
        report('y bicgstab', mulcnt, vec, res, rv0);
    catch
        fprintf(2, 'dropped incomplete Cholesky failed!\n');
    end
    
    L = chol(DRY, 'lower');
    fprintf(2, 'complete Cholesky has %i elements\n', nnz(L));
    
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    report('y bicgstab', mulcnt, vec, res, rv0);
end

%% a smaller sparse symmetric preconditionner
% We project the previous preconditionner to make it isotropic in X, Y Z

if 0
    
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
    [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    report('s bicgstab', mulcnt, vec, res, rv0);
    
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, L, L');
    report('s bicgstab(1)', mulcnt, vec, res, rv0);
end


%% incomplete LU factorization

% the WET matrix is not necessarily symmetric because of the diagonal matrix
% may not commute with 'ela'. However, if all point drags are equal, then
% surely WET is symmetric and we can use incomplete Cholesky factorization

if 0
    
    WET = eye(dim) - diag(time_step./drg) * ela;
    
    sWET = sparse(WET);
    
    fprintf(2, '--- WET has %i elements; ', nnz(sWET));
    fprintf(2, 'norm(WET-transpose(WET)) = %.2f; ', norm(WET-WET',1));
    
    if ( norm(WET-WET',1) < 1 )
        
        clear OPT;
        OPT.michol = 'off';
        OPT.type = 'nofill';
        L = ichol(sWET, OPT);
        fprintf(2, 'incomplete Cholesky has %i elements\n', nnz(L));
        
        mulcnt = 0;
        [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
        report('I bicgstab', mulcnt, vec, res, rv0);
        
    else
        
        clear OPT;
        OPT.type = 'nofill';
        OPT.milu = 'off';
        [L, U] = ilu(sWET, OPT);
        fprintf(2, '    incomplete LU has %i + %i elements\n', nnz(L), nnz(U));
        
        mulcnt = 0;
        [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, U);
        report('i bicgstab', mulcnt, vec, res, rv0);
        
    end   
end


%% check preconditionning with iLU on the full matrix

if 0
    
    % preconditionner = incomplete LU factorization
    [L, U] = ilu(system, OPT);
    fprintf(2, 'incomplete LU      has %i + %i elements\n', nnz(L), nnz(U));
    
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, U);
    report('iLU bicgstab', mulcnt, vec, res, rv0);

end

%% Functions

    function y = multiply(x, mode)
        mulcnt = mulcnt + 1;
        if nargin < 2 || strcmp(mode, 'notransp')
            y = system * x;
        else
            y = system' * x;
        end
    end

    function y = precondition(x, mode)
        if nargin < 2 || strcmp(mode, 'notransp')
            y = con * x;
        else
            y = con' * x;
        end
    end

    function y = preconditionPRC(x, mode)
        if nargin < 2 || strcmp(mode, 'notransp')
            y = PRC * x;
        else
            y = PRC' * x;
        end
    end


    function convergence_plot(str, mvs, data, lin)
        %fprintf(1, '%s    %4.1f %4i\n', txt, mvs, length(data));
        mvs = (0:(length(data)-1)) * ( mvs / length(data) );
        % crop data:
        up = min(256,length(data));
        dat = data(1:up);
        mvs = mvs(1:up);
        if isempty(convergence_axes)
            figure('Name', 'Convergence');
            p = semilogy(mvs, dat,'DisplayName',str);
            convergence_axes = gca;
            xlabel('Number of MAT.vec');
            ylabel('Relative residual');
            title('Convergence');
            legend();
            hold on;
        else
            p = semilogy(convergence_axes, mvs, dat,'DisplayName',str);
        end
        %pick a random color
        col = rand(1,3);
        while sum(col) < 1
            col = rand(1,3);
        end
        p.Color = col;
        p.LineWidth = 2;
        p.LineStyle = lin;
    end

    function report(s, mv, v, r, rv0, lin)
        if nargin < 6
            lin = '-';
        end
        tr = norm(system*v-rhs);
        fprintf(1, '    %-14s     converged after %4i matvecs residual %f %f error %e\n', s, mv, tr, r, norm(v-solution));
        convergence_plot(s, mv, rv0./rv0(1), lin);
        mulcnt = 0;
    end

if nargin < 1 && dim > 0
    sys = [];
end
end
