function [SYS, rhs, CON] = cytosim_dump(path)

% This is used to explore Cytosim's linear system in matlab
% - load the matrices and vector from Cytosim's dump
% - plot convergence pattern of BICGstab, with and without preconditionning
%
% F. Nedelec, 16.10.2014, 03.2018, 06.2018, 26.01.2019, 30.06.2019, 11.08.2019, 
% 17.08.2019, 7.01.2020, 03.06.2020, 19.06.2020, 27.12.2020, 13.01.2021, 
% 21.08.2021

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
    obj = fread(fopen('obj.bin'), dim, 'uint32');  % object id
    mob = fread(fopen('mob.bin'), dim, precision);
    SYS = fread(fopen('sys.bin'), [dim, dim], precision);
    ela = fread(fopen('ela.bin'), [dim, dim], precision);  % elasticity matrix
    PRJ = fread(fopen('prj.bin'), [dim, dim], precision);  % projection matrix
    CON = fread(fopen('con.bin'), [dim, dim], precision);  % preconditionner
    %pts = fread(fopen('pts.bin'), dim, precision);
    rhs = fread(fopen('rhs.bin'), dim, precision);
    sol = fread(fopen('sol.bin'), dim, precision);
    
    cd(cwd);
else
    error(['cannot find dump directory ',path]);
end

fprintf(1, '----------------------- loaded system of size %i with time_step %f -----------------------\n', dim, time_step);

%% Report some

if 0
    fprintf(2, '    norm(rhs) = %f\n', norm(rhs));
    fprintf(2, '    norm(sol) = %f\n', norm(sol));
    fprintf(2, '    norm(ela) = %f\n', norm(ela, 2));
    fprintf(2, '    norm(PRJ) = %f\n', norm(PRJ, 2));
    mag = time_step * max(max(abs(sparse(PRJ) * sparse(ela))));
    fprintf(1, '    norm8(system-eye) : %e\n', mag);
end

if 0
    figure('name', 'System matrix'); imshow(abs(SYS));
    figure('name', 'Projection matrix'); imshow(abs(PRJ));
    figure('name', 'Elasticity matrix'); imshow(abs(ela));
    return;
end

%% Check matrix

if ( 0 )
    
    MAT = eye(dim) - time_step * PRJ * ela;
    err0 = norm(MAT-SYS, 1);
    
    nbo = 0;
    nbv = 0;
    for o = 0:max(obj)
        n = sum(obj==o);
        if n > 0
            nbo = nbo + 1;
            nbv = nbv + n^2;
        end
    end
    
    fprintf(1, '%i mecables with %i block scalars (%s)\n', nbo, nbv, precision);
    fprintf(2, '    norm(ela-transpose(ela)) = %f\n', norm(ela-ela',1));

    fprintf(2, '    norm8(matrix - reconstituted_matrix) : %e\n', err0);
    if ( err0 > 1e-5 )
        figure;
        imshow(abs(MAT));
        set(gcf, 'name', 'Reconstituted matrix');
    end
end

if ( 0 )
    figure('Position', [50 50 1000 1000]);
    plot(reshape(MAT,1,dim*dim), reshape(SYS,1,dim*dim), '.')
    xl = xlim;
    ylim(xl);
    xlabel('Reconstituted matrix');
    ylabel('cytosim matrix');
end
if ( 0 )
    figure('name', 'System matrix values');
    plot(abs(SYS), '^b'); hold on;
    plot(abs(MAT), 'vr');
end
if ( 0 )
    % checking the formula to calculate transposed(SYS)*X:
    MAT = speye(dim) - time_step * sparse(PRJ) * sparse(ela);
    TAM = speye(dim) - time_step * sparse(ela) * sparse(PRJ);
    fprintf(2, ' I - tau * PRJ * ela  has %9i elements\n', nnz(MAT));
    fprintf(2, ' I - tau * ela * PRJ  has %9i elements\n', nnz(TAM));
    figure('name', 'System structure', 'Position', [50 100 2000 1000]);
    subplot(1,2,1);
    imshow(abs(MAT));
    subplot(1,2,2);
    imshow(abs(TAM));
end

if ( 0 )
    uuu = speye(dim) - time_step * sparse(ela);
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
end

drawnow;

%% RECALCULATE SOLUTION

reltol = abstol / norm(rhs);

system = sparse(SYS);

mulcnt = 0;
solution = bicgstab(@multiply, rhs, reltol*0.001, dim);
maxit = mulcnt;

fprintf(2, '    %i matvecs; norm(sol - matlab_sol) = %f\n', mulcnt, norm(sol-solution));
    
if 0
    figure('Name', 'Validation of solution');
    plot(solution, sol, 'k.');
    xlabel('matlab solution');
    ylabel('cytosim solution');
    xl = xlim;
    ylim(xl);
    drawnow;
end

if 1
    % Check equivalent system obtained by Woodbury's identity
    TAM = speye(dim) - time_step * sparse(ela) * sparse(PRJ);
    [vec,~,~,itr] = bicgstab(TAM, ela*rhs, reltol*norm(rhs)/norm(ela*rhs), maxit);
    wood = rhs + time_step * PRJ * vec;
    fprintf(2, '    %i matvecs; norm(woodbury_sol - sol) = %f\n', 2*itr, norm(wood-solution));
end

%% Try different iterative solvers

convergence_axes = [];

mulcnt = 0;
[vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit);
report('bicgstab', mulcnt, vec, res, rv0, '-');

if 0
    % Matlab BiCG
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicg(@multiply, rhs, reltol, maxit);
    report('bicg', mulcnt, vec, res, rv0, '-');
end
if 1
    % Matlab BiCGStab(L)
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstabl(@multiply, rhs, reltol, maxit);
    report('bicgstab(1)', mulcnt, vec, res, rv0, '-');
end
if 0
    % Matlab CGS
    mulcnt = 0;
    [vec,~,res,itr,rv0] = cgs(@multiply, rhs, reltol, maxit);
    report('CGS', mulcnt, vec, res, rv0, '-');
end
if 0
    % BiCORSTAB
    mulcnt = 0;
    [vec,~,res,itr,rv0] = BiCORSTAB(@multiply, rhs, reltol, maxit);
    report('BiCORstab', mulcnt, vec, res, rv0, ':');
end
if 0 %% no preconditionning
    for i = 1:4
        OPT.Tol = reltol;
        OPT.ell = min(2^i, dim);
        OPT.MaxIt = maxit;
        [vec,rv0,itr] = cgstab(system, rhs, [], [], OPT);
        report(sprintf('cgstab(%i)', OPT.ell), 2*itr, vec, rv0(end,1), rv0(:,1), ':');
    end
end
if 1
    % GMRES
    for i = 4:6
        mulcnt = 0;
        RS = min(2^i, dim);
        [vec,~,res,itr,rv0] = gmres(@multiply, rhs, RS, reltol, maxit);
        report(sprintf('GMRES %03i', RS), mulcnt, vec, res, rv0, '-.');
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
if 0
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

title('Convergence');
drawnow;


%% Check convergence on each subspace
fprintf(1, '  --  --  --  --  --  --  --  -- REDUCED SYSTEMS --  --  --  --  --  --  --  --  --\n');

if 0
    % solve a system using the block corresponding to each Mecable
    off = 0;
    for o = 0:max(obj)
        n = sum(obj==o);
        if n > 1
            fprintf(2, '  Block %i size %3i:', o, n);
            M = SYS(off+(1:n), off+(1:n));
            [vec,~,res,itr] = bicgstab(M, rhs(off+(1:n)), reltol, maxit);
            err = norm(vec-solution(off+(1:n)));
            fprintf(1, ' bCGstab reached residual %f within %2.0f matvecs (error %f)\n', res, 2*itr, err);
        end
        off = off + n;
    end
    return;
end


%% USING LOADED PRECONDITIONNER CALCULATED BY CYTOSIM
fprintf(1, '  --  --  --  --  --  --  --  -- PRECONDITIONNED --  --  --  --  --  --  --  --  --\n');

% Calculate block-diagonal preconditionner
BDP = zeros(dim);
for o = 0:max(obj)
    x = find(obj==o);
    if ~isempty(x)
        BDP(x,x) = inv(SYS(x,x));
    end
end

%figure('name', 'Mecable footprint'); imshow(BDP);
fprintf(2, '    norm8(preconditionner - reconstituted_preconditionner) : %e\n', norm(BDP-CON,1));

fprintf(2, '    Elasticity            has %9i elements\n', nnz(ela));
fprintf(2, '    Mobility/Projection   has %9i elements\n', nnz(PRJ));
fprintf(2, '    Given Preconditionner has %9i elements\n', nnz(CON));
fprintf(2, '    Block preconditionner has %9i elements\n', nnz(BDP));

convergence_axes = [];

mulcnt = 0;
[vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @precondition);
report('P bicgstab', mulcnt, vec, res, rv0);

if 0
    % checking the reconstituted block preconditionner:
    mulcnt = 0;
    [vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @preconditionBDP);
    report('R bicgstab', mulcnt, vec, res, rv0);
end
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
if 0
    % CORS
    mulcnt = 0;
    [vec,~,res,itr,rv0] = CORS(@multiply, rhs, reltol, maxit, @precondition);
    report('P CORS', mulcnt, vec, res, rv0);
    % BiCOR
    mulcnt = 0;
    [vec,~,res,itr,rv0] = BiCOR(@multiply, rhs, reltol, maxit, @precondition);
    report('P BiCOR', mulcnt, vec, res, rv0);
    % BiCGCR2
    [vec,~,res,itr,rv0] = BiCGCR2(system, rhs, reltol, maxit, @precondition);
    report('P BiCGCR2', 2*itr, vec, res, rv0);
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

fprintf(2, 'Mecable point mobility coefficients : %f +/- %f\n', mean(mob), var(mob));

if 0
    
    avm = time_step * mean(mob);
    DRY = eye(dim) - diag(avm.*ones(length(ela), 1)) * ela;
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
    
    avm = time_step * mean(mob);
    SML = eye(dim) - diag(avm.*ones(length(ela), 1)) * ela;
    
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
    
    L = diagonal_expand(SML, ord(2));
    
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, L');
    report('s bicgstab', mulcnt, vec, res, rv0);
    
    mulcnt = 0;
    [vec,~,res,itr,rv0] = bicgstabl(@multiply, rhs, reltol, maxit, L, L');
    report('s bicgstab(1)', mulcnt, vec, res, rv0);
end


%% incomplete LU factorization

% the WET matrix is not necessarily symmetric because the diagonal matrix
% may not commute with 'ela'. However, if all point drags are equal, then
% surely WET is symmetric and we can use incomplete Cholesky factorization

if 0
    
    % without the projection:
    WET = speye(dim) - time_step .* spdiags(mob) * sparse(ela);
        
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
    clear OPT;
    OPT.type = 'nofill';
    OPT.milu = 'off';
    [L, U] = ilu(system, OPT);
    fprintf(2, 'incomplete LU    has %i + %i elements\n', nnz(L), nnz(U));
    
    mulcnt = 0;
    [vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, L, U);
    report('iLU bicgstab', mulcnt, vec, res, rv0);
    
    LU = L*U;
    fprintf(2, 'incomplete L*U   has %i elements\n', nnz(LU));

end
if 0
    
    % try to threshold the preconditionner:
    aval = abs(nonzeros(LU));
    fprintf(2, ' %i LU values in [ %f %f ]\n', numel(aval), min(aval), max(aval));
    
    for i = -2 : 3
        l = 10^i;
        u = 10^(i+1);
        c = sum((l<aval).*(aval<u));
        fprintf(2, '    %9i LU values in [ %.0f %.0f ]\n', c, l, u);
    end
    
    %figure('Name', 'iLU elements'); histogram(nL, 1000); 
    
    rLU = LU .* ( abs(LU) > 0.1 );
    fprintf(2, 'thresholded iLU (>0.1)     has %i elements\n', nnz(rLU));
    mulcnt = 0;
    [vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, rLU);
    report('tiLU bicgstab', mulcnt, vec, res, rv0);
    
    rLU = LU .* ( abs(LU) > 0.2 );
    fprintf(2, 'thresholded iLU (>0.2)     has %i elements\n', nnz(rLU));
    mulcnt = 0;
    [vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, rLU);
    report('tiLU bicgstab', mulcnt, vec, res, rv0);

end
if 0
    
    % try to reduce the preconditionner:    
    rLU = dimension_collapse(LU, ord(2));
    fprintf(2, 'collapsed  rLU   has %i elements\n', nnz(rLU));

    if 0
        xLU = diagonal_expand(rLU, ord(2));
        mulcnt = 0;
        [vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, xLU);
        report('pilu bicgstab', mulcnt, vec, res, rv0);
    end
    xLU = diagonal_expand(rLU, ord(2));
    % replace block diagonal by full preconditionner:
    xLU(logical(CON)) = 0;
    xLU = xLU + sparse(CON);
    mulcnt = 0;
    [vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, xLU);
    report('xilu bicgstab', mulcnt, vec, res, rv0);

end

title('Convergence (preconditionned)');
drawnow;

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
            y = CON * x;
        else
            y = CON' * x;
        end
    end

    function y = preconditionBDP(x, mode)
        if nargin < 2 || strcmp(mode, 'notransp')
            y = BDP * x;
        else
            y = BDP' * x;
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
            legend();
            hold on;
        else
            p = semilogy(convergence_axes,mvs,dat,'DisplayName',str);
        end
        %pick a random color
        col = rand(1,3);
        while sum(col) < 1
            col = rand(1,3);
        end
        p.Color = col;
        p.LineWidth = 4;
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


    function res = dimension_collapse(mat, ord)
        res = sparse(size(mat,1)/ord, size(mat,2)/ord);
        if ( ord == 3 )
            for d = 1:3
                iii = d:3:size(mat,1);
                res = res + mat(iii, 1:3:end) + mat(iii, 2:3:end) + mat(iii, 3:3:end);
            end
            res = (1/9) * res;
        elseif ( ord == 2 )
            for d = 1:2
                iii = d:2:size(mat,1);
                res = res + mat(iii, 1:2:end) + mat(iii, 2:2:end);
            end
            res = (1/4) * res;
        end
    end

    function res = diagonal_expand(mat, D)
        res = sparse(size(mat,1)*D, size(mat,2)*D);
        if ( D == 3 )
            res(1:3:end, 1:3:end) = mat;
            res(2:3:end, 2:3:end) = mat;
            res(3:3:end, 3:3:end) = mat;
        elseif ( D == 2 )
            res(1:2:end, 1:2:end) = mat;
            res(2:2:end, 2:2:end) = mat;
        end
    end

    function res = block_expand(mat, D)
        zis = size(mat,1)*D;
        res = sparse(zis, size(mat,2)*D);
        if ( D == 3 )
            for d = 1:3
                iii = d:3:zis;
                res(iii, 1:3:end) = mat;
                res(iii, 2:3:end) = mat;
                res(iii, 3:3:end) = mat;
            end
        elseif ( D == 2 )
            for d = 1:2
                iii = d:2:zis;
                res(iii, 1:2:end) = mat;
                res(iii, 2:2:end) = mat;
            end
        end
    end

if nargin < 1 && dim > 0
    SYS = [];
end
end
