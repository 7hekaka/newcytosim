function [sys, rhs, con] = check_dump(path)

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


%% Loading System

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

mulcnt = 0;
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
convergence_axes = [];

mulcnt = 0;
[vec,~,res,itr,rv0] = bicgstab(@multiply, rhs, reltol, maxit);
report('bicgstab', mulcnt, vec, res, rv0, '-.');

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

mulcnt = 0;
[vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @precondition);
report('P bicgstab', mulcnt, vec, res, rv0);

% checking the reconstituted block preconditionner:
mulcnt = 0;
[vec,~,res,~,rv0] = bicgstab(@multiply, rhs, reltol, maxit, @preconditionPRC);
report('P bicgstab', mulcnt, vec, res, rv0);

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
