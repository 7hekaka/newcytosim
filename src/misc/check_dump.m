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

abstol = 0.001;


%% Loading

if isfolder(path)

    cwd = pwd;
    cd(path);
    
    ord = load('ord.txt');
    dim = ord(1);
    stp = load('stp.txt');
    time_step = stp(1);
    if ( length(stp) > 1 )
        abstol = stp(2);
    end
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


[vec,~,~,itr,rv0] = bicgstabl(system, rhs, reltol, maxit);
report('bicgstabl', 2*itr, vec, rv0(end));
convergence_plot(rv0/rv0(1), 'k:');


%% WITH TRADITIONAL PRECONDITIONNING

[vec,~,~,itr,rv0] = bicgstab(system, rhs, reltol, maxit, @precondition);
report('P BCGS', 2*itr, vec, rv0(end));
convergence_plot(rv0/rv0(1),'k');

%% subfunctions

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
        % crop data:
        dat = data(1:min(256,length(data)));
        if isempty(convergence_axes)
            semilogy(dat, col, 'Linewidth', 2);
            convergence_axes = gca;
        else
            semilogy(convergence_axes, dat, col, 'Linewidth', 2);
        end
    end

    function report(s, i, v, r)
        tr = norm(system*v-rhs);
        fprintf(1, '    %-14s     converged after %4i vecmuls residual %f %f error %e\n', s, i, tr, r, norm(v-solution));
    end

end
