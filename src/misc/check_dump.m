function [sys, rhs, con] = check_dump(path)

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
    %drg = fread(fopen('drg.bin'), dim, 'double');
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
    fprintf(1, '    norm(MAT-transpose(MAT)) = %f\n', vecnorm(ela-ela'));

    fprintf(1, '    norm8(matrix - reconstituted_matrix) : %e\n', err0);
    if ( err0 > 1e-8 )
        imshow(abs(mat));
        set(gcf, 'name', 'Reconstituted matrix');
    end

    fprintf(1, '    norm8(preconditionner - reconstituted_preconditionner) : %e\n', max(max(abs(PRC-con))));
    fprintf(1, '    norm(rhs) = %f\n', norm(rhs));
    fprintf(1, '    norm(sol) = %f\n', norm(sol));
    fprintf(1, '    norm(con*rhs) = %f\n', norm(con*rhs));


end

figure('name', 'Mecable footprint'); imshow(PRC);
figure('name', 'System matrix structure'); spy(sys);


%% Getting the solution with high precision

reltol = abstol / norm(rhs);

maxit = dim;
system = sparse(sys);

solution = bicgstab(system, rhs, abstol*0.001, maxit);
     
fprintf(1, '    norm(sol - matlab_sol) = %f\n', norm(sol-solution));
    
if 0
    figure;
    plot(solution, sol, 'k.');
    xlabel('matlab solution');
    ylabel('cytosim solution');
    xl = xlim;
    ylim(xl);
    drawnow;
end

%% BCGS without preconditionning
	
[vec,~,~,itr,rv0] = bicgstab(system, rhs, reltol, maxit);
report('bicgstab', 2*itr, vec, rv0(end));

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
        fprintf(1, '%-14s     converged after %4i vecmuls residual %f %f error %e\n', s, i, tr, r, norm(v-solution));
    end

end
