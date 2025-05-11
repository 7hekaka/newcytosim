function netsolve

% simulate two flexible fibers, connected
% by a Hookean link that moves along both fibers.
% Test for implicit integration 
% 
% Francois Nedelec, 2000-2001
%
% This method was later described in:
% Collective Langevin dynamics of flexible cytoskeletal fibers
% Nedelec F and Foethke D
% New Journal of Physics 9 (11) 427, Nov 2007.
% Note that this MATLAB code also runs with the free clone Octave


M        = 16;         %number of points for first fiber;
P        = 16;         %number of points for second Fiber
L        = 0.5;        %length of segments
rigid    = 5;          %rigidity of the polymer

stiff    = 100;        %stiffness of link
v        = 1;          %speed of dispacement of the links
ab       = [ M-1, P-1] * L/2;   % initial positions of the links: middle


fb       = 0.5;        %feed back coefficient to re-establish length constraint

dt       = 1e-2;       %time-step
nbsteps  = 300;
rec      = 100;

%mobility of the segment, of length L
mu    = 10/L;
LN    = max(M,P) * L;
mudt  = mu * dt;

figure('Position',[ 100, 100, 600, 600], 'MenuBar','None', 'Name','Implicit');
%set(gca,'Position',[0 0 10 10]);


%vector x contains [x,y] of all successive points

x = zeros((M+P)*2, 1);          %positions of points

%initial position is an <X>:
d = L / 2 * [ sqrt(2), sqrt(2) ];
for i=1:M
   x(2*i-1:2*i) = ( i - (M+1)/2 ) * d;
end
d = L / 2 * [ sqrt(2), -sqrt(2) ];
for i=M+1:M+P
   x(2*i-1:2*i) = ( i - (M+(P+1)/2) ) * d;
end


%building the matrix for bending rigidity:
rigidLLL = rigid / ( L ^ 3 );
R1 = rigidLLL * RigidityMatrix(M);
R2 = rigidLLL * RigidityMatrix(P);

R = zeros( 2*(M+P), 2*(M+P) );
R(1:2:2*M, 1:2:2*M) = R1;
R(2:2:2*M, 2:2:2*M) = R1;
R(2*M+1:2:2*(M+P), 2*M+1:2:2*(M+P)) = R2;
R(2*M+2:2:2*(M+P), 2*M+2:2:2*(M+P)) = R2;


%integrating the motion:
for t=0:nbsteps
   
   %display the position
   if ( 1 )
      plot( x(1:2:2*M-1), x(2:2:2*M), '-ko', 'LineWidth', 2);
      hold on;
      plot( x(2*M+1:2:2*(M+P)-1), x(2*M+2:2:2*(M+P)), '-ks', 'LineWidth', 2);
      axis( [ -LN/2 LN/2 -LN/2 LN/2 ]);
      hold off;
      pause(0.01);
   end
   
   %matrix B for Hookean connections:
   B   = zeros(2*(M+P), 2*(M+P));
   
   % abscissa where the link is attached:
   ab  = min( ab + v * dt, [L*(M-1)-0.01, L*(P-1)-0.01] );
   abN = floor(ab/L);
   abf = ab - L*abN;
   abN = abN + [ 1, 1+M ];

   indx = [ abN(1), abN(1)+1, abN(2), abN(2)+1 ];
   F = [ -1+abf(1)/L, -abf(1)/L,  1-abf(2)/L,  abf(2)/L ];
   w = [  1-abf(1)/L;  abf(1)/L; -1+abf(2)/L; -abf(2)/L ];
   
   B(2*indx-1, 2*indx-1) = w * F;
   B(2*indx,   2*indx  ) = w * F;
   B = stiff * B;
   
   % columns and lines should summ up to zero:
   %sum(B,1)
   %sum(B,2)

   %building the jacobian of length constraints:
   J = zeros(M+P-2, 2*(M+P));
   c = 1;
   for r = [ 1:M-1, M+1:M+P-1]
      J(c, 2*r-1) = x(2*r-1) - x(2*r+1);
      J(c, 2*r+0) = x(2*r+0) - x(2*r+2);
      J(c, 2*r+1) = x(2*r+1) - x(2*r-1);
      J(c, 2*r+2) = x(2*r+2) - x(2*r+0);
      c = c + 1;
   end   
   
   %J is rectangular, and J*J' is symetric definite positive
   
   JJJJ = ( eye( 2*(M+P) ) - J' * inv( J * J' ) * J );
   
   % final matrix:
   mat = eye(2*(M+P)) - mudt .* JJJJ * ( R + B );
   
   %implicit, first order method:
   x  = mat \ x;
   
   if ( fb > 0 )
      %Feed-back on the length constraints:
      dx = zeros(2*(M+P), 1);
      for i = [ 1:M-1, M+1:M+P-1]
         d  = sqrt( ( x(2*i+1) - x(2*i-1) )^2 + ( x(2*i+2) - x(2*i) )^2 );
         dl = ( d - L ) / d;
         dx(2*i-1) = dx(2*i-1) + dl * ( x(2*i+1) - x(2*i-1) );
         dx(2*i)   = dx(2*i)   + dl * ( x(2*i+2) - x(2*i) );
         dx(2*i+1) = dl * ( x(2*i-1) - x(2*i+1) );
         dx(2*i+2) = dl * ( x(2*i) - x(2*i+2) );
      end      
      x = x + fb * dx;
   end
   
end


return;


function R = RigidityMatrix( N )

R = zeros(N,N);

if ( N < 3 ) return; end

if ( N == 3 )
   R = [ -1 2 -1; 2 -4 2; -1 2 -1 ]; 
   return;
end

if ( N == 4 )
   R = [ -1 2 -1 0; 2 -5 4 -1; -1 4 -5 2; 0 -1 2 -1 ]; 
   return;
end

for i=1:N;     R(i, i)   = -6;  end
for i=1:N-1;   R(i, i+1) =  4; R(i+1, i) =  4;  end
for i=1:N-2;   R(i+2, i) = -1; R(i, i+2) = -1;  end
R(1:2, 1:2)     = [-1 2; 2 -5];
R(N-1:N, N-1:N) = [-5 2; 2 -1];
return;

      
