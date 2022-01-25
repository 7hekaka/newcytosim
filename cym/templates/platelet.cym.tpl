% Formation of a ring of microtubules
% Initial simulation by F. Nedelec May 2011
% First platelet simulation by Antonio Politi, Karin Sadoul, Romain Gibeaux, Aastha Mathur in Summer 2012
% Revised, F. Nedelec 4 Oct 2016; 17 Jan 2017; 20 Feb 2017
% Scaling implemented 31 Aug 2017, +25% 3 Sep 2017
% New Kinesin-Binder model, Jan. 2021

% The averaged sized platelet:
% length of ring         Lmb = 6.82 um      % mean of 4 EM +25% = 8.53
% radius of ring         R = 1.08 um        % mean of 4 EM +25% = 1.35
% Polymer                P = 102 um         % mean of 4 EM +25% = 127
% Number of microtubules N = 24             % mean of 4 EM 
% Microtubule length     L = 5.2 um         % mean of 4 EM +25% = 6.5
% Number of nucleators   G ~ 25             % super-resolution
% radius of cell         Rc = R*1.1 + 0.3   % light
% half-thickness         T = Rc * 0.125 um  % wisdom
% Monomer+Polymer        M = 2*P            % wisdom


% Scaling:

[[Lmb = round(random.gauss(10.2, 1.8), 3)]]% preconfig.perimeter = [[Lmb]]
[[R = round(Lmb/(2*math.pi), 3)]]% preconfig.radius = [[R]]
[[Rcell = round(R*1.1+0.1, 3)]]% preconfig.radius_cell = [[Rcell]]
[[S = Lmb/(8.53)]]% preconfig.scale = [[S]]
[[P = round(S*S*127,4)]]% preconfig.polymer = [[P]]
[[M = P*2.0]]% preconfig.tubulin = [[M]]
[[G = int(S*S*(random.gauss(30, 5.5)))]]% preconfig.nucleator = [[G]]
[[hydrol = 0.17]]% preconfig.hydrolysis = [[hydrol]]
[[K=int(4000*S*S)]]% preconfig.kinesin = [[K]]

set simul system
{
    steric = 1, 600
    time_step = 0.001
    viscosity = 0.5
    precondition = 1
    binding_grid_step = 0.256
    display = ( style=3; back_color=white; point_value=0.01; point_size=1; line_width=0.5;)
}

set space cell
{
    % volume ellipsoid = 4/3 PI * A * B * C
    shape = ellipse
    dimensions = [[Rcell]] [[Rcell]] [[Rcell*0.25]]
    display = ( color=0x00000044, white; visible=2; )
}

set fiber microtubule
{
    rigidity = 20
    steric = 1, 0.024
    segmentation = 0.160
    confine = inside, 1000, cell
    
    min_length = 0.032
    
    mesh = 1, 0.032
    mesh_aging_rate = 0.01;

    activity = dynamic
    unit_length = 0.008
    growing_speed = 0.128
    shrinking_speed = -0.56
    hydrolysis_rate = [[round(hydrol,4)]]
    growing_force = inf
    total_polymer = [[round(M)]]
    
    display = ( color=light_blue; line=2.5; )
}

%----------------------------------

set hand kinesin
{
    binding = 10, 0.07
    unbinding = 0.4, 3
    activity = move
    unloaded_speed =  0.4
    variable_speed = -0.2
    stall_force = 5
    display = ( color = green; )
}

set hand binder
{
    binding = 10, 0.07
    unbinding = 0.1, 5
    activity = bind
    display = ( color = gray; )
}

set couple complex
{
    hand1 = kinesin
    hand2 = binder
    length = 0.05
    stiffness = 200
    diffusion = 5
}

set hand docker
{
    binding = 10, 0.080
    unbinding = 0.1, 5
    activity = slide
    mobility = 0.1
    display = ( color=red; size=3; width=2; )
}

set hand nucleator
{
    binding = 0, 0.06
    unbinding = 0, inf
    hold_growing_end = 1;
    hold_shrinking_end = 0;

    activity = nucleate
    nucleate = 0.05, microtubule, ( length=0.064; plus_end=grow; )
    display = ( color=yellow; size=3; width=2; )
}

set couple gamma
{
    hand1 = docker
    hand2 = nucleator
    diffusion = 0.2
    stiffness = 500
    length = 0.06
}

%----------------------------------

new space cell
new [[G]] couple gamma
new [[K]] couple complex

run 250000 system
{
    nb_frames = 0
}

run 500000 system 
{
    nb_frames = 100
}



