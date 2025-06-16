%% Things to add:
% give an approximate and find with newton the true %% easiest
% ask for all periodic points of certain period %% a tiny bit harder

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    %--- Define coordinates of the periodic orbit in original coordinates
    % % % orientation-reversing period-three orbit
    p3_x = [2.620266326709355, -2.2163325979982256, -1.4982100829623703];
    p3_y = [-2.2163325979982256, -1.4982100829623703, 2.620266326709355];
    p3_z = [3.2542647169820524, 1.6887850623802372, 0.5283319918939143];

    % % orientation-preserving period-three orbit
    % p3_x = [2.269075165256071, -1.408946374686311, 1.5341475636814805];
    % p3_y = [-1.408946374686311, 1.5341475636814805, 2.269075165256071];
    % p3_z = [-2.8587492611614835, -4.8394454880800915, -4.273187022014629];

    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength=100; % What is the approximate arclength of the manifold
    opts.per_orbit.name ='p3';
    opts.per_orbit.coord = struct('x',p3_x,'y',p3_y,'z',p3_z);
    opts.stability='Smanifold';

    %--- Number of iterations used to compute the manifold
    opts.max_funditer=30; % how many times (max) the algorithm iterates the fundamental domain

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-7;

%% Computing the manifold

    opts.branch = 'pos'; %which branch: 'pos', 'neg' or '' to consider sign of initial step.
    manif = GrowFundCurv1D_periodic(opts);
    manif = add_branch_periodic(manif, opts, 'neg');

    %% Computing intersection points
    angle=pi/2; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    manif=inter_plane_periodic(manif,angle);

%% Plot
    manifplot_periodic(manif);

%% Epsilon pseudo orbit (orientation-preserving)
    
    branch    = 'pos';
    idxPO     = 1;
    idxpoint = manif.inter.points{idxPO}.(branch).idx(end);
    orbit     = eps_pseudo_orbit_periodic(manif, idxPO, idxpoint, branch);
    
    % plot the epsilon orbit. Starting and end point are colored in solid red.
    hold on
    plot3(orbit.x([1,end]),orbit.y([1,end]),orbit.z([1,end]),'r.','MarkerSize',27) %epsilon orbit
    plot3(orbit.x,orbit.y,orbit.z,'ko--','LineWidth',1.2) %epsilon orbit