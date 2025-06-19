%% Things to add:
% give an approximate and find with newton the true %% easiest
% ask for all periodic points of certain period %% a tiny bit harder

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    % %--- Define coordinates of the periodic orbit in original coordinates
    % orientation-reversing
    p10_x = [1.80883613823861	1.46213628176954	1.51950665206156	1.45245864980980	1.63451187497422	1.09263333562531	2.51579883138804	-2.45703376070103	-2.59175455064107	-1.78008152255837];
    p10_y = [1.46213628176954	1.51950665206156	1.45245864980980	1.63451187497422	1.09263333562531	2.51579883138804	-2.45703376070103	-2.59175455064107	-1.78008152255837	1.80883613823861];
    p10_z = [-4.58293490032142	-4.03738559861616	-3.32535606627784	-2.53796862972360	-1.41105048069410	-0.600627241207610	1.79504614193891	-0.302978390374339	-2.95532861909027	-5.32647586546669];

    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength = 10; % What is the approximate arclength of the manifold
    opts.per_orbit.name ='p10';
    opts.per_orbit.coord = struct('x',p10_x,'y',p10_y,'z',p10_z);
    opts.stability='Smanifold';

    %--- Number of iterations used to compute the manifold
    opts.max_funditer=100; % how many times (max) the algorithm iterates the fundamental domain

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-7; % pos value is positive branch and viceversa

%% Computing the manifold: Ws(pmin) orientation-preserving
     opts.branch = 'pos'; %which branch: 'pos', 'neg' or '' to consider sign of initial step.
     manif=GrowFundCurv1D_periodic(opts);

%% Computing intersection points
    angle=pi/2; %the angle of the plane from [-pi, pi]. (angle=pi/2: x==0 (y>0), angle=0: y==0 (x>0))
    manif=inter_plane_periodic(manif,angle);
%% Plot
    manifplot_periodic(manif);

%% Epsilon pseudo orbit (orientation-preserving)
    
branch    = 'pos';
idxPO     = 1;

if numel(manif.inter.points{idxPO}.(branch).idx)>0

    idxpoint  = manif.inter.points{idxPO}.(branch).idx(end);
    orbit     = eps_pseudo_orbit_periodic(manif, idxPO, idxpoint, branch);
    
    % plot the epsilon orbit. Starting and end point are colored in solid red.
    hold on
    plot3(orbit.x([1,end]),orbit.y([1,end]),orbit.z([1,end]),'r.','MarkerSize',27) %epsilon orbit
    plot3(orbit.x,orbit.y,orbit.z,'ko--','LineWidth',1.2) %epsilon orbit
end