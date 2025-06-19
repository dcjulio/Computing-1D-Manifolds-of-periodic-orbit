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
    p6_x = [1.474386617060952 2.448982035466422 -2.239828995155547 -1.551528538179433 2.464707893761456 -1.409326440116204];
    p6_y = [2.448982035466422 -2.239828995155547 -1.551528538179433 2.464707893761456 -1.409326440116204 1.474386617060952];
    p6_z = [-0.644888819046554 1.417159924991936 0.027626884831551 -1.507325522448952 0.052987057843133 -1.324547147567191];

    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.6); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength = 10; % What is the approximate arclength of the manifold
    opts.per_orbit.name ='p6';
    opts.per_orbit.coord = struct('x',p6_x,'y',p6_y,'z',p6_z);
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