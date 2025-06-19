%% Things to add:
% give an approximate and find with newton the true %% easiest
% ask for all periodic points of certain period %% a tiny bit harder

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    % %--- Define coordinates of the periodic orbit in original coordinates
    % orientation-preserving
    p7_x = [1.77019027758592	1.60375145067744	1.09692420117424	2.51563186167503	-2.45748092382682	-2.59390204947524	-1.79108356512382];
    p7_y = [1.60375145067744	1.09692420117424	2.51563186167503	-2.45748092382682	-2.59390204947524	-1.79108356512382	1.77019027758592];
    p7_z = [-1.69319655834972	-0.428084419342221	0.583222897963574	3.21549933923131	1.40111828325075	-0.912560109574340	-2.88615569661303];

    % orientation-preserving
    % p7_x = [1.57324256130362	2.27122793209917	-1.43044908793855	1.47244702718601	2.46103447851263	-2.29842481258375	-1.82106696265445];
    % p7_y = [2.27122793209917	-1.43044908793855	1.47244702718601	2.46103447851263	-2.29842481258375	-1.82106696265445	1.57324256130362];
    % p7_z = [-2.55754432734339	-0.797825260712904	-2.38783940079404	-1.39296025376683	0.789482173992440	-1.35104620379283	-3.44232240720584];

    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength = 10; % What is the approximate arclength of the manifold
    opts.per_orbit.name ='p7';
    opts.per_orbit.coord = struct('x',p7_x,'y',p7_y,'z',p7_z);
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