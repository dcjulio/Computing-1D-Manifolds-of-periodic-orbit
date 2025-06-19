%% Things to add:
% give an approximate and find with newton the true %% easiest
% ask for all periodic points of certain period %% a tiny bit harder

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    % %--- Define coordinates of the periodic orbit in original
    
    % orientation-reversing
    p8_x = [1.76888928611210	1.46784036101440	1.51477788874350	1.46509583946946	1.59906081454622	1.20347575954193	2.27192785183111	-1.32269889178852];
    p8_y = [1.46784036101440	1.51477788874350	1.46509583946946	1.59906081454622	1.20347575954193	2.27192785183111	-1.32269889178852	1.76888928611210];
    p8_z = [-6.75095401588474	-6.63330445804729	-6.44518746091324	-6.26912911362643	-5.92389412180549	-5.90519718662466	-4.81430877211848	-7.09986941833069];

    % %orientation-preserving
    % p8_x = [1.80084852721732	1.48378904419337	1.45811551416654	1.62876243408885	1.10969827905101	2.47994099924457	-2.28301684344943	-1.75614820724719];
    % p8_y = [1.48378904419337	1.45811551416654	1.62876243408885	1.10969827905101	2.47994099924457	-2.28301684344943	-1.75614820724719	1.80084852721732];
    % p8_z = [-5.06567824050984	-4.59502484441845	-4.05591429913560	-3.23833472487387	-2.77630339079763	-0.851623069712588	-3.30496452710454	-5.72210563977264];


    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength = 10; % What is the approximate arclength of the manifold
    opts.per_orbit.name ='p8';
    opts.per_orbit.coord = struct('x',p8_x,'y',p8_y,'z',p8_z);
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