%% Things to add:
% give an approximate and find with newton the true %% easiest
% ask for all periodic points of certain period %% a tiny bit harder

%% Adding the path
    clear all
    close all
    addpath('./GrowFundCurv1D_functions');

%% Options

    % %--- Define coordinates of the periodic orbit in original coordinates
    % orientation-preserving zero eigenvalue
    p20_x = [1.76627593644908	1.47149249950178	1.50482704297529	1.49404782087972	1.51637299603198	1.45239859064108	1.63562643509420	1.08900658762872	2.52337672157299	-2.49413205526508	-2.77770772557269	-2.76742059212671	-2.62530441605514	-1.86199709932061	1.52055812693817	2.44650211239848	-2.24154002405168	-1.55845231314513	2.44368839486809	-1.30407727726945];
    p20_y = [1.47149249950178	1.50482704297529	1.49404782087972	1.51637299603198	1.45239859064108	1.63562643509420	1.08900658762872	2.52337672157299	-2.49413205526508	-2.77770772557269	-2.76742059212671	-2.62530441605514	-1.86199709932061	1.52055812693817	2.44650211239848	-2.24154002405168	-1.55845231314513	2.44368839486809	-1.30407727726945	1.76627593644908];
    p20_z = [-4.41437368736107	-3.82575592533150	-3.08608006742252	-2.20924826002730	-1.13472491600078	0.0907286914401455	1.74450086482237	3.18240762541557	6.34226587207167	5.11658699122093	3.36219666389242	1.26721540454420	-1.10464593060211	-3.18757221604314	-2.30452853231360	-0.318932126377843	-2.62425857570509	-4.70756260399124	-3.20538672992139	-5.15054135317512];


    %--- Information of the system
    opts.thesystem=StdHenon3D_periodic; % What is the name of the system file
    opts.par=struct('a', 4.2,'b', 0.3, 'xi', 1.2); % The parameter values and names (has to match with the names defined in StdHenon3D)
    opts.user_arclength = 10; % What is the approximate arclength of the manifold
    opts.per_orbit.name ='p20';
    opts.per_orbit.coord = struct('x',p20_x,'y',p20_y,'z',p20_z);
    opts.stability='Smanifold';

    %--- Number of iterations used to compute the manifold
    opts.max_funditer=100; % how many times (max) the algorithm iterates the fundamental domain

    %--- Accuracy parameters (default)
    %opts.accpar.alphamax=0.3;
    %opts.accpar.deltalphamax=0.001; 
    %opts.accpar.deltamin=0.000001;
    %opts.accpar.deltamax=0.01;  

    %--- Initial step (default)
    %opts.accpar.init_step=10^-8; % pos value is positive branch and viceversa

%% Computing the manifold: Ws(pmin) orientation-preserving

    opts.branch = 'pos'; %which branch: 'pos', 'neg' or '' to consider sign of initial step.     
    manif=GrowFundCurv1D_periodic(opts);

    %% Add another branch
    manif = add_branch_periodic(manif, opts, 'neg');

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