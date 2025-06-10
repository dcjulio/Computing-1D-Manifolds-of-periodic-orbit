function manif = init_manif_periodic(opts)

%---- manif.name: Name of the manifold---% 
% for example Ws_pmin_pos is the stable manifolds of the fixed point pmin, the branch to the positive values
%
%---- manif.orientability: Orientability of the manifold---% 
%
%---- manif.fixp: Information of the fixed points associated to the manifold---% 
%
%---- manif.stability: Stability of the manifold---% 
%
%---- manif.points: Coordinates of the manifold ---% 
%
%---- manif.system_info: contains general info about the map ---% 
% system_info.par: The parameter values
% system_info.fixp: the fixed points with their  eigensystem, orientability, stability, etc
%
%---- manif.grow_info: information for the algorithm---% 
%
%%
%the map function where the system is defined (example StdHenon3D)
thesystem=opts.thesystem;   
%% Initializate field names and the structure 'manif'

%names of the fields
names = {
    'name'
    'par' % parameters
    'type' %periodic orbit or fixed point
    'orientability'
    'stability'   % stability
    'dimension'
    'per_orbit'   % periodic orbit
    'grow_info'    % options for computation
    }; 

manif=struct();
n=numel(names);

for k=1:n
    manif.(names{k})=[];
end

manif.type = sprintf('period-%i orbit', numel(opts.per_orbit.coord.x));
%% Default accuracty conditions

% default acc conditions
manif.grow_info.alphamax=0.3;
manif.grow_info.deltalphamax=10^(-3); 
manif.grow_info.deltamin=10^(-6);
manif.grow_info.deltamax=10^(-2);    
manif.grow_info.init_step=10^(-7);  

%% Update accuracy conditions if needed
%rewrite them if we other values are defined
manif.grow_info.thesystem = thesystem;

if isfield(opts,'accpar')
    if isfield(opts.accpar,'alphamax')
        manif.grow_info.alphamax=opts.accpar.alphamax;
    end
    if isfield(opts.accpar,'deltalphamax')
        manif.grow_info.deltalphamax=opts.accpar.deltalphamax;
    end
    if isfield(opts.accpar,'deltamin')
        manif.grow_info.deltamin=opts.accpar.deltamin;
    end
    if isfield(opts.accpar,'deltamax')
        manif.grow_info.deltamax=opts.accpar.deltamax;
    end
    if isfield(opts.accpar,'init_step')
        manif.grow_info.init_step=opts.accpar.init_step;
    end
end


%% Parameter values
manif.par=opts.par;  % parameters

%% General info of periodic orbit
manif.per_orbit.name=opts.per_orbit.name;
manif.per_orbit.coord_original=opts.per_orbit.coord;
manif.per_orbit.coord_compactified=thesystem.compactify(opts.per_orbit.coord); %periodic orbit in compactified coordinates

% eigensystem
[manif.per_orbit.eigval,manif.per_orbit.eigvec_x0]=eigensystem(opts.per_orbit.coord,opts);
eigval=manif.per_orbit.eigval;
eigvec_x0=manif.per_orbit.eigvec_x0;

% computes the dimension and orientation properties of the stable manifold
manif.per_orbit.Smanifold.dimension=sum(abs(eigval)<1);
manif.per_orbit.Smanifold.orientability=orientability(eigval,'Smanifold');
manif.per_orbit.Smanifold.eigval=eigval(abs(eigval)<1);
manif.per_orbit.Smanifold.eigvec_x0=eigvec_x0(abs(eigval)<1,:);

% computes the dimension and orientation properties of the unstable manifold
manif.per_orbit.Umanifold.dimension=sum(abs(eigval)>1);
manif.per_orbit.Umanifold.orientability=orientability(eigval,'Umanifold');
manif.per_orbit.Umanifold.eigval=eigval(abs(eigval)>1);
manif.per_orbit.Umanifold.eigvec=eigvec_x0(abs(eigval)>1,:);

% computes the dimension and orientation properties of the stable manifold
manif.per_orbit.Smanifold.dimension=sum(abs(eigval)<1);
manif.per_orbit.Smanifold.orientability=orientability(eigval,'Smanifold');
if strcmp(opts.stability,'Smanifold')
    manif.dimension=manif.per_orbit.Smanifold.dimension;
    manif.orientability=manif.per_orbit.Smanifold.orientability;
end

% computes the dimension and orientation properties of the unstable manifold
manif.per_orbit.Umanifold.dimension=sum(abs(eigval)>1);
manif.per_orbit.Umanifold.orientability=orientability(eigval,'Umanifold');
if strcmp(opts.stability,'Umanifold')
    manif.dimension=manif.per_orbit.Umanifold.dimension;
    manif.orientability=manif.per_orbit.Umanifold.orientability;
end


%% Stability and orientability and dimension of the manif to compute

manif.stability=opts.stability;

%if the field branch is defined, then follow that definition to know which branch to compute
if isfield(opts,'branch')
   if strcmp(opts.branch,'pos')
       manif.grow_info.init_step=abs(manif.grow_info.init_step);
   elseif strcmp(opts.branch,'neg')
       manif.grow_info.init_step=-abs(manif.grow_info.init_step);
   end
end

%% Name of the manifold. Example: Ws_pmin_pos

% defining the name of the manifold
manif.name = sprintf('W%s_%s', lower(manif.stability(1)),opts.per_orbit.name);  

%% algorithm information
manif.grow_info.stability=manif.stability; % stability of the manifold
manif.grow_info.orientability=manif.orientability; % orientability of the manifold
manif.grow_info.dimension=manif.dimension; % dimension of the manifold

manif.grow_info.eigval=manif.per_orbit.(manif.stability).eigval; 
manif.grow_info.eigvec=manif.per_orbit.(manif.stability).eigvec_x0; 

manif.grow_info.max_funditer=opts.max_funditer; % number of iteration of the algorithm
manif.grow_info.user_arclength=opts.user_arclength;

   
%----------------------------------------------
%-------------- FUNCTIONS ---------------------
%----------------------------------------------

% > --------  eigenvalue and eigenvector in compactified coordinates
    function [eigval,eigvec_x0]=eigensystem(per_orbit,opts)
    % Input: 
                % per_orbit: .x, .y and .z coordinatesof the periodic orbit
    % Output: 
                % eigval: the eigenvalue of the periodic orbit
                % eigenval: eigenvalues associated with the first point of the orbit

    syms x y z
    points=struct('x',x,'y',y,'z',z);
   
    system=opts.thesystem;
    period=numel(per_orbit.x);
    
    %Jacobian J_f of the original system (uncompactified)
    F  = system.ff(points,opts);
    JF = jacobian([F.x, F.y, F.z], [x, y, z]);

    %Jacobian J_t of the compactification
    T  = system.compactify(points);
    JT = jacobian([T.x, T.y, T.z], [x, y, z]);

    % The Jacobian of a periodic point x0 of period k is
    %J_fk(x0)= J_f(x(k-1)) * J_f(x(k-2)) * ... * J_f(x1) * J_f(x0) 

    %compute J_fk(x0)
    JFk = 1; 
    for i=1:period
        x   = per_orbit.x(i);
        y   = per_orbit.y(i);
        z   = per_orbit.z(i);
        JFk = double(subs(JF))*JFk;
    end

    %Compute the eigenval and eigenvec of the original system JFk(x0)
    [eigvecF,D]=eig(JFk);
    eigval=diag(D)';

    % Compute eigenvec in compactified coordinates of first point in the orbit (x0)
    % compactified coordinates at point x0
    x   = per_orbit.x(1);
    y   = per_orbit.y(1);
    z   = per_orbit.z(1);
    JTp = double(subs(JT)); %compactified coordinates at point x0

    %Transform the eigenvec to the compactified coordinates
    eigvec_x0=JTp*eigvecF;
    eigvec_x0=normc(eigvec_x0)'; %normalize each column
end

%----------------------------------------------
%----------------------------------------------
%----------------------------------------------


% > -------- orientability
function [orientability]=orientability(eigval,Stab)

    if strcmp(Stab,'Smanifold')
        if prod(eigval(abs(eigval) < 1))>0
            orientability='orientation-preserving';
        else
            orientability='orientation-reversing';
        end
    end


    if strcmp(Stab,'Umanifold')
        if prod(eigval(abs(eigval) > 1))>0
            orientability='orientation-preserving';
        else
            orientability='orientation-reversing';
        end
    end
    
end  

% > -------- oarclength
function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)

arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength
    
end