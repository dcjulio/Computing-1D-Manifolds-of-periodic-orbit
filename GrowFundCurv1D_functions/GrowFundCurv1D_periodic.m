function Manif=GrowFundCurv1D_periodic(opts,eigenval_idx)
% Input:
% - opts: all the options needed to compute the manifold
% - eigenval_idx: optional argument. It is which eigenvector to consider for 2D manifolds
if nargin < 2
    eigenval_idx = [];
else
    if eigenval_idx~=1 && eigenval_idx~=2
        fprintf('\n ERROR! The optional second argument eigenval_idx can be either 1 or 2\n')
        return;
    end
end
    
    %% Initializing the manifold structure
    manif = init_manif_periodic(opts);
    % figure
    % hold on

    tic
    %%
    %saves the information of the system stored in the structure 'manif' and 'opts'
    par=manif.par;
    thesystem=opts.thesystem;
    per_orbit=manif.per_orbit.coord_compactified;
    period = numel(per_orbit.x);
    name = opts.per_orbit.name;

%% Error in case the given period-k orbit is not a periodic orbit
    
PO_uncomp=manif.per_orbit.coord_original;
PO.x=PO_uncomp.x(1); PO.y=PO_uncomp.y(1); PO.z=PO_uncomp.z(1);
orbit=PO;

for k=2:period
    PO = thesystem.ff(PO,opts);
    orbit.x(k)=PO.x;
    orbit.y(k)=PO.y;
    orbit.z(k)=PO.z;
end    
PO = thesystem.ff(PO,opts);
orbit.x(1)=PO.x;
orbit.y(1)=PO.y;
orbit.z(1)=PO.z;

distance=sum(sqrt((PO_uncomp.x-orbit.x).^2+(PO_uncomp.y-orbit.y).^2+(PO_uncomp.z-orbit.z).^2));
if distance > 10^-6
    fprintf('\n----  Error! The periodic orbit appears to be incorrect.\n The distances between the given orbit and the iterations from the first point are:\n');
    for k=1:period
        disp(sqrt((PO_uncomp.x(k)-orbit.x(k)).^2+(PO_uncomp.y(k)-orbit.y(k)).^2+(PO_uncomp.z(k)-orbit.z(k)).^2))
    end
    return
end

    %% Warnings in case the specifications are not correct for the computation of this manifold

    if strcmp(manif.orientability,'orientation-reversing')
        fprintf('\n----  Manifold is non-orientable. We will proceed computing both branches\n');
        mapiter=2;
    else
        mapiter=1;
    end

    % if computing unstable manifold, use image
    % if computing stable manifold, use preimage
    if strcmp(manif.stability,'Smanifold')
        sign = -1;
    elseif strcmp(manif.stability,'Umanifold')
        sign = 1;
    end

    if manif.grow_info.init_step>0
        branch1='pos';
        branch2='neg';
    else
        branch1='neg';
        branch2='pos';
    end
    

    %% Warnings if the manifold we want to compute is not 1D
    
    if manif.grow_info.dimension>1
        fprintf('\n----  The dimension of the manifold is %i and this algorithm is designed to compute one-dimensional manifolds. -----\n',manif.grow_info.dimension)
        fprintf(' To choose automatically which of the two eigenvalues to use, include a second argument to GrowFundCurv1D(opts,#) where #: 1 or 2\n\n') 
        fprintf('The eigenvectors are:\n')

        for i=1:manif.grow_info.dimension
            fprintf('(%i) Eigenvalue: %.3f, Eigenvector: (%.3f, %.3f, %.3f)\n',i, manif.grow_info.eigval(i), manif.grow_info.eigvec(i,1),manif.grow_info.eigvec(i,2),manif.grow_info.eigvec(i,3))
        end

        if ~isempty(eigenval_idx)
            x = eigenval_idx;
            fprintf('\nChoice: eigenval_idx = %i\n\n', eigenval_idx)
        else
            fprintf('\n Which eigenvector do you want to use?') 
            prompt = "\n... Press 1 or 2\n... Press 0 in case you want to exit \n\n";
            x = input(prompt);
        end

        if x == 0
            return
        else
            manif.grow_info.eigval=manif.grow_info.eigval(x);
            manif.grow_info.eigvec=manif.grow_info.eigvec(x,:);
        end

    end


%% Printing the general information of this run

    names=fieldnames(par);
    
    fprintf('\n----');
    for k=1:length(names)
        fprintf(' %s:%0.2f ',names{k},par.(names{k}));
    end
    fprintf('----');

    if strcmp(manif.orientability,'orientation-preserving')
        fprintf('\n----  Computing W%s(%s), %s branch up to arclength %.0f each ---- ',lower(opts.stability(1)),opts.per_orbit.name, branch1, opts.user_arclength);
    else
        fprintf('\n----  Computing W%s(%s), both branches up to arclength %.0f each ---- ',lower(opts.stability(1)),opts.per_orbit.name,opts.user_arclength);
    end
    
    fprintf('\n----  Total maximum iterations of the fundamental domain: %i ---- ',manif.grow_info.max_funditer);
    fprintf('\n\n')
    fprintf('----- Acc. Conditions ----- ')
    fprintf('\n| AlphaMax:%.2e',manif.grow_info.alphamax)
    fprintf('\n| DeltaAlphaMax:%.2e ',manif.grow_info.deltalphamax)
    fprintf('\n| DeltaMin:%.2e ' ,manif.grow_info.deltamin)
    fprintf('\n| DeltaMax:%.2e \n' ,manif.grow_info.deltamax)
    fprintf('---------- \n')
    
    
    %-------------------------------------------------------------------
    %-------------------------------------------------------------------
    % initializing the fields for the initial information
    manif.grow_info.runinf.rem_deltamin=0; %how many points are removed because of deltamin
    manif.grow_info.runinf.rem_nan=0; %how many points are removed because of NaN values
    manif.grow_info.runinf.rem_inf=0; %how many points are removed at infinity because of duplication
    manif.grow_info.runinf.add_alphamax=0; %how many points are added because of alpha_max
    manif.grow_info.runinf.add_deltamax=0;%how many points are added because of delta_max
    manif.grow_info.runinf.add_deltalphamax=0;  %how many points are added because of (delta alpha)_max

%% computes first fundamental domain

% by default, we choose the eigenvector going to positive x
if manif.grow_info.eigvec(1) < 0 % is it is going to negative x, then consider the other side
    manif.grow_info.eigvec = -manif.grow_info.eigvec;
end

% Choosing the initial distance correctly:
% at what distance from the fix point the fundamental domain starts so 
% the end of the fundamental domain is approx at distance: opts.accpar.init_step
% we assume the linearlize system to do the calculation

% distance scaled by the expansion/contraction-rate (eigenvalue)
    if strcmp(manif.stability,'Smanifold')
        dist = manif.grow_info.init_step*abs(manif.grow_info.eigval^mapiter);
    elseif strcmp(manif.stability,'Umanifold')
        dist = manif.grow_info.init_step/abs(manif.grow_info.eigval^mapiter);
    end

% starting point of the fundamental domain
fund_init.x= per_orbit.x(1)+dist*manif.grow_info.eigvec(1);
fund_init.y= per_orbit.y(1)+dist*manif.grow_info.eigvec(2);
fund_init.z= per_orbit.z(1)+dist*manif.grow_info.eigvec(3);

% mapping the initial point to obtain the first fundamental domain
fund_end = thesystem.mapping(fund_init,opts,sign*mapiter*period);


% interpolate N=4 points linearly
fund.points.x = linspace(fund_init.x, fund_end.x, 6);
fund.points.y = linspace(fund_init.y, fund_end.y, 6);
fund.points.z = linspace(fund_init.z, fund_end.z, 6);


fund.points.arc = arclength(fund.points);
fund.points.idx_fund_dom = [1 numel(fund.points.x)];


manif.points=cell(1,period);



% from the manifold associated to which periodic point the data comes from (which one is their "preimage")
if strcmp(manif.stability,'Umanifold') % if we are computing Ustable, we use f^1
    preimage_idx = [period 1:period-1]; %[k, 1,...,k-1] -> 1 comes from 2..., etc
    iteration_idx = 1:period; % [1, 2, ..., 6] -> 1 maps to 2....
elseif strcmp(manif.stability,'Smanifold')
    preimage_idx = [2:period 1]; %[2, 3,...,k, 1] -> 1 comes from 2..., etc
    iteration_idx = [1 period:-1:2]; % [1, 6, ..., 2] -> 1 maps to 6....
end 


if strcmp(manif.orientability,'orientation-preserving') 
    total_arc=zeros(1,period);
    % initialize the cells
    [manif.points{1:period}] = deal(struct('name',[], branch1, []));
    for k=1:period
        manif.points{k}.name = sprintf('%s_%s',name,num2str(k));
        manif.points{k}.(branch1) = struct('x', [], 'y', [], 'z', [], 'arc', [], 'idx_fund_dom', [], 'idx_preimages', []);
    end

    % initialize names
    for k=1:period
        manif.points{k}.(branch1).branch_preimage = branch1; %name of the branch
        manif.points{k}.(branch1).name_preimage = manif.points{preimage_idx(k)}.name; %index of the periodic orbit
    end

    manif.points{1}.(branch1).x=fund.points.x;
    manif.points{1}.(branch1).y=fund.points.y;
    manif.points{1}.(branch1).z=fund.points.z;
    manif.points{1}.(branch1).arc=fund.points.arc;
    manif.points{1}.(branch1).idx_fund_dom=fund.points.idx_fund_dom;
    manif.points{1}.(branch1).idx_preimages=zeros(1,numel(fund.points.x));

    total_arc(1) = fund.points.arc(end);
   

elseif strcmp(manif.orientability,'orientation-reversing') 
    total_arc_branch1=zeros(1,period);
    total_arc_branch2=zeros(1,period);

    % initialize the cells
    [manif.points{1:period}] = deal(struct('name',[], branch1, [], branch2, []));
    for k=1:period
        manif.points{k}.name = sprintf('%s_%s',name,num2str(k));
        manif.points{k}.(branch1) = struct('x', [], 'y', [], 'z', [], 'arc', [], 'idx_fund_dom', [], 'idx_preimages', []);
        manif.points{k}.(branch2) = struct('x', [], 'y', [], 'z', [], 'arc', [], 'idx_fund_dom', [], 'idx_preimages', []);
    end

    % initialize names
    for k=1:period
        if k==1
            manif.points{k}.(branch1).branch_preimage = branch2; %name of the branch
            manif.points{k}.(branch1).name_preimage = manif.points{preimage_idx(k)}.name; %index of the periodic orbit

            manif.points{k}.(branch2).branch_preimage = branch1; %name of the branch
            manif.points{k}.(branch2).name_preimage = manif.points{preimage_idx(k)}.name; %index of the periodic orbit
        else
            manif.points{k}.(branch1).branch_preimage = branch1; %name of the branch
            manif.points{k}.(branch1).name_preimage = manif.points{preimage_idx(k)}.name; %index of the periodic orbit

            manif.points{k}.(branch2).branch_preimage = branch2; %name of the branch
            manif.points{k}.(branch2).name_preimage = manif.points{preimage_idx(k)}.name; %index of the periodic orbit
        end
    end

    manif.points{1}.(branch1).x=fund.points.x;
    manif.points{1}.(branch1).y=fund.points.y;
    manif.points{1}.(branch1).z=fund.points.z;
    manif.points{1}.(branch1).arc=fund.points.arc;
    manif.points{1}.(branch1).idx_fund_dom=fund.points.idx_fund_dom;
    manif.points{1}.(branch1).idx_preimages=zeros(1,numel(fund.points.x));

    total_arc_branch1(1) = fund.points.arc(end);
end

% plot3(fund.points.x,fund.points.y,fund.points.z,'LineWidth',2)
% plot3(per_orbit.x,per_orbit.y,per_orbit.z,'*')
%%

fprintf(' Starting fundamental domain arclength  %f\n',fund.points.arc(end))

Manif=manif;

% progress table
varNames = arrayfun(@(i) sprintf('p%d_%d', period, i), 1:period, 'UniformOutput', false);
if strcmp(manif.orientability,'orientation-preserving')
    data = num2cell(zeros(1, period), 1);
    T = table(data{:},'VariableNames',varNames,'RowName',{sprintf('%s arc%f\n',branch1)});
else
    data = num2cell(zeros(2, period), 1);
    T = table(data{:},'VariableNames',varNames,'RowName',{sprintf('%s arc%f\n',branch1),sprintf('%s arc%f\n',branch2)});
end

T{1,1} = round(fund.points.arc(end), 2, 'significant');
disp(T)

%----------- Starting the loop
iter=0;
stop_arc=zeros(1,period); % flag to stop iterating if we reached the desired arclength
stop_arc_branch2 = zeros(1,period); % flag to stop iterating if we reached the desired arclength (neg branch)
stop_arc_branch1 = zeros(1,period); % flag to stop iterating if we reached the desired arclength (pos branch)

% iter: number of iterations of the fundamental domain
while iter < manif.grow_info.max_funditer && sum(stop_arc) ~= period && (sum(stop_arc_branch2) + sum(stop_arc_branch1)) ~= 2*period
       % we haven't reached the max number of iterations
       % and, we haven't reached the arclength (orientation-preserving case)
       % and, we haven't reached the arclength for the two branches of the manifolds (orientation-reversing case)
 
         iter=iter+1;
         idx_seg = iteration_idx(mod(iter,period)+1); %index if current segment of manifold. It is which segment we are computing from 1 to period.
          % For instance, if is iter=2, k=period and we are computing the stable manifold, then idx_seg = k (we are computing the segment asociated to pk_k)
        
          %Chech in which branch are we at
          if strcmp(manif.orientability,'orientation-preserving')
              branch=branch1;
          elseif strcmp(manif.orientability,'orientation-reversing')
            if mod(iter,2*period) < period
                branch=branch1;
            else
                branch=branch2;
            end
          end


         %mapping the points
         mappoints = thesystem.mapping(fund.points,opts,sign);

         idx_eps_preimages=fund.points.idx_fund_dom(1):fund.points.idx_fund_dom(2);
 

     %% STARTING THE ALGORITHM
 
    %% ----------- Removing points at infinity or NaN values. This is only because issues from compactification

    %Delete NaN values
    nan_idx=union(union(find(isnan(mappoints.x)),find(isnan(mappoints.y))),find(isnan(mappoints.z)));
    Manif.grow_info.runinf.rem_nan=Manif.grow_info.runinf.rem_nan+numel(nan_idx);

      if numel(nan_idx)==numel(mappoints.x)
        fprintf('\n\n ALL THE POINTS HAVE BEEN MAPPED TO NAN VALUES. LAST POINT IS AT (%f,%f,%f)\n\n',fund.points.x(end),fund.points.y(end),fund.points.z(end))

        break;
       end

    mappoints.x(nan_idx)=[]; 
    mappoints.y(nan_idx)=[]; 
    mappoints.z(nan_idx)=[]; 

    fund.points.x(nan_idx)=[]; 
    fund.points.y(nan_idx)=[]; 
    fund.points.z(nan_idx)=[]; 

    idx_eps_preimages(nan_idx)=[]; 


    %Delete points at infinity
    inf_idx=union(find((mappoints.x.^2+mappoints.y.^2)==1),find(abs(mappoints.z)==1));
    Manif.grow_info.runinf.rem_inf=Manif.grow_info.runinf.rem_inf+numel(inf_idx);



    if numel(inf_idx)==numel(mappoints.x)
        fprintf('\n\n THE MANIFOLD HAS REACHED INFINITY AT (%f,%f,%f)\n\n',fund.points.x(end),fund.points.y(end),fund.points.z(end))
        break;
    end



    mappoints.x(inf_idx)=[]; 
    mappoints.y(inf_idx)=[]; 
    mappoints.z(inf_idx)=[]; 

    fund.points.x(inf_idx)=[]; 
    fund.points.y(inf_idx)=[]; 
    fund.points.z(inf_idx)=[]; 

    idx_eps_preimages(inf_idx)=[];


     %% Replace first point of the current mapped segment by the last point of the previous segment (continuous manifold)
     % and monitoring the distance between first point of the mapped points and the last point of previous fundamental domain on that branch

    if numel(Manif.points.(branch).x) > 0 % there has to be a segment already on the branch to check for the distance
        dist = sqrt((Manif.points.(branch).x(end) - mappoints.x(1))^2 +...
                    (Manif.points.(branch).y(end) - mappoints.y(1))^2 +...
                    (Manif.points.(branch).z(end) - mappoints.z(1))^2);

         fprintf('dist %.e)\n', dist)

         %replace the first point of the mapped points by the last point of the previous segment
         mappoints.x(1) = Manif.points{idx_seg}.(branch).x(end);
         mappoints.y(1) = Manif.points{idx_seg}.(branch).y(end);
         mappoints.z(1) = Manif.points{idx_seg}.(branch).z(end);

         if dist > Manif.grow_info.deltamin
             fprintf('Warning! The distance between the first point of the current mapped segment and the last point of the previous segment exceeds Deltamin. Current distance: %.e %.e\n', dist)
         end
     end

    %%
    fprintf('\n ITERATION OF THE FUNDAMENTAL DOMAIN %i: ',iter)


    arc_mappoints = arclength(mappoints); %arclength of the mapped fundamental domain

    % Check how much arclength is left to compute
    if strcmp(manif.orientability,'orientation-preserving')
        needed_arc = opts.user_arclength - total_arc(idx_seg);
     
    elseif strcmp(manif.orientability,'orientation-reversing')
        if mod(iter,2*period) < period
            needed_arc = opts.user_arclength - total_arc_branch1(idx_seg);
        else
            needed_arc = opts.user_arclength - total_arc_branch2(idx_seg);
        end
    end

    fprintf('%s %s branch \n', manif.points{idx_seg}.name, branch);

    % if in this iteration we exceed the desired arclength, then we chop
    % the fundamental domain up to the desired arclength
    if needed_arc < 0
        fprintf('Warning! We have looped back to a branch that already have the desired arclength. The process will terminate... \n');
        return;
    end

    if arc_mappoints(end) > needed_arc
        
         idx_arc = find(arc_mappoints > needed_arc,1); %where we exceed the extra needed arc
         % chop the fund domain and the mappoints up to there
         
         fund.points.x=fund.points.x(1:idx_arc);
         fund.points.y=fund.points.y(1:idx_arc);
         fund.points.z=fund.points.z(1:idx_arc);
     
         mappoints.x=mappoints.x(1:idx_arc);
         mappoints.y=mappoints.y(1:idx_arc);
         mappoints.z=mappoints.z(1:idx_arc);
     
         idx_eps_preimages=idx_eps_preimages(1:idx_arc);
     
         if strcmp(manif.orientability,'orientation-preserving')
             stop_arc(idx_seg)=1; % this is the last iteration
         elseif strcmp(manif.orientability,'orientation-reversing')
             if mod(iter,2*period) < period % we are on pos
                 stop_arc_branch1(idx_seg)=1;
             else
                 stop_arc_branch2(idx_seg)=1;
             end
         end
    end

% % print the progression of which branch has reached the arclength. 0 false, 1 true
% if strcmp(manif.orientability,'orientation-preserving')
%     stop_arc
% elseif strcmp(manif.orientability,'orientation-reversing')
%     stop_arc_branch1
%     stop_arc_branch2
% end

%------------------------------------------------------- 
%---%----------- Adding points depending on Acc. Cond.  

    %initializing the structures to add points
    add_acc=struct(); 
    newpoints=struct();
    mapnewpoints=struct();


	add_acc.iter=0;
    add_acc.failed=[]; % points that failed acc cond in last loop
    add_acc.loop=true; % still doing the while loop 


    % Interpolation mesh
    fund_initial=fund.points; %starting mesh of the fundamental domain
    t_initial=0:1/(numel(fund_initial.x)-1):1; % parametrization for meshpoints

%---%-------------- Loop of the same mesh checking acc cond (this adds points)
    while add_acc.loop 
%---%--------------
        add_acc.loop=false; %to stop while loop % if at least one point is added this turns true
        add_acc.iter=add_acc.iter+1; %counter of iterations


        % Interpolating points from previous fundamental domain.
         if add_acc.iter==1
            tt=t_initial;
            t_interp=tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of interpolated points
            interp = makima3D(fund_initial,t_initial,t_interp); % compute interpolated preimage
            mapinterp = thesystem.mapping(interp,opts,sign); % interpolated image
         else
             tt=sort([tt t_interp(add_acc.add)]); %parametrization of (new) mesh points
             t_interp=tt(1:end-1)+(tt(2:end)-tt(1:end-1))/2; % parametrization of (new) interpolated points
             interp = makima3D(fund_initial,t_initial,t_interp); % compute interpolated preimage
             mapinterp = thesystem.mapping(interp,opts,sign); % interpolated image
         end

%%
        add_acc.add=[]; %points we are going to add

        % idx of points to check acc cond
        if add_acc.iter==1
            for_idx=2:(numel(mappoints.x)-1);
        else
            for_idx=add_acc.failed;
        end      

        fprintf('  loop number %i (points to check %i...',add_acc.iter,numel(for_idx));

%-------%---------- Going through the points that failed
        million=0;
        for k=for_idx
%-------%----------      
            % a flag for when # million points have been checked
            million=million+1;

            if floor(million/1000000)==ceil(million/1000000) %is integer?
                fprintf(' -checkpoint %i million points checked-...',floor(million/1000000));
            end


            % coordinates of mapped points
            add_acc.p0=[mappoints.x(k-1), mappoints.y(k-1), mappoints.z(k-1)];
            add_acc.p1=[mappoints.x(k), mappoints.y(k), mappoints.z(k)]; % the point we are actually looking at
            add_acc.p2=[mappoints.x(k+1), mappoints.y(k+1), mappoints.z(k+1)];

            % Distance btw points
            add_acc.delta0=norm(add_acc.p1-add_acc.p0); % before
            add_acc.delta2=norm(add_acc.p1-add_acc.p2); % after
            add_acc.alpha = angles(add_acc.p0,add_acc.p1,add_acc.p2); % angle btw points

            % points btw p0p1 and p1p2 in the interpolated points
            add_acc.p0_new=[mapinterp.x(k-1), mapinterp.y(k-1), mapinterp.z(k-1)];
            add_acc.p2_new=[mapinterp.x(k), mapinterp.y(k), mapinterp.z(k)];


%-----------%------ Adding points


%-----------%------ If delta > deltamax 

            %------ If it is the second point in the mesh
            %------ Check the first delta and add a point before
            if k==2 && add_acc.delta0>manif.grow_info.deltamax  

                %add point p01
                add_acc.add    =[add_acc.add k-1]; %idx of the point we are going to add
                add_acc.loop=true; % We have to check if we need to put more points in the mesh
                Manif.grow_info.runinf.add_deltamax=Manif.grow_info.runinf.add_deltamax+1;

            %------ Check the second delta and add a point after 
            elseif add_acc.delta2>manif.grow_info.deltamax %k>2 && add_acc.delta2>manif.grow_info.deltamax
                %add point p12
                add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                add_acc.loop=true; %idx of the point we are going to add
                Manif.grow_info.runinf.add_deltamax=Manif.grow_info.runinf.add_deltamax+1;



%-----------%------ If alpha > alphamax  or   Delta*alpha > Delta*alpha max

            %------ If alpha fails or BOTH Delta*alpha fail
            %------ Choose where to add a point.
            elseif add_acc.alpha>=manif.grow_info.alphamax || (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) %1 

            %------ Add point where Delta>Deltamin

                %-- If only Delta0>Deltamin
                %-- Add a point btw p0 and p1
                if add_acc.delta0>manif.grow_info.deltamin && add_acc.delta2<manif.grow_info.deltamin 
                    % Add a point if we didnt added the point in the previous acc cond checks
                    if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                        add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                        add_acc.loop=true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                            Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                        else
                            Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                        end
                    end
                end

                %-- If only Delta2>Deltamin
                %-- Add a point btw p1 and p2
                if add_acc.delta2>manif.grow_info.deltamin && add_acc.delta0<manif.grow_info.deltamin
                    %add point p12
                    add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                    add_acc.loop=true; %idx of the point we are going to add
                    if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                        Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                    else
                        Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                    end
                end

                %-- If both Delta0 and Delta 2 > Deltamin
                %-- Choose where to add point
                if add_acc.delta2>manif.grow_info.deltamin && add_acc.delta0>manif.grow_info.deltamin
                    add_acc.alpha0_new = angles(add_acc.p0_new,add_acc.p1,add_acc.p2); % angle btw points
                    add_acc.alpha2_new = angles(add_acc.p0,add_acc.p1,add_acc.p2_new); 

                    if add_acc.alpha0_new < add_acc.alpha2_new 
                        % Add a point if we didnt added the point in the previous acc cond checks
                        if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                            add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                            add_acc.loop=true; %idx of the point we are going to add
                            if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                                Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                            else
                                Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                            end
                        end
                    else
                        %add point p12
                        add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                        add_acc.loop=true; %idx of the point we are going to add
                        if (add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax) 
                            Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                        else
                            Manif.grow_info.runinf.add_alphamax=Manif.grow_info.runinf.add_alphamax+1;
                        end
                    end
                end



%-----------%------ If just one Delta*alpha > Delta*alpha max ( and alpha < alphamax (previous elseif is when alpha > alphamax )

            %------ Delta0*alpha > max, and Delta0 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta0*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta0>manif.grow_info.deltamin
                % Add a point if we didnt added the point in the previous acc cond checks
                if numel(add_acc.add)==0 || add_acc.add(end)~=k-1
                    add_acc.add    =[add_acc.add k-1];  %idx of the point we are going to add
                    add_acc.loop=true; %idx of the point we are going to add
                    Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;
                end

            %------ Delta2*alpha > max, and Delta2 > Deltamin
            %------ Add a point btw p0 and p1
            elseif add_acc.delta2*add_acc.alpha>=manif.grow_info.deltalphamax && add_acc.delta2>manif.grow_info.deltamin
                %add point p12
                add_acc.add    =[add_acc.add k]; %idx of the point we are going to add
                add_acc.loop=true; %idx of the point we are going to add
                Manif.grow_info.runinf.add_deltalphamax=Manif.grow_info.runinf.add_deltalphamax+1;


%-----------%------              
            end    % (if loop) Adding points     
%-----------%------      
%-------%---------- 
        end       % (for loop) Going through the points that failed 
%-------%----------

        newpoints.x=interp.x(add_acc.add);
        newpoints.y=interp.y(add_acc.add);
        newpoints.z=interp.z(add_acc.add);


        mapnewpoints.x=mapinterp.x(add_acc.add);
        mapnewpoints.y=mapinterp.y(add_acc.add);
        mapnewpoints.z=mapinterp.z(add_acc.add);

        newidx_preimages=idx_eps_preimages(add_acc.add);

       fprintf(' added points: %i) \n',numel(add_acc.add));


        % Add the points 
        if ~isempty(add_acc.add)

            % get updated idx of failed points
            plus=0:(numel(add_acc.add)-1);
            add_acc.failed = unique(sort([add_acc.add+plus add_acc.add+plus+1 add_acc.add+plus+2]));
            add_acc.failed=add_acc.failed(add_acc.failed>1);
            add_acc.failed=add_acc.failed(add_acc.failed<numel(mappoints.x)+numel(add_acc.add));

            % add points in the mapped manifold and in the old manifold
            mappoints.x=insert(mappoints.x,mapnewpoints.x,add_acc.add);
            mappoints.y=insert(mappoints.y,mapnewpoints.y,add_acc.add);
            mappoints.z=insert(mappoints.z,mapnewpoints.z,add_acc.add); 

            fund.points.x=insert(fund.points.x,newpoints.x,add_acc.add);
            fund.points.y=insert(fund.points.y,newpoints.y,add_acc.add);
            fund.points.z=insert(fund.points.z,newpoints.z,add_acc.add); 

            idx_eps_preimages=insert(idx_eps_preimages,newidx_preimages,add_acc.add);

        end

if add_acc.iter>30
    add_acc.loop=false;
end

%---%--------------       
    end           % (while loop) Checking acc cond 
%---%--------------

%---%----------- Final section: save info

    fund.points.x=mappoints.x;
    fund.points.y=mappoints.y;
    fund.points.z=mappoints.z;
    fund.points.arc = arclength(fund.points);

    %chop the manifold up to the desired arclength
    if fund.points.arc(end) > needed_arc
         idx_arc = find(fund.points.arc < needed_arc, 1, 'last'); %find last point that is less than the needed arc

         % chop the fund domain and the mappoints up to the needed arclength
         fund.points.x = spline(fund.points.arc,fund.points.x,[fund.points.arc(1:idx_arc) needed_arc]);
         fund.points.y = spline(fund.points.arc,fund.points.y,[fund.points.arc(1:idx_arc) needed_arc]);
         fund.points.z = spline(fund.points.arc,fund.points.z,[fund.points.arc(1:idx_arc) needed_arc]);
         fund.points.arc = arclength(fund.points);
    end

    fprintf(' Fund domain arclength %.2f \n', fund.points.arc(end));




    % obtain the starting index of the fund domain to store. start=1 if is
    % the first segment store on that branch, start = 2
    if numel(Manif.points{idx_seg}.(branch).x) == 0 % this is not the first time storing data on that branch
        start=1;
    else
        start=2;
    end

    iter_fund=numel(Manif.points{idx_seg}.(branch).idx_fund_dom)/2 + 1;

    % Add new segment to the entire manifold
    N = numel(Manif.points{idx_seg}.(branch).x);
    Manif.points{idx_seg}.(branch).x=[Manif.points{idx_seg}.(branch).x fund.points.x(start:end)];
    Manif.points{idx_seg}.(branch).y=[Manif.points{idx_seg}.(branch).y fund.points.y(start:end)];
    Manif.points{idx_seg}.(branch).z=[Manif.points{idx_seg}.(branch).z fund.points.z(start:end)];
    Manif.points{idx_seg}.(branch).arc = arclength(Manif.points{idx_seg}.(branch));

    Manif.points{idx_seg}.(branch).idx_preimages=[Manif.points{idx_seg}.(branch).idx_preimages idx_eps_preimages(1:numel(fund.points.x(start:end)))];
    Manif.points{idx_seg}.(branch).idx_fund_dom(iter_fund,:)=[max(N,1) numel(Manif.points{idx_seg}.(branch).x)];
    fund.points.idx_fund_dom = Manif.points{idx_seg}.(branch).idx_fund_dom(iter_fund,:);

    
    %Check once again if we reached the desired arclength
    if strcmp(manif.orientability,'orientation-preserving') 

        total_arc(idx_seg) = Manif.points{idx_seg}.(branch).arc(end);
        fprintf(' Total arclength of manifold of %s %s: %.2f \n\n', Manif.points{idx_seg}.name, branch, total_arc(idx_seg));
        T{1,idx_seg} = round(total_arc(idx_seg), 2, 'significant');
        disp(T)
        %Check once again if we reached the desired arclength
        if total_arc(idx_seg) > opts.user_arclength 
            stop_arc(idx_seg) = 1;
        end
        
    elseif strcmp(manif.orientability,'orientation-reversing')
        if mod(iter,2*period)<period

            total_arc_branch1(idx_seg) = Manif.points{idx_seg}.(branch).arc(end);
            fprintf(' Total arclength of manifold of %s %s: %.2f \n\n', Manif.points{idx_seg}.name, branch, total_arc_branch1(idx_seg));
            T{1,idx_seg} = round(total_arc_branch1(idx_seg), 2, 'significant');
            disp(T)
            %Check once again if we reached the desired arclength
            if total_arc_branch1(idx_seg) > opts.user_arclength 
                stop_arc_branch1(idx_seg) = 1;
            end

        else

            total_arc_branch2(idx_seg) = Manif.points{idx_seg}.(branch).arc(end);
            fprintf(' Total arclength of manifold of %s %s: %.2f \n\n',Manif.points{idx_seg}.name, branch, total_arc_branch2(idx_seg));
            T{2,idx_seg} = round(total_arc_branch2(idx_seg), 2, 'significant');
            disp(T)
            %Check once again if we reached the desired arclength
            if total_arc_branch2(idx_seg) > opts.user_arclength 
                stop_arc_branch2(idx_seg) = 1;
            end

        end
    end
% plot3(fund.points.x,fund.points.y,fund.points.z,'LineWidth',2)
%---%----------- END Final section: save info
end

% for k=1:period
%     text(Manif.per_orbit.coord_compactified.x(k),Manif.per_orbit.coord_compactified.y(k),Manif.per_orbit.coord_compactified.z(k),sprintf('%s', strrep(Manif.points{k}.name,'_','\_')))
% end
% xlabel('x')
% ylabel('y')
% zlabel('z')


%--------------- END adding points
% 
% Erase last fundamental domain if the computation stopped chopping the last part of the manifold
if strcmp(manif.orientability,'orientation-preserving') && sum(stop_arc) == period
    for k=1:period
        Manif.points{k}.(branch1).idx_fund_dom(end, :) = [];
    end
elseif strcmp(manif.orientability,'orientation-reversing') && sum(stop_arc_branch1) + sum(stop_arc_branch2) == 2*period
    for k=1:period
        Manif.points{k}.(branch1).idx_fund_dom(end, :) = [];
        Manif.points{k}.(branch2).idx_fund_dom(end, :) = [];
    end
end

 
Manif.grow_info.runinf.time=toc;

fprintf('\n elapsed time is %.3f seconds\n\n',Manif.grow_info.runinf.time)
fprintf('   * %i points removed \n',Manif.grow_info.runinf.rem_deltamin+Manif.grow_info.runinf.rem_nan+Manif.grow_info.runinf.rem_inf) 
fprintf('   * %i points added from deltamax \n',Manif.grow_info.runinf.add_deltamax) 
fprintf('   * %i points added from alpha \n',Manif.grow_info.runinf.add_alphamax) 
fprintf('   * %i points added from delta*alpha \n',Manif.grow_info.runinf.add_deltalphamax) 
%% FUNCTIONS

function interp = makima3D(points,t,tt)
% get interpolation points
% t: parametrization of the points
% tt: parametrization of the interpolated points
interp=struct();

interp.x = interp1(t,points.x,tt,'makima','extrap');
interp.y = interp1(t,points.y,tt,'makima','extrap');
interp.z = interp1(t,points.z,tt,'makima','extrap');
end 

%----------------

function arclen = arclength(points)
%arclength between each point of a vector (px,py,pz)
arclen=((points.x(1:end-1)-points.x(2:end)).^2 + (points.y(1:end-1)-points.y(2:end)).^2 + (points.z(1:end-1)-points.z(2:end)).^2).^(1/2);
arclen=[0 cumsum(arclen)];
end % function arclength

%----------------

function alpha=angles(p0,p1,p2)
%angle between p0p1 and p1p2
n1 = (p1 - p2) / norm(p1 - p2);  % Normalized vectors
n2 = (p0 - p1) / norm(p0 - p1);
alpha = atan2(norm(cross(n1, n2)), dot(n1, n2)); %gives value from 0 to pi
end

%----------------

function Anew=insert(A,B,ind)
% Anew: new vector with the new values
% A: Old vector
% B: vector with new values
%ind: index to insert values after this row


    % Preallocate output
    Anew = zeros(1,numel(A)+numel(B));

    % Find indices for old data
    addRows = ismember(1:numel(A), ind);
    oldDataInd = (1:numel(A)) + cumsum([0, addRows(1:end-1)]);

    % Add in old data
    Anew(oldDataInd) = A;

    % Find indices for new data
    newDataInd = (1:length(ind)) + ind;

    % Add in new data
    Anew(newDataInd) = B;
end

end