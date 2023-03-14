function [maxMat,maxRet] = AMEDA_eddy_centers(LOW,LNAM,X,Y,U,V,K)
%[centers0,centers] = mod_eddy_centers(source,stp,fields)
%
%  Detect the potential eddy centers present in the domain,
%  for each time step of the time serie of the 2-D velocity fields
%  defined by U and V, using the LNAM and Local Okubo-Weiss fields 
%  calculated by mod_fields.m and stored as {detection_fields(t)} in
%  [path_out,'fields_inter',runname]:
%
% - 'source' is the type of netcdf data (AVISO, NEMO, ROMS,...) that
%   determine which load_field_'source'is  used in the routine
% - 'fields' is the step 'stp' of the detection_fields computed with
%   mod_eddy_fields.m
%
% This routine will use:
%  - K: the abs(LNAM(LOW<0)) threshold to delimit the contour of
%       the potential eddy centers (one per contour)
%  - bx: used to define smaller areas around each center and find which
%       ones are included (and alone) in a closed streamline. Must take into
%       account the resolution factor!
%  - DH: delta isocontour of ssh
%
%  For a description of the input parameters see mod_eddy_param.m.

%
%  There are two step to select potential center:
%
%  First the max(|LNAM(LOW<0)>K|) are computed.
%
%  Then the potential centers saved are the max LNAM surrounded by at least
%  two closed contour of ssh (or psi by using compute_psi)
%  and have a minimal and a maximal size defined by nRmin and nR_lim.
%
%  Potential eddy centers are saved/updated as the structure array
%  {center(t)} in [path_out,'eddy_centers_',runname'] with followings
%  fields (type and coordinates):
%  - centers(t).step : step when the eddy was detected
%  - centers(t).type(n) : eddy type (1 => cyclonic; -1 => anticyclonic)
%  - centers(t).X(n) : eddy center X coordinate
%  - centers(t).Y(n) : eddy center Y coordinate
%  - centers(t).i(n) : eddy center column index
%  - centers(t).j(n) : eddy center row index
%
%
%  Max of |LNAM(LOW<0)>K| with the same fields as {centers} are
%  also saved in {centers0}
%  
%  (t is the time step index; n is the indice of eddy detected at t)
%
%-------------------------
%   Ver. 3.2 Apr 2015 Briac Le Vu
%   Ver. 3.1 2014 LMD from Nencioli et al. routines
%-------------------------
%
%=========================

%Set K param
Rmin = 25; Rmax = 1000;

%This is gridded
grid_ll = 1;

%Set mask
mask = ~isnan(U);

% to calculate psi extrapole U and V to 0 in the land
U(isnan(U)) = 0;
V(isnan(V)) = 0;

%----------------------------------------
% initialise centers as structure
centers0 = struct('step',nan,'type',[],'X',[],'Y',[],'i',[],'j',[]);
centers = centers0;

%---------------------------------------------
% Max LNAM 'centers0' for a given step k above lat_min (basically 5°)
%---------------------------------------------

%---------------------------------------------
% LNAM n OW criteria define contour including potential centers
OW = LOW; % Okubo Weiss
LOW = abs(LNAM);
LOW(OW<=0 | isnan(OW)) = 0;

%Calculate f
f_i = zeros(size(OW));
for latLoop = 1:(length(Y))
    % Calculate coriolis parameter across the entire grid
    f_i(:,latLoop) = gsw_f(Y(latLoop));
end

% Create regular grid
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%Get dX from lat and lon
disty = sw_dist( Xgrid(1:2,1), Xgrid(1:2,1), 'km');
distx = nan*Xgrid;
for i=2:size(Xgrid,1)-1
    for j=2:size(Xgrid,2)-1
        distx(i,j) = sw_dist( Ygrid(i,[j-1 j+1]), Xgrid(i,[j-1 j+1]), 'km')/2;
    end
end
Dxi = ( distx + disty )/2;
Dxi(1,:)  = Dxi(2,:)-diff(Dxi(2:3,:));
Dxi(end,:)= Dxi(end-1,:)+diff(Dxi(end-2:end-1,:));
Dxi(:,1)  = Dxi(:,2);
Dxi(:,end)= Dxi(:,end-1);

%get contours
HF = figure('visible','off');
CS = contour(Xgrid,Ygrid,LOW,[K K]);
close(HF)
X = Xgrid;
Y = Ygrid;


% Initialization
j = 1; % first coordinates of the contour scan
k = 1; % first contour
n_min = 4; %minimum 10 grid points (4 inside)
lat_min = -10;

% scan each LNAM contour
while j < size(CS,2)

    n = CS(2,j); % number of coordinates for the contour(j)
    xv = CS(1,j+1:j+n); % X values serie for the contour(j) coordinates
    yv = CS(2,j+1:j+n); % Y values serie for the contour(j) coordinates

    % validate only bigger contour
    if n >= n_min

        % make a mask outside the contour
        in = InPolygon(X,Y,xv,yv);
        maskin = mask;
        maskin(~in) = 0;
        maskin(in)  = 1;
        %Lm = LNAM.*maskin;
        Lm = LNAM;
        Lm(~in) = NaN;

        % L maximum value inside the contour(j)  
        if any(mask(:).*maskin(:)>0) && max(abs(Lm(:)))~=0

            LC = Lm(abs(Lm)==max(abs(Lm(:))));

            % save coordinates of the L maximum inside the contour(j)
            if mask(find(Lm==LC(1),1))==1
                
                xLmax = X(find(Lm==LC(1),1));
                yLmax = Y(find(Lm==LC(1),1));

                if ~grid_ll || (grid_ll && abs(yLmax) > lat_min)
                    
                    if ~any(centers0.X==xLmax & centers0.Y==yLmax)
                        centers0.type(k) = sign(LC(1));
                        centers0.X(k)    = xLmax;
                        centers0.Y(k)    = yLmax;
                        [centers0.j(k),centers0.i(k)] = find(Lm==LC(1),1);
                        MaxL(k) = LC(1);

                        % increment the counter
                        k = k + 1; % next contour
                    end
                end
            end
        end
    end

    % increment the counter
    j = j + n + 1; % series of coordinates of the next contour 
end
% 
% disp(['  -> ',num2str(k-1),' max LNAM found step ',num2str(stp)])
% 
% if k==1
%     disp(['!!! WARNING !!! No LNAM extrema found - check the LNAM computation step ',num2str(stp)])
% end
%% Set up return matrices
%Loop on struct
f1 = []; f1count = 1;
f2 = []; f2count = 1;
%Loop on found centers
for i = 1:length(centers0.i)
    thisX = centers0.j(i);
    thisY = centers0.i(i);
    thisInd = sub2ind(size(LNAM),thisX,thisY);
    %Slot into max or min category
    if (centers0.type(i) >= 0)
        f1(f1count) = thisInd;
        f1count = f1count + 1;
    else
        f2(f2count) = thisInd;
        f2count = f2count + 1;
    end
end
%DEBUG PLOT
% figure; pcolor(Xgrid,Ygrid,LOW); hold on; shading interp; contour(Xgrid,Ygrid,LOW,[K K], 'k'); hold on; scatter(centers0.X, centers0.Y, '+m');
maxMat=sort([f1(:);f2(:)]);
maxRet={f1,f2};
