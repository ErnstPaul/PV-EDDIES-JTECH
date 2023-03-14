function [maxMat,maxRet] = Threshold_eddy_centers(FIELD,STATS,X,Y,FACTOR)
%  Potential eddy centers are saved/updated as the structure array
%  {center(t)} in [path_out,'eddy_centers_',runname'] with followings
%  fields (type and coordinates):
%  - centers(t).step : step when the eddy was detected
%  - centers(t).type(n) : eddy type (1 => cyclonic; -1 => anticyclonic)
%  - centers(t).X(n) : eddy center X coordinate
%  - centers(t).Y(n) : eddy center Y coordinate
%  - centers(t).i(n) : eddy center column index
%  - centers(t).j(n) : eddy center row index
%=========================

%This is gridded
grid_ll = 1;

%Set mask
mask = ~isnan(FIELD);

%----------------------------------------
% initialise centers as structure
centers0 = struct('step',nan,'type',[],'X',[],'Y',[],'i',[],'j',[]);

%---------------------------------------------
% Max FIELD 'centers0' for a given step k above lat_min (basically 5°)
%---------------------------------------------

%---------------------------------------------
% Define contour, including potential centers
LOW = abs(FIELD);
LOW(isnan(FIELD)) = 0;

% Create regular grid
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end

%The contour factor is calculated here
K = STATS(1,2)*FACTOR;

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

% scan each FIELD contour
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
        Lm = FIELD;
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
                        centers0.type(k) = sign(LC(1)-STATS(1,1));
                        centers0.X(k)    = xLmax;
                        centers0.Y(k)    = yLmax;
                        [centers0.j(k),centers0.i(k)] = find(Lm==LC(1),1);

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

%% Set up return matrices
%Loop on struct
f1 = []; f1count = 1;
f2 = []; f2count = 1;
%Loop on found centers
for i = 1:length(centers0.i)
    thisX = centers0.j(i);
    thisY = centers0.i(i);
    thisInd = sub2ind(size(FIELD),thisX,thisY);
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
% figure; pcolor(Xgrid,Ygrid,FIELD); hold on; shading interp; contour(Xgrid,Ygrid,LOW,[K K], 'k'); hold on; scatter(centers0.X, centers0.Y, '+m');
maxMat=sort([f1(:);f2(:)]);
maxRet={f1,f2};
