%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%             Cette fonction calcule des statistiques                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Originally created, unchanged from Cori Pegliasco's code

% Compute the trajectory statistics used in the Cost Function
% Compute the differences of Amplitude, Radius, EKE between each
% timestep, and the distance between the consecutive eddy centers.
% The trajectories build on the first 100 days seems adequate.

% Inputs:
%  - path to '/EXTRACTION/EDDY_PROPERTIES/AE_Properties.mat' and
% '/EXTRACTION/EDDY_PROPERTIES/CE_Properties.mat'
%  _ dist_lim: Maximum distance (in degree) of the closest center
% Outputs:
%  - mean and std relative to the differences of Amplitude, Radius,
% EKE between each timestep, and the distance between the consecutive
% eddy centers.
%
% Last modification June 2016
function Paul_OPT_Step4_SS(inputDir,propDir,depthlevel)

load('DefaultData_OPT.mat') %See Paul_FSS_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
%months as a list of strings, yearsS as a list of strings and yearsS
addpath(strcat(basepath, 'FUNCTIONS'));

% Input directories and configfile paths.
%input_dir_interim = strcat(basepath, '/EXTRACTION/INDIVIDUAL_FILES/NEM_DLEV_');
input_dir_config = inputDir;
close all
clc
input_dir= propDir; %EDDY_PROPERTIES'];
input_dir= strcat(input_dir, num2str(depthlevel));

%% Statistics on AEs trajectories

filename = [input_dir '/AE_Properties.mat'];
load(filename)
filename_cont = [input_dir '/Eddy_Contours.mat'];
load(filename_cont,'AEs')
ndays = length(Nanti);
dist_lim = 2;
NNN = [10 100 500 1000:1000:20000 25000:5000:50000];
Nestimate_lim = 10000;

%% Trajectories construction
% Reorganize properties by decreasing amplitude:
% Reorganize eddy properties by decreasing amplitude and store
% the resulting indices in the variable
% (Stronger Amplitude of Amplitudd(IND)= First line of non-NAN data ))
[~,IND] = sort(Amplitude,1,'descend');
listvar = whos; % Create a structure list of variables name
[nr,nc] = size(Amplitude); % Store dimensions of a "full" variable

% Loop to reorganize each variable following the decreasing amplitude order
% specified by IND.
for k = 1:nc
    for i=1:length(listvar)
        N = listvar(i).size;
        if N(1)==nr && N(2)==nc  && strcmp(listvar(i).name,'IND')~=1 % To reorganize variable of nr,nc size but not IND
            eval(['' listvar(i).name '(:,k,:)=' listvar(i).name '(IND(:,k),k,:);']);
        end
    end
end

% Initialization:
NN = isfinite(Xcenter);         % Matrix = 1 when there is an eddy
numtraj = 0;			% ID of the trajectories
idtraj = NaN(nr,nc);            % Matrix with the associated trajectory id of every eddy evolved in a trajectory
Anticyclonic_Trajectories = NaN(1e4,2e3,6); % array to store Amplitude, Radius, EKE, Xcenter and Ycenter
date_num = double(date_num);
h = waitbar(0,'Part I.Anticyclones');
warning off

Nestimate = 0;
for j=1:ndays-1 % Want to get eddies annually evenly distributed 
    disp(int2str(j))
    waitbar(j/(ndays-1),h)
    idx = find(NN(:,j)==1 & isnan(idtraj(:,j))==1)'; % Eddies non already detected
    disp([ 'Map: ' datestr(date_num(j)) ])
    if Nestimate>Nestimate_lim
        break
    end
    for i=idx % Loop on each eddy (from biggest to smallest)
        numtraj = numtraj + 1; 	% Each new eddy become a trajectory
        if Nestimate>Nestimate_lim
            break
        end
        % Counter on time (indy=collumn position) and eddy studied (indx= row position)
        indy=j; 	% time (collumn position)
        indx=i;		% eddy studied (row position)
        idtraj(i,j) = numtraj; 		% Fill the matrix of trajectory id
        %disp(['Nestimate:   ' Nestimate  ' idx:  ' indx])
        lifetime = 1;  	% Initialization of the eddy lifetime
        
        % Fill the trajectory parameters whith the first eddy characteristics
        Anticyclonic_Trajectories(numtraj,lifetime,1) = Amplitude(i,j);
        Anticyclonic_Trajectories(numtraj,lifetime,2) = Radius(i,j);
        Anticyclonic_Trajectories(numtraj,lifetime,3) = EKE(i,j);
        Anticyclonic_Trajectories(numtraj,lifetime,4) = 0;
        Anticyclonic_Trajectories(numtraj,lifetime,5) = Xcenter(i,j);
        Anticyclonic_Trajectories(numtraj,lifetime,6) = Ycenter(i,j);
        
        % While eddy = 1, each eddy found takes part of the trajectory
        % A new eddy is defined as the first eddy of a trajectory
        eddy = 1;
        
        % Load contour and center position of detected eddy
        Xcon0 =double(AEs{indx,indy,1});
        Ycon0 =double(AEs{indx,indy,2});
        Xc = Xcenter(indx,indy);
        Yc = Ycenter(indx,indy);
        
        % search for the closest eddies on the next map
        [next_map_eddies,tot_dist,found_centroid_at] = cori_find_eddy_next_map(Xcentroid,Ycentroid,indx,indy,1,dist_lim,nc);
        
        if isempty(next_map_eddies)==0 % if a close eddy is found
            % testing the contours intersections
            [matching_contour,area_intersection] = cori_matching_contours(AEs,Area,indx,indy,1,next_map_eddies,tot_dist);
            if isempty(matching_contour)==0
                indd = matching_contour;
                lifetime = lifetime +1;
                Nestimate = Nestimate + 1;
                if Nestimate>Nestimate_lim
                    break
                end
                %disp(['Line 145; Nestimate:   ' int2str(Nestimate)  ' idx:  ' int2str(indx)])
                
                if length(matching_contour)==1 % only one contour
                elseif  length(matching_contour) >1 % more than 1 contour -> cost function
                    
                    cf_radius=((Radius(indd,indy+1)-Radius(indx,indy))/Radius(indx,indy)).^2;
                    cf_eke=((EKE(indd,indy+1)-EKE(indx,indy))/EKE(indx,indy)).^2;
                    cf_amp=((Amplitude(indd,indy+1)-Amplitude(indx,indy))/Amplitude(indx,indy)).^2;
                    
                    cf =  cf_radius + cf_eke + cf_amp;
                    [~,ii] = min(cf); % Find the the minimum of the cost function
                    indd = indd(ii);  % Define the corresponding eddy as the same eddy at the previous time
                end
                
                % Fill the trajectory
                Anticyclonic_Trajectories(numtraj,lifetime,1)=Amplitude(indd,indy+1);
                Anticyclonic_Trajectories(numtraj,lifetime,2)=Radius(indd,indy+1);
                Anticyclonic_Trajectories(numtraj,lifetime,3)=EKE(indd,indy+1);
                Anticyclonic_Trajectories(numtraj,lifetime,4)=ac_distance(Yc,Xc,Ycenter(indd,indy+1),Xcenter(indd,indy+1));
                Anticyclonic_Trajectories(numtraj,lifetime,5)=Xcenter(indd,indy+1);
                Anticyclonic_Trajectories(numtraj,lifetime,6)=Ycenter(indd,indy+1);
                
                % If this eddy was not detected as taking part of a trajectory before
                if isnan(idtraj(indd,indy+1))==1
                    
                    % Store the ID of the trajectory associated with this eddy
                    idtraj(indd,indy+1) = numtraj;
                    
                end
                % Go to next map
                indy = indy+1;
                indx = indd;
                
                if indy ==  nc % The last map is reached
                    %disp('Trajectory might continue longer than the period')
                    eddy = 0;
                end
                
            else % No Eddy -> trajectory stops -> break while
                %eddy = 0;
                %disp('No eddy intersection')
            end
            
        else % NO center near enough -> trajectory stops
            %eddy = 0;
            %disp('Too far')
        end
    end
end
close(h);

% Statistics

x = Anticyclonic_Trajectories(1:numtraj,:,5)'; x = x(:);
y = Anticyclonic_Trajectories(1:numtraj,:,6)'; y = y(:);
delt_dist_center = ac_distance(y(2:end),x(2:end),y(1:end-1),x(1:end-1));
delt_Ampl = diff(Anticyclonic_Trajectories(1:numtraj,:,1),1,2);
delt_Rad = diff(Anticyclonic_Trajectories(1:numtraj,:,2),1,2);
delt_EKE = diff(Anticyclonic_Trajectories(1:numtraj,:,3),1,2);

mean_delt_Ampl = nanmean(delt_Ampl(:));
std_delt_Ampl = nanstd(delt_Ampl(:));

mean_delt_Rad = nanmean(delt_Rad(:));
std_delt_Rad = nanstd(delt_Rad(:));

mean_delt_EKE = nanmean(delt_EKE(:));
std_delt_EKE = nanstd(delt_EKE(:));

mean_delt_dist = nanmean(delt_dist_center(:));
std_delt_dist = nanstd(delt_dist_center(:));

% Save data in previous matrices
eval([' save -append ' filename ' Metadata mean* std* numtraj'])
keep input_dir ndays  dist_lim  Nestimate_lim

%% Statistics on CEs trajectories
% Open Anticyclonic matrix
filename = [input_dir '/CE_Properties.mat'];
load(filename)

% Open contours matrices
filename_cont = [input_dir '/Eddy_Contours.mat'];
load(filename_cont,'CEs')

%% Trajectories construction
% Reorganize properties by decreasing amplitude:
% Reorganize eddy properties by decreasing amplitude and store
% the resulting indices in the variable
% (Stronger Amplitude of Amplitudd(IND)= First line of non-NAN data ))
[~,IND] = sort(Amplitude,1,'descend');
listvar = whos; % Create a structure list of variables name
[nr,nc] = size(Amplitude); % Store dimensions of a "full" variable


% Loop to reorganize each variable following the decreasing amplitude order
% specified by IND.
for k = 1:nc
    for i=1:length(listvar)
        N = listvar(i).size;
        if N(1)==nr && N(2)==nc  && strcmp(listvar(i).name,'IND')~=1 % To reorganize variable of nr,nc size but not IND
            eval(['' listvar(i).name '(:,k,:)=' listvar(i).name '(IND(:,k),k,:);']);
        end
    end
end

% Initialization:
NN = isfinite(Xcenter);         % Matrix = 1 when there is an eddy
numtraj = 0;			% ID of the trajectories
idtraj = NaN(nr,nc);            % Matrix with the associated trajectory id of every eddy evolved in a trajectory
Cyclonic_Trajectories = NaN(1e4,2e3,6); % array to store Amplitude, Radius, EKE, Xcenter and Ycenter
date_num = double(date_num);
h = waitbar(0,'Part II.Cyclones');

Nestimate = 0;
Nestimate_lim=ndays-1;
% I think this is the limit of max number of eddies used. Want to get all
% eddies, so if use this line in the future, make Nestimate_lim VERY BIG
warning off
for j=1:ndays-1
    idx = find(NN(:,j)==1 & isnan(idtraj(:,j))==1)'; % Eddies non already detected
    disp([ 'Map: ' datestr(date_num(j)) ])
    if Nestimate>Nestimate_lim
        break
    end
    checklim=Nestimate-Nestimate_lim;
    disp(int2str(checklim))
    for i=idx % Loop on each eddy (from biggest to smallest)
        if Nestimate>Nestimate_lim
            break
        end
        numtraj = numtraj + 1; 	% Each new eddy become a trajectory
        % Counter on time (indy=collumn position) and eddy studied (indx= row position)
        indy=j; 	% time (collumn position)
        indx=i;		% eddy studied (row position)
        idtraj(i,j) = numtraj; 		% Fill the matrix of trajectory id
        lifetime = 1;  	% Initialization of the eddy lifetime
        % Fill the trajectory parameters whith the first eddy characteristics
        Cyclonic_Trajectories(numtraj,lifetime,1) = Amplitude(i,j);
        Cyclonic_Trajectories(numtraj,lifetime,2) = Radius(i,j);
        Cyclonic_Trajectories(numtraj,lifetime,3) = EKE(i,j);
        Cyclonic_Trajectories(numtraj,lifetime,4) = 0;
        Cyclonic_Trajectories(numtraj,lifetime,5) = Xcenter(i,j);
        Cyclonic_Trajectories(numtraj,lifetime,6) = Ycenter(i,j);
        % While eddy = 1, each eddy found takes part of the trajectory
        % A new eddy is defined as the first eddy of a trajectory
        eddy = 1;
        
        %  while eddy==1 % While an eddy is detected and it is not the last day of the data time coverage
        
        % Load contour and center position of detected eddy
        Xcon0 =double(CEs{indx,indy,1});
        Ycon0 =double(CEs{indx,indy,2});
        Xc = Xcenter(indx,indy);
        Yc = Ycenter(indx,indy);
        
        % search for the closest eddies on the next map
        [next_map_eddies,tot_dist,found_centroid_at] = cori_find_eddy_next_map(Xcentroid,Ycentroid,indx,indy,1,dist_lim,nc);
        
        if isempty(next_map_eddies)==0 % if a close eddy is found
            % testing the contours intersections
            [matching_contour,area_intersection] = cori_matching_contours(CEs,Area,indx,indy,1,next_map_eddies,tot_dist);
            if isempty(matching_contour)==0
                indd = matching_contour;
                lifetime = lifetime +1;
                Nestimate = Nestimate + 1;
                waitbar(Nestimate/(Nestimate_lim),h)
                
                if Nestimate>Nestimate_lim
                    break
                end
                disp(['Line 334; CE Nestimate:   ' int2str(Nestimate)  ' idx:  ' int2str(indx)])
                
                if length(matching_contour)==1 % only one contour
                elseif  length(matching_contour) >1 % more than 1 contour -> cost function
                    
                    cf_radius=((Radius(indd,indy+1)-Radius(indx,indy))/Radius(indx,indy)).^2;
                    cf_eke=((EKE(indd,indy+1)-EKE(indx,indy))/EKE(indx,indy)).^2;
                    cf_amp=((Amplitude(indd,indy+1)-Amplitude(indx,indy))/Amplitude(indx,indy)).^2;
                    
                    
                    cf =  cf_radius + cf_eke + cf_amp;
                    [~,ii] = min(cf); % Find the the minimum of the cost function
                    indd = indd(ii);  % Define the corresponding eddy as the same eddy at the previous time
                end
                
                % Fill the trajectory
                Cyclonic_Trajectories(numtraj,lifetime,1)=Amplitude(indd,indy+1);
                Cyclonic_Trajectories(numtraj,lifetime,2)=Radius(indd,indy+1);
                Cyclonic_Trajectories(numtraj,lifetime,3)=EKE(indd,indy+1);
                Cyclonic_Trajectories(numtraj,lifetime,4)=ac_distance(Yc,Xc,Ycenter(indd,indy+1),Xcenter(indd,indy+1));
                Cyclonic_Trajectories(numtraj,lifetime,5)=Xcenter(indd,indy+1);
                Cyclonic_Trajectories(numtraj,lifetime,6)=Ycenter(indd,indy+1);
                
                % If this eddy was not detected as taking part of a trajectory before
                if isnan(idtraj(indd,indy+1))==1
                    
                    % Store the ID of the trajectory associated with this eddy
                    idtraj(indd,indy+1) = numtraj;
                    
                end
                % Go to next map
                indy = indy+1;
                indx = indd;
                
                if indy ==  nc % The last map is reached
                    %disp('Trajectory might continue longer than the period')
                    eddy = 0;
                end
                
            else % No Eddy -> trajectory stops -> break while
                eddy = 0;
                %                     disp('No eddy intersection')
            end
            
        else % NO center near enough -> trajectory stops
            eddy = 0;
            %                 disp('Too far')
        end
    end
    
end
close(h);
disp(numtraj)

%% Statistics

x = Cyclonic_Trajectories(1:numtraj,:,5)'; x = x(:);
y = Cyclonic_Trajectories(1:numtraj,:,6)'; y = y(:);
delt_dist_center = ac_distance(y(2:end),x(2:end),y(1:end-1),x(1:end-1));
delt_Ampl = diff(Cyclonic_Trajectories(1:numtraj,:,1),1,2);
delt_Rad = diff(Cyclonic_Trajectories(1:numtraj,:,2),1,2);
delt_EKE = diff(Cyclonic_Trajectories(1:numtraj,:,3),1,2);

mean_delt_Ampl = nanmean(delt_Ampl(:));
std_delt_Ampl = nanstd(delt_Ampl(:));

mean_delt_Rad = nanmean(delt_Rad(:));
std_delt_Rad = nanstd(delt_Rad(:));

mean_delt_EKE = nanmean(delt_EKE(:));
std_delt_EKE = nanstd(delt_EKE(:));

mean_delt_dist = nanmean(delt_dist_center(:));
std_delt_dist = nanstd(delt_dist_center(:));

% Save data in previous matrices
eval([' save -append ' filename ' Metadata mean* std* numtraj'])

%% Output values
clc
display('OBTAINED STATISTICS')
display(' ')
display('ANTICYCLONES: ')
filename = [input_dir '/AE_Properties.mat'];
load(filename)
display(['     Delta Amplitudes :    ' num2str([100*mean_delt_Ampl 100*std_delt_Ampl])])
display(['     Delta Radii :    ' num2str([mean_delt_Rad std_delt_Rad])])
display(['     Delta EKE :    ' num2str([10000*mean_delt_EKE 10000*std_delt_EKE])])
display(['     Delta Distances :    ' num2str([mean_delt_dist std_delt_dist])])

display(' ')
display('CYCLONES: ')
filename = [input_dir '/CE_Properties.mat'];
load(filename)
display(['     Delta Amplitudes :    ' num2str([100*mean_delt_Ampl 100*std_delt_Ampl])])
display(['     Delta Radii :    ' num2str([mean_delt_Rad std_delt_Rad])])
display(['     Delta EKE :    ' num2str([10000*mean_delt_EKE 10000*std_delt_EKE])])
display(['     Delta Distances :    ' num2str([mean_delt_dist std_delt_dist])])

end