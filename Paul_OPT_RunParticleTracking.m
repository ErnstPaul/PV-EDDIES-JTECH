%% Paul_OPT_RunParticleTracking.m (version 4.1)
%Author: Yves Morel
%Date Created: 02/2020
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 7/22/2022 - Continued work, finished RETRO_PT loops
%PE 6/23/2022 - Generalized and revamped portability
%--------------------------------------
%Purpose: Runs the particle tracking suite in regards to the particles flowing into an eddy
%Inputs: Processed NEMO files for a date and time integration window, lats/lons in question
%Outputs: Particle tracking outputs as per Assene et al., 2020
%--------------------------------------
%% Inputs
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
partPath = [basepath '/PARTICLE_TRACKING/'];
addpath(partPath)

%String to use for the run
titlestring = 'NEM_RSW_CHAGOS_BACKWARD';

%Temporal variables
date2start  =  '20190529'; %YYYYMMDD day to start on in single quotes
days2integrate  =  100; %number of days to integrate forward or backward
forwardIntegrate  =  false; %integrate forwards in time (dispersion) - true or false
backwardIntegrate  =  true; %integrate backwards in time (accumulation) - opposite of above

%Spatial variables
lon0 = 74.;%central position in °
lat0 = -8;%central position in °

%Reference profile variables
refLon = 72;
refLat = 0;
%refyear is the one used above
refmonth = 5;
refday = 1;

%For the following, check the eddy in question against density (use Paul_OPT_3DColumn.m)
zed0  =  736.6; %central position in m (positive)
dzedvar0  =  100; %maximum vertical distribution range (in m)

%Radius of the eddy in km
dradius  =  112; %maximum horizontal distribution range (in km)

%Generally these are set but you can tune these
Npart  =  500; %number of particles to consider (try with 1 at first ... but for good diag use 500)
vortlim  =  1e-6; %limit of vorticity for particles to be selected!!!ATTENTION TO SIGN!!! USE POS FOR CE, NEG FOR AE

%% CALCULATING DERIVED PARAMETERS & INITIALIZING
output_dir = [partPath 'RESULTS/' titlestring];
[status,message,messageid] = mkdir(output_dir);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Area of calculation (ZOOM)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Eddy Radius to degrees, roughly
dradius  =  dradius/111.32;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of time parameters to pre-initialize matrices for speed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
minYear = 2015; %NEMO model v3.3.1
% Time of integration of trajectory
DT = 1; %frequency of outputs from NEMO model in days
DTSEC = DT*86400; % number of seconds per day
delt = 3600; %timestep in seconds for trajectory calculation (this is assumed to be an hour)
Ntot = floor(DT*86400/delt); %total hours per day
delt = DT*86400./Ntot; %seconds per hours
Nday = floor(86400/delt); %number of hours per day (if DT  =  1, same as Ntot)
if forwardIntegrate
    %Calculation of date params
    day2use = str2double(date2start(7:8));
    month2use = str2double(date2start(5:6));
    year2use = str2double(date2start(1:4));
    Tinit = datetime(year2use,month2use,day2use);
    Tend = Tinit+days2integrate;
    length2use = length(Tinit:Tend)*length(1:Ntot)+1;
elseif backwardIntegrate
    %Calculation of date params
    day2use = str2double(date2start(7:8));
    month2use = str2double(date2start(5:6));
    year2use = str2double(date2start(1:4));
    Tinit = datetime(year2use,month2use,day2use);
    Tend = Tinit-days2integrate;
    length2use = length(Tend:Tinit)*length(1:Ntot)+1;
end

% Initializing final matrices
lonalltraj = zeros(length2use,Npart);
latalltraj= zeros(length2use,Npart);
zedalltraj= zeros(length2use,Npart);
ualltraj= zeros(length2use,Npart);
valltraj= zeros(length2use,Npart);
walltraj= zeros(length2use,Npart);
temperaturealltraj= zeros(length2use,Npart);
salinityalltraj= zeros(length2use,Npart);
densityalltraj= zeros(length2use,Npart);
PValltraj= zeros(length2use,Npart);
vortalltraj= zeros(length2use,Npart);
Rialltraj= zeros(length2use,Npart);

%% Acquire Reference Profile
refyear = year2use-minYear;
[kMaxRhoRef, rhoRef, zRef] = refProfPV(refLon, refLat, refyear, refmonth, refday);
refProf = {zRef, rhoRef, kMaxRhoRef};

%% TRAJEC STEP 1: VORTICITY CALCULATION AND DETERMINATION
%Time step on each particle
for ipart = 1:Npart

    tic %timing

    %Initializing these matrices for this particle
    lontraj = nan(length2use,1);
    lattraj = nan(length2use,1);
    zedtraj = nan(length2use,1);
    utraj = nan(length2use,1);
    vtraj = nan(length2use,1);
    wtraj = nan(length2use,1);
    temperaturetraj = nan(length2use,1);
    salinitytraj = nan(length2use,1);
    densitytraj = nan(length2use,1);
    PVtraj = nan(length2use,1);
    vorttraj = nan(length2use,1);
    Ritraj = nan(length2use,1);

    % Make sure the selected particle reaches a minimum vorticity value
    % NOTE : the sign matters
    choosepart = 0; %no particle chosen
    while choosepart == 0

        %Determine particle starting position
        theta = 2*pi*rand; %random distribution around a point
        rdist = dradius*rand; %calculate away from radius
        dzedvar = dzedvar0*(2.*rand-1.); %calculate away from radius vertically
        lontraj(1) = lon0+rdist*cos(theta);%in °
        lattraj(1) = lat0+rdist*sin(theta);%in °
        zedtraj(1) = zed0+dzedvar;

        %Determine particle initial vorticity
        vorttraj(1) = Paul_OPT_Step1_PT(lontraj(1),lattraj(1),zedtraj(1),date2start);

        %Determine if particle vorticity meets minimum value
        if vortlim>0
            if vorttraj(1)<vortlim
                choosepart = 0; % do not keep particle, vorticity too low
            else
                choosepart = 1; % keep particle
            end
        else
            if vorttraj(1)>vortlim
                choosepart = 0;
            else
                choosepart = 1;
            end %end vorticity determination (1)
        end %end vorticity determination (1)
    end %end choosing particle loop

    disp(["Particle " + num2str(ipart) + " found in " +  toc +  " seconds. Proceeding to tracking."])

    %% STEPS 2 & 3: PARTICLE CHOSEN, GO FORWARD OR BACKWARD IN TIME TO TRACK THE PARTICLE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now loop over time and calculate trajectory
    % and tracer evolution for particle
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    tic
    %From start to finish
    if forwardIntegrate
        % Throw into trajectory calculation function
        [lontraj,lattraj,zedtraj,utraj,vtraj,wtraj,temperaturetraj,salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj,timetraj]...
            = Paul_OPT_Step2FORWARD_PT(lontraj,lattraj,zedtraj,utraj,vtraj,wtraj,temperaturetraj,...
            salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj,ipart,Tinit,Tend,refProf);
    %From finish to start
    elseif backwardIntegrate
        % Throw into trajectory calculation function
        [lontraj,lattraj,zedtraj,utraj,vtraj,wtraj,temperaturetraj,salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj,timetraj]...
            = Paul_OPT_Step2RETRO_PT(lontraj,lattraj,zedtraj,utraj,vtraj,wtraj,temperaturetraj,...
            salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj,ipart,Tinit,Tend,refProf);
    end%ENDINTEGRATION
    
    disp(["Particle " + num2str(ipart) + " trajectory finished in " +  toc +  " seconds. Proceeding to save."])

    % Take care of initial values (not calculated except for trajectory)
    utraj(1) = utraj(2);
    vtraj(1) = vtraj(2);
    wtraj(1) = wtraj(2);
    temperaturetraj(1) = temperaturetraj(2);
    salinitytraj(1) = salinitytraj(2);
    densitytraj(1) = densitytraj(2);
    PVtraj(1) = PVtraj(2);
    vorttraj(1) = vorttraj(2);
    Ritraj(1) = Ritraj(2);

    % Save fields back into master matrices
    lonalltraj(:,ipart) = lontraj;
    latalltraj(:,ipart) = lattraj;
    zedalltraj(:,ipart) = zedtraj;
    ualltraj(:,ipart) = utraj;
    valltraj(:,ipart) = vtraj;
    walltraj(:,ipart) = wtraj;
    temperaturealltraj(:,ipart) = temperaturetraj;
    salinityalltraj(:,ipart) = salinitytraj;
    densityalltraj(:,ipart) = densitytraj;
    PValltraj(:,ipart) = PVtraj;
    vortalltraj(:,ipart) = vorttraj;
    Rialltraj(:,ipart) = Ritraj;

    % save each particle traj in file
    file4save = [partPath 'RESULTS/' titlestring '/PARTICLE' num2str(ipart)];
    save(file4save,'lontraj','lattraj','zedtraj','utraj','vtraj','wtraj','temperaturetraj','salinitytraj','densitytraj','PVtraj','vorttraj','Ritraj','Tinit','Tend')
    disp(["Particle " + num2str(ipart) + " saved in " +  toc +  " seconds. Proceeding to next particle."])

end %ENDPARTICLELOOP

%% Save master matrices into complete file
for ipart = 1:Npart
    tic
    % Save fields back into master matrices
    file4load = [partPath 'RESULTS/' titlestring '/PARTICLE' num2str(ipart)];
    load(file4load)
    lonalltraj(:,ipart) = lontraj;
    latalltraj(:,ipart) = lattraj;
    zedalltraj(:,ipart) = zedtraj;
    ualltraj(:,ipart) = utraj;
    valltraj(:,ipart) = vtraj;
    walltraj(:,ipart) = wtraj;
    temperaturealltraj(:,ipart) = temperaturetraj;
    salinityalltraj(:,ipart) = salinitytraj;
    densityalltraj(:,ipart) = densitytraj;
    PValltraj(:,ipart) = PVtraj;
    vortalltraj(:,ipart) = vorttraj;
    Rialltraj(:,ipart) = Ritraj;
    disp("Particle " + num2str(ipart) + " Placed Into All Matrix in " + toc + " Seconds.")
end
file4save = [partPath 'RESULTS/' titlestring '/ALL'];
save(file4save,'lonalltraj','latalltraj','zedalltraj','ualltraj','valltraj','walltraj','temperaturealltraj','salinityalltraj','densityalltraj','PValltraj','vortalltraj','Rialltraj','Tinit','Tend','timetraj')
disp(["All particles saved."])


