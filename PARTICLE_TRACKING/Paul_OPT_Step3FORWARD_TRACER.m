%% Paul_OPT_Step3FORWARD_TRACER.m (version 2.1)
%Author: Yves Morel
%Date Created: 02/2020
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 7/22/2022 - Generalized and revamped portability
%--------------------------------------
%Purpose: Runs calculation of parameters for a given particle
%Inputs: Inputs from Paul_OPT_Step2RETRO_PT
%Outputs: PV, Ri, etc. for a given particle
%--------------------------------------
function [temperaturetraj,salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj]...
    = Paul_OPT_Step3FORWARD_TRACER(lontraj,lattraj,zedtraj,icel,jcel,kcel,day2use,month2use,year2use,dt,refProf,u4D,v4D)

%% Inputs
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m

%unpacking reference profile variables
zref = refProf{1};
rhoref = refProf{2};
kmaxrhoref = refProf{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Static Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time of integration of trajectory
DT = 1; %frequency of outputs from NEMO model in days
DTSEC = DT*86400; % number of seconds per day
delt = 3600; %timestep in seconds for trajectory calculation (this is assumed to be an hour)
Ntot = floor(DT*86400/delt); %total hours per day
delt = DT*86400./Ntot; %seconds per hours
dtsec=(Ntot+1-dt)*delt-0.5/Ntot;%at new position and new time, offset = - 0.5/Ntot

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static Misc. Parameters, do not change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Pressure for potential density
Pref = 0.; %surface
rho0 = 1024.; % reference density (used for Ri)
g = 9.81; % gravity
%constants for derivation
np2=4; np1=3; n=2; nm1=1; nder=1; np1der=2;

%% Load variables
%Load S and T
file2use = input_final_NEM{day2use(1),month2use(1),year2use(1)};
s3D = ncread(file2use,'so',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
t3D = ncread(file2use,'thetao',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
[x, y, z] = size(s3D);
%create master timewise matrix
temp4DpreCT = zeros(x,y,z,4);
sal4D = zeros(x,y,z,4);
temp4DpreCT(:,:,:,1) = reshape(inpaint_nans(t3D),[2,2,2]);
sal4D(:,:,:,1) = reshape(inpaint_nans(s3D),[2,2,2]);

%Step on dates for our full read of the situation
for dateStep = 2:4
    file2use = input_final_NEM{day2use(dateStep),month2use(dateStep),year2use(dateStep)}; %file determination
    %Load U and V
    t3D = ncread(file2use,'thetao',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
    s3D = ncread(file2use,'so',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
    temp4DpreCT(:,:,:,dateStep) = reshape(inpaint_nans(t3D),[2,2,2]);
    sal4D(:,:,:,dateStep) = reshape(inpaint_nans(s3D),[2,2,2]);
end%endimmediatedateloop

%Convert T to CT
pNEMgrid = gsw_p_from_z(depthsNEM(kcel-1:kcel)*-1,latsNEM(jcel-1:jcel));
pNEMgrid = [pNEMgrid; pNEMgrid];
temp4D = zeros(x,y,z,4);
for latLoop = 1:2
    for timeLoop = 1:4
        temp4D(:,latLoop,:,timeLoop) = gsw_CT_from_t(squeeze(sal4D(:,latLoop,:,timeLoop)),...
            squeeze(temp4DpreCT(:,latLoop,:,timeLoop)),pNEMgrid);
    end
end

% TEMPERATURE
%First interpolate temporally
derivt(:,:,:,nder)=(temp4D(:,:,:,np1)-temp4D(:,:,:,nm1))/(2*1*86400);
derivt(:,:,:,np1der)=(temp4D(:,:,:,np2)-temp4D(:,:,:,n))/(2*1*86400);
deltat1=(temp4D(:,:,:,np1)-temp4D(:,:,:,n))/(1*1*86400);
deltat2=(3*deltat1-2*derivt(:,:,:,nder)-derivt(:,:,:,np1der))/DTSEC;
deltat3=(-2*deltat1+derivt(:,:,:,nder)+derivt(:,:,:,np1der))/(DTSEC^2);
% interpolation
temp2=temp4D(:,:,:,n)+dtsec*derivt(:,:,:,nder)+(dtsec^2)*deltat2+(dtsec^3)*deltat3;
temperature=temp2;
% Interpolate at particle position
% cell for tracer is icel, icel+1 ; jcel, jcel+1 ; kcel, kcel+1
Ttraj=0.;
ptottraj=0.;
dlon=abs(lonsNEM(icel)-lonsNEM(icel+1));
dlat=abs(latsNEM(jcel)-latsNEM(jcel+1));
dzed=abs(depthsNEM(kcel)-depthsNEM(kcel+1));
for i=1:2
    for j=1:2
        for k=1:2
            ptraj=abs(1-abs(lonsNEM(icel+i-1)-lontraj)/dlon);
            ptraj=ptraj*abs(1-abs(latsNEM(jcel+j-1)-lattraj)/dlat);
            ptraj=ptraj*abs(1-abs(depthsNEM(kcel+k-1)-zedtraj)/dzed);
            Ttraj=Ttraj+ptraj*temp2(i,j,k);
            ptottraj=ptottraj+ptraj;
        end
    end
end
temperaturetraj=Ttraj/ptottraj;

% SALINITY
derivt(:,:,:,nder)=(sal4D(:,:,:,np1)-sal4D(:,:,:,nm1))/(2*1*86400);
derivt(:,:,:,np1der)=(sal4D(:,:,:,np2)-sal4D(:,:,:,n))/(2*1*86400);
deltat1=(sal4D(:,:,:,np1)-sal4D(:,:,:,n))/(1*1*86400);
deltat2=(3*deltat1-2*derivt(:,:,:,nder)-derivt(:,:,:,np1der))/DTSEC;
deltat3=(-2*deltat1+derivt(:,:,:,nder)+derivt(:,:,:,np1der))/(DTSEC^2);
% interpolation
temp2=sal4D(:,:,:,n)+dtsec*derivt(:,:,:,nder)+(dtsec^2)*deltat2+(dtsec^3)*deltat3;
salinity=temp2;
% Interpolate at particle position
% cell for tracer is icel, icel+1 ; jcel, jcel+1 ; kcel, kcel+1
Ttraj=0.;
ptottraj=0.;
dlon=abs(lonsNEM(icel)-lonsNEM(icel+1));
dlat=abs(latsNEM(jcel)-latsNEM(jcel+1));
dzed=abs(depthsNEM(kcel)-depthsNEM(kcel+1));
for i=1:2
    for j=1:2
        for k=1:2
            ptraj=abs(1-abs(lonsNEM(icel+i-1)-lontraj)/dlon);
            ptraj=ptraj*abs(1-abs(latsNEM(jcel+j-1)-lattraj)/dlat);
            ptraj=ptraj*abs(1-abs(depthsNEM(kcel+k-1)-zedtraj)/dzed);
            Ttraj=Ttraj+ptraj*temp2(i,j,k);
            ptottraj=ptottraj+ptraj;
        end
    end
end
salinitytraj=Ttraj/ptottraj;

% DENSITY
densitytraj = gsw_rho(salinitytraj,temperaturetraj,Pref); %potential density
density=gsw_rho(salinity(:,:,:),temperature(:,:,:),Pref);
Grho=nan(size(density));
for i=1:2
    for j=1:2
        for k=1:2
            Grho(i,j,k)=rho2z4pv_new(kmaxrhoref,rhoref,zref,squeeze(density(i,j,k)));%rescaled PV
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Calculation of potential vorticity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%Only one reading here (compare to FilePreparation)
i=1; j=1; k=1;

%initialize coriolis parameter
fcorioloc = zeros(2,2);
for latLoop = -1:0
    % Calculate coriolis parameter across this area
    fcorioloc(:,latLoop+2) = gsw_f(latsNEM(kcel+latLoop,1));
end

%DX, DY, and DZ diffs for below
dZDiff(1:2) = depthsNEM(kcel)-depthsNEM(kcel-1);
dXDiff = ones(2,2)*(111.32/12)*1000; %% TIMES 1000 = m
for lonCalib = -1:0
    lonCalibrator = cosd(latsNEM(jcel+lonCalib));
    dXDiff(:,lonCalib+2) = dXDiff(:,lonCalib+2) .* lonCalibrator;
end
dYDiff = ones(2,2)*(110.573/12)*1000; %%TIMES 1000 = m

% Grabbing U and V from previously defined values in Step2
udef=zeros(size(density));
vdef=udef;
if lontraj>=lonsNEM(icel)
    udef(i,:,:)=u4D(i,:,:);
else
    udef(i,:,:)=u4D(i+1,:,:);
end
if lattraj>=latsNEM(jcel)
    vdef(:,j,:)=v4D(:,j,:);
else
    vdef(:,j,:)=v4D(:,j+1,:);
end

% East and West sides PV_x
zetax1 = -(vdef(i,j,k+1)-vdef(i,j,k))/dZDiff(k+1);
rhox1 = 0.25*(Grho(i,j,k)+Grho(i,j+1,k)+...
    Grho(i,j,k+1)+Grho(i,j+1,k+1));
zetax2 = -(vdef(i+1,j,k+1)-vdef(i+1,j,k))/dZDiff(k+1);
rhox2 = 0.25*(Grho(i+1,j,k)+Grho(i+1,j+1,k)+...
    Grho(i+1,j,k+1)+Grho(i+1,j+1,k+1));
pvx = (zetax2*rhox2-zetax1*rhox1)/dXDiff(i+1,j);
% South and North sides PV_y
zetay1 = (udef(i,j,k+1)-udef(i,j,k))/dZDiff(k+1);
rhoy1 = 0.25*(Grho(i,j,k)+Grho(i+1,j,k)+...
    Grho(i,j,k+1)+Grho(i+1,j,k+1));
zetay2 = (udef(i,j+1,k+1)-udef(i,j+1,k))/dZDiff(k+1);
rhoy2 = 0.25*(Grho(i,j+1,k)+Grho(i+1,j+1,k)+...
    Grho(i,j+1,k+1)+Grho(i+1,j+1,k+1));

pvy = (zetay2*rhoy2-zetay1*rhoy1)/dYDiff(i,j+1);
% Up and Down sides PV_z
zetaz1=(vdef(i+1,j,k)-vdef(i,j,k))/dXDiff(i+1,j);
zetaz1=zetaz1-(udef(i,j+1,k)-udef(i,j,k))/dYDiff(i,j+1);
fc=0.25*(fcorioloc(i,j)+fcorioloc(i+1,j)+...
    fcorioloc(i,j+1)+fcorioloc(i+1,j+1));
rhoz1= 0.25*( Grho(i,j,k)+Grho(i+1,j,k)+...
    Grho(i,j+1,k)+Grho(i+1,j+1,k));
zetaz2=(vdef(i+1,j,k+1)-vdef(i,j,k+1))/dXDiff(i+1,j);
zetaz2=zetaz2-(udef(i,j+1,k+1)-udef(i,j,k+1))/dYDiff(i,j+1);
rhoz2= 0.25*( Grho(i,j,k+1)+Grho(i+1,j,k+1)+...
    Grho(i,j+1,k+1)+Grho(i+1,j+1,k+1));
pvz = ((zetaz2+fc)*rhoz2-(zetaz1+fc)*rhoz1)/dZDiff(k+1);

% Total potential vorticity
PVtraj = pvx+pvy+pvz;

% Vorticity
vorttraj=0.5*(zetaz1+zetaz2);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% calculation of Richardson number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dzUpv=0.5*(udef(i,j,k+1)+udef(i,j+1,k+1)-udef(i,j,k)-udef(i,j+1,k))/dZDiff(k+1);
dzVpv=0.5*(vdef(i,j,k+1)+vdef(i+1,j,k+1)-vdef(i,j,k)-vdef(i+1,j,k))/dZDiff(k+1);
dzrhopv=0.25*( density(i,j,k)+density(i+1,j,k)+density(i,j+1,k)+density(i+1,j+1,k)-...
    density(i,j,k+1)-density(i+1,j,k+1)-density(i,j+1,k+1)-density(i+1,j+1,k+1) )/dZDiff(k+1);
Ritraj=-(g*dzrhopv/rho0)/(dzUpv^2+dzVpv^2);
end %ENDFUNCTION