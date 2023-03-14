%% Paul_OPT_Step2FORWARD_PT.m (version 2.2)
%Author: Yves Morel
%Date Created: 02/2020
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/24/2022 - Generalized and revamped portability
%--------------------------------------
%Purpose: Runs trajectory calculation for a given particle
%Inputs: Inputs from Paul_OPT_RunParticleTracking
%Outputs: Previous position for a given particle
%--------------------------------------
function [lontraj,lattraj,zedtraj,utraj,vtraj,wtraj,temperaturetraj,salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj,timetraj]...
    = Paul_OPT_Step2FORWARD_PT(lontraj,lattraj,zedtraj,utraj,vtraj,wtraj,temperaturetraj,salinitytraj,densitytraj,PVtraj,vorttraj,Ritraj,ipart,Tinit,Tend,refProf)

%% Inputs
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculation of Static Time parameters
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Time of integration of trajectory
DT = 1; %frequency of outputs from NEMO model in days
DTSEC = DT*86400; % number of seconds per day
delt = 3600; %timestep in seconds for trajectory calculation (this is assumed to be an hour)
Ntot = floor(DT*86400/delt); %total hours per day
delt = DT*86400./Ntot; %seconds per hours

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Static Misc. Parameters, do not change
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference Pressure for potential density
Rearth = 6378137.; % Earth radius used for lon/lat calculation
ntraj = 1;
%constants for derivation
np2=4; np1=3; n=2; nm1=1; nder=1; np1der=2;

%% Prior-to-trajectory determinations
%find cell indices within which particle is postionned before step
icel = nearestpoint(lontraj(1),lonsNEM);
jcel = nearestpoint(lattraj(1),latsNEM);
kcel = nearestpoint(zedtraj(1),depthsNEM);

%Calculate these maxima for derivation later
maxlonsNEM = max(lonsNEM);
maxlatsNEM = max(latsNEM);
maxdepthsNEM = max(depthsNEM);

%% Loop for trajectory calculation (days)
length2use = length(Tinit:Tend)*length(1:Ntot);
timetraj = zeros(1,length2use);
timetraj(1)=0.5/Ntot;
for time = Tinit:Tend % loop for trajectory on days
    tic
    disp("Performing forwards calculations for particle " + num2str(ipart) + " on " + datestr(time))

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %% Now loop with time step of hours between outputs
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for dt=1:Ntot
        % First load velocity fields at t=time/-1/+1/+2 and calculate time
        % derivatives for interpolation
        %find indices within which particle is postioned before step
        icel = nearestpoint(lontraj(ntraj),lonsNEM);
        jcel = nearestpoint(lattraj(ntraj),latsNEM);
        kcel = nearestpoint(zedtraj(ntraj),depthsNEM);

        %Handling extreme edge cases at grid edge (should not happen most of the time)
        if icel == 1
            icel = icel+1;
        end
        if jcel == 1
            icel = icel+1;
        end
        if kcel == 1
            icel = icel+1;
        end
        if isnan(icel)
            return;
        end
        if isnan(jcel)
            return;
        end
        if isnan(kcel)
            return;
        end

        %grab the file to use
        date2start = [time-1, time, time+1, time+2];
        day2use = day(date2start);
        month2use = month(date2start);
        year2use = year(date2start);
        minYear = 2015; %NEMO model v3.3.1
        year2use = year2use - minYear;
        file2use = input_final_NEM{day2use(1),month2use(1),year2use(1)}; %file determination
        %Load U and V
        u3D = ncread(file2use,'uo',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
        v3D = ncread(file2use,'vo',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
        [x, y, z] = size(u3D);
       %create master timewise matrix
        u = zeros(x,y,z,4);
        v = zeros(x,y,z,4);
        u(:,:,:,1) = reshape(inpaint_nans(u3D),[2,2,2]);
        v(:,:,:,1) = reshape(inpaint_nans(v3D),[2,2,2]);

        %Step on dates for our full read of the situation
        for dateStep = 2:4
            file2use = input_final_NEM{day2use(dateStep),month2use(dateStep),year2use(dateStep)}; %file determination
            %Load U and V
            u3D = ncread(file2use,'uo',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
            v3D = ncread(file2use,'vo',[icel-1,jcel-1,kcel-1,1],[2,2,2,1]);
            u(:,:,:,dateStep) = reshape(inpaint_nans(u3D),[2,2,2]);
            v(:,:,:,dateStep) = reshape(inpaint_nans(v3D),[2,2,2]);
        end%endimmediatedateloop

        %Calculate W from continuity equation because we don't have its raw output here
        w = zeros(size(u));
        for i = 1:2
            for j = 1:2
                for t = 1:4
                    w(i,j,:,t) = (gradient(squeeze(u(i,j,:,t))) + gradient(squeeze(v(i,j,:,t)))).*-1;
                end
            end
        end

        %DERIVATIVE STEPS - U
        derivu(:,:,:,nder)=(u(:,:,:,np1)-u(:,:,:,nm1))/(2*DTSEC);
        derivu(:,:,:,np1der)=(u(:,:,:,np2)-u(:,:,:,n))/(2*DTSEC);
        deltau1=(u(:,:,:,np1)-u(:,:,:,n))/(DTSEC);
        deltau2=(3*deltau1-2*derivu(:,:,:,nder)-derivu(:,:,:,np1der))/DTSEC;
        deltau3=(-2*deltau1+derivu(:,:,:,nder)+derivu(:,:,:,np1der))/(DTSEC^2);
        %DERIVATIVE STEPS - V
        derivv(:,:,:,nder)=(v(:,:,:,np1)-v(:,:,:,nm1))/(2*DTSEC);
        derivv(:,:,:,np1der)=(v(:,:,:,np2)-v(:,:,:,n))/(2*DTSEC);
        deltav1=(v(:,:,:,np1)-v(:,:,:,n))/(DTSEC);
        deltav2=(3*deltav1-2*derivv(:,:,:,nder)-derivv(:,:,:,np1der))/DTSEC;
        deltav3=(-2*deltav1+derivv(:,:,:,nder)+derivv(:,:,:,np1der))/(DTSEC^2);
        %DERIVATIVE STEPS - W
        derivw(:,:,:,nder)=(w(:,:,:,np1)-w(:,:,:,nm1))/(2*DTSEC);
        derivw(:,:,:,np1der)=(w(:,:,:,np2)-w(:,:,:,n))/(2*DTSEC);
        deltaw1=(w(:,:,:,np1)-w(:,:,:,n))/(DTSEC);
        deltaw2=(3*deltaw1-2*derivw(:,:,:,nder)-derivw(:,:,:,np1der))/DTSEC;
        deltaw3=(-2*deltaw1+derivw(:,:,:,nder)+derivw(:,:,:,np1der))/(DTSEC^2);
        %DERIVATIVE STEPS - TIME
        dtsec=(dt-1)*delt;
        u4D=u(:,:,:,n)+dtsec*derivu(:,:,:,nder)+(dtsec^2)*deltau2+(dtsec^3)*deltau3;
        v4D=v(:,:,:,n)+dtsec*derivv(:,:,:,nder)+(dtsec^2)*deltav2+(dtsec^3)*deltav3;
        w4D=w(:,:,:,n)+dtsec*derivw(:,:,:,nder)+(dtsec^2)*deltaw2+(dtsec^3)*deltaw3;

        % Now advance position
        ntraj=ntraj+1;
        timetraj(ntraj)=timetraj(ntraj-1)+1/Ntot;%in days, remember it is offset by 0.5/Ntot

        % Now interpolate fields spatially for each component lam, phi, z
        
        % U
        % Interpolate U field at particle position
        % cell for U is icel, icel+1 ; jcel, jcel+1 ; kcel, kcel+1
        Utraj=0.;
        ptottraj=0.;
        dlon=abs(lonsNEM(icel)-lonsNEM(icel+1));
        dlat=abs(latsNEM(jcel)-latsNEM(jcel+1));
        dzed=abs(depthsNEM(kcel)-depthsNEM(kcel+1));
        for i=1:2
            for j=1:2
                for k=1:2
                    ptraj=abs(1-abs(lonsNEM(icel+i-1)-lontraj(ntraj-1))/dlon);
                    ptraj=ptraj*abs(1-abs(latsNEM(jcel+j-1)-lattraj(ntraj-1))/dlat);
                    ptraj=ptraj*abs(1-abs(depthsNEM(kcel+k-1)-zedtraj(ntraj-1))/dzed);
                    Utraj=Utraj+ptraj*u4D(i,j,k);
                    ptottraj=ptottraj+ptraj;
                end
            end
        end
        Utraj=Utraj/ptottraj;
        utraj(ntraj)=Utraj;
        dlondt=(Utraj/(Rearth*cos(lattraj(ntraj-1)*pi/180)))*(180/pi);%convert velocity into dlon/dt "speed" (in degrees)
        % Now calculate new position at timetraj
        lontraj(ntraj)=lontraj(ntraj-1)-dlondt*delt;
        
        % V
        % Interpolate V field at particle position
        % cell for V is icel, icel+1 ; jcel, jcel+1 ; kcel, kcel+1
        Vtraj=0.;
        ptottraj=0.;
        dlon=abs(lonsNEM(icel)-lonsNEM(icel+1));
        dlat=abs(latsNEM(jcel)-latsNEM(jcel+1));
        dzed=abs(depthsNEM(kcel)-depthsNEM(kcel+1));
        for i=1:2
            for j=1:2
                for k=1:2
                    ptraj=abs(1-abs(lonsNEM(icel+i-1)-lontraj(ntraj-1))/dlon);
                    ptraj=ptraj*abs(1-abs(latsNEM(jcel+j-1)-lattraj(ntraj-1))/dlat);
                    ptraj=ptraj*abs(1-abs(depthsNEM(kcel+k-1)-zedtraj(ntraj-1))/dzed);
                    Vtraj=Vtraj+ptraj*v4D(i,j,k);
                    ptottraj=ptottraj+ptraj;
                end
            end
        end
        Vtraj=Vtraj/ptottraj;
        vtraj(ntraj)=Vtraj;
        dlatdt=(Vtraj/(Rearth))*(180/pi);%convert velocity into dlat/dt "speed" (in degrees)
        % Now calculate new position at timetraj
        lattraj(ntraj)=lattraj(ntraj-1)-dlatdt*delt;
        
        % W
        % Interpolate W field at particle position
        % cell for W is icel, icel+1 ; jcel, jcel+1 ; kcel, kcel+1
        Wtraj=0.;
        ptottraj=0.;
        dlon=abs(lonsNEM(icel)-lonsNEM(icel+1));
        dlat=abs(latsNEM(jcel)-latsNEM(jcel+1));
        dzed=abs(depthsNEM(kcel)-depthsNEM(kcel+1));
        for i=1:2
            for j=1:2
                for k=1:2
                    ptraj=abs(1-abs(lonsNEM(icel+i-1)-lontraj(ntraj-1))/dlon);
                    ptraj=ptraj*abs(1-abs(latsNEM(jcel+j-1)-lattraj(ntraj-1))/dlat);
                    ptraj=ptraj*abs(1-abs(depthsNEM(kcel+k-1)-zedtraj(ntraj-1))/dzed);
                    Wtraj=Wtraj+ptraj*w4D(i,j,k);
                    ptottraj=ptottraj+ptraj;
                end
            end
        end
        Wtraj=Wtraj/ptottraj;
        wtraj(ntraj)=Wtraj;
        dzeddt=Wtraj;%convert velocity into dZ/dt "speed"
        % Now calculate new position at timetraj
        zedtraj(ntraj)=zedtraj(ntraj-1)-dzeddt*delt;
        
        %Calculate different tracer values at new particle position
        [temperaturetraj(ntraj),salinitytraj(ntraj),densitytraj(ntraj),PVtraj(ntraj),vorttraj(ntraj),Ritraj(ntraj)]...
    = Paul_OPT_Step3FORWARD_TRACER(lontraj(ntraj),lattraj(ntraj),zedtraj(ntraj),icel,jcel,kcel,day2use,month2use,year2use,dt,refProf,u4D,v4D);

        disp(["Particle " + num2str(ipart) + " trajectory step " + num2str(ntraj) + " of " + num2str(length2use) +  " PV calculated in " +  toc +  " seconds. Proceeding to next time step."])
    end
end
