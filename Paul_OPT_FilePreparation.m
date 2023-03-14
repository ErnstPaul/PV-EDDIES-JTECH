%% Paul_OPT_FilePreparation.m (version 3.2)
%Author: Paul Ernst
%Date Created: 5/11/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 9/23/2022 - Cleaned for final release
%PE 6/2/2022 - Changed for OPT
%PE 5/24/2022 - Added Yves Morel fixes (Dirac surface PV adjustments)
%PE 5/11/2022 - Created
%--------------------------------------
%Purpose: Grabs all NEMO data for a given area and saves it with the following variables:
%   Isopycnally Averaged PV at Surface
%   Isopycnally Averaged VORT at Surface
%   Surface: U, V, S, T, RHO, O-W Parameter, SSH, Vorticity, LNAM, EKE
%   X and Y
%MATLAB DateNum
%Inputs: NEMO Files, path information, time/spatial specifications (see MethodOptimization)
%Outputs: One output file for each day in the time period containing above variables
%--------------------------------------
function [output_dir] = Paul_OPT_FilePreparation(latlonbounds, yearmetadata, yearsMAT, monthsMAT, daysMAT, depthlevel, rhoBounds, outputDir,refLon,refLat)
%% Initialize variables and paths.
% Loading defaults and functions.
close all; clc
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
%months as a list of strings, yearsS as a list of strings and yearsS
addpath(strcat(basepath, 'FUNCTIONS'));
% Setting bounds based on mainfunction latlonbounds.
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);
NLnem = nearestpoint(NL, latsNEM);   SLnem = nearestpoint(SL, latsNEM);
ELnem = nearestpoint(EL, lonsNEM);   WLnem = nearestpoint(WL, lonsNEM);
% Setting number of Z coords to jump (max of model or of choice - NEMO max is 50)
dZ = 50;

% Make and confirm existence of output directory.
output_dir_interim = strcat(basepath, ['/EXTRACTION/INDIVIDUAL_FILES/' outputDir]);
output_dir = strcat(output_dir_interim, [num2str(depthlevel) + "/"]);
[status,message,messageid] = mkdir(output_dir);

% Create configfile struct and save to output directory.
Directory_Of_NetCDF_Files = input_dir_NEM; Directory_Of_Extracted_Files = output_dir;
Northern_Limit = NL; Southern_Limit = SL; Eastern_Limit = EL; Western_Limit = WL;
Beginning_Date = yearmetadata(1); Final_Date = yearmetadata(2);
Date_Of_Extraction = datestr(now);
savedest = strcat(output_dir, '/ConfigFile_fulltime.mat');
save(savedest, 'Directory_Of_NetCDF_Files','Directory_Of_Extracted_Files',...
    'Northern_Limit', 'Southern_Limit', 'Eastern_Limit','Western_Limit',...
    'Date_Of_Extraction', 'Beginning_Date', 'Final_Date');

%Set up mean and standard deviation variables for OW and VORT
OWSTATS = zeros(1,2); %MEAN, STD
VORTSTATS = zeros(1,2);
statsCount = 0;

%% Loop through each year to grab all relevant data and slot it appropriately

% Create waitbar, should work as intended
h = waitbar(0,'Extraction of Variables from NEMOv3.1');

% Datewise triple loop, setting up anomaly loop here as well
count = 0;
steps = length(yearsMAT)*length(monthsMAT)*length(daysMAT);

%% /END RHOBOUNDS

file2use = input_final_NEM{1,1,1};
b = [];

for yearMark=1:length(yearsMAT)
    year2use = find(yearsMAT(yearMark)==yearsSNEM);
    for monthNum = 1:length(monthsMAT)
        month2use = monthsMAT(monthNum);
        for dayNum = 1:length(daysMAT)
            day2useNEM = daysMAT(dayNum);
            tic
            %ensure correct index found here
            year = yearsSNEM(year2use); month = monthsS(month2use);
            if day2useNEM < 10
                dayS = insertBefore(num2str(day2useNEM),1,"0");
            else
                dayS = num2str(day2useNEM);
            end
            datecat = [year + month + dayS];
            count = count + 1;
            waitbar(count/steps,h) % Update waitbar
            % What is our single file that we want to load right now?
            file2use = input_final_NEM{day2useNEM,month2use,year2use};

            % Make sure file exists, some files may not (months with less than 31 days)
            if (isempty(file2use))
                continue
            end
            % If this data already exists, no need to recalculate
            filename_out = [output_dir + "pv_uv_" + datecat + ".mat"];
            if isfile(filename_out)
                disp(filename_out + " already exists, skipping.")
                continue
            end

            disp("Extracting: " + datecat + " [" + file2use + "]");
            % Load relevant variables from the file (T, S, U, V, Z)
            % Reads only the specified section of file and the depth levels above/below for PV
            % gradient calculations
            ptemps = ncread(file2use, 'thetao',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
            sals = ncread(file2use, 'so',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
            uNEM = ncread(file2use, 'uo',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
            vNEM = ncread(file2use, 'vo',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
            zosNEM = ncread(file2use, 'zos',...
                [WLnem,SLnem,1], [ELnem-WLnem, NLnem-SLnem, 1]);
            dNEM = transpose(double(ncread(file2use, 'depth', depthlevel, dZ)));

            maxDepthLev = length(dNEM);
            % Calculate difference in depths
            dZDiff(2:length(dNEM)) = diff(dNEM);
            dZDiff(1) = dNEM(1)*0.5;
            dXDiff = ones();

            %% Calculate PV
            % Create PV cell
            sizemat = size(ptemps);
            pvNEM = zeros(sizemat);

            %X and Y diffs
            dXDiff = ones(sizemat(1),sizemat(2))*(111.32/12)*1000; %% TIMES 1000?
            lonCalibrator = cosd(latsNEM(SLnem:NLnem));
            for i = 1:length(lonCalibrator)
                dXDiff(:,i) = dXDiff(:,i) .* lonCalibrator(i);
            end
            dYDiff = ones(sizemat(1),sizemat(2))*(110.573/12)*1000; %%TIMES 1000????

            % Convert depth to pressure across the entire latitude band
            pNEMgrid = zeros((NLnem-SLnem+1),length(dNEM));
            fNEM = zeros(sizemat);
            for latLoop = 1:(NLnem-SLnem+1)
                pNEMgrid(latLoop,:) = gsw_p_from_z(dNEM.*-1,latsNEM(SLnem+latLoop,1));

                % Calculate coriolis parameter across the entire grid
                fNEM(:,latLoop,:) = gsw_f(latsNEM(SLnem+latLoop,1));
            end

            % Uniformity across depth levels
            fNEMSingleLevel = fNEM(:,:,1);

            % Convert salinity to absolute salinity across the grid (!doubles proc time!)
            saNEM = zeros(sizemat);
            for la = 1:(NLnem-SLnem+1)
                for lo = 1:(ELnem-WLnem+1)
                    saNEM(lo,la,:) = gsw_SA_from_SP(transpose(squeeze(sals(lo,la,:))),...
                        pNEMgrid(la,:), lonsNEM(WLnem+lo,1), latsNEM(SLnem+la,1));
                end
            end

            % Comment above and uncomment this block if you use normal salinity
%             saNEM = sals;

            % Convert temperature to conservative temperature across the entire latitude band
            ctNEM = nan(sizemat);
            for latLoop = 1:(NLnem-SLnem+1)
                ctNEM(:,latLoop,:) = gsw_CT_from_t(squeeze(saNEM(:,latLoop,:)),squeeze(ptemps(:,latLoop,:)),pNEMgrid(latLoop,:));
            end

            % Calculate in-situ density across the entire longitude band
            rhoNEM = nan(sizemat);
            for lonLoop = 1:(ELnem-WLnem+1)
                rhoNEM(lonLoop,:,:) = gsw_rho(squeeze(saNEM(lonLoop,:,:)),squeeze(ctNEM(lonLoop,:,:)),0);
            end

            % Acquire new reference profile every monsoon -- adjust this if you want different ones
            % Feel free to specify different time intervals to grab different reference profiles
            if (day2useNEM == 1)
                if ((month2use == 1) || (month2use == 11))
                    refyear = year2use;
                    refmonth = 1;
                    refday = 1;
                    [kMaxRhoRef, rhoRef, zRef] = refProfPV(refLon, refLat, refyear, refmonth, refday);
                    disp("New winter monsoon: new reference profile of density acquired.")
                elseif (month2use == 7)
                    refyear = year2use;
                    refmonth = 7;
                    refday = 1;
                    [kMaxRhoRef, rhoRef, zRef] = refProfPV(refLon, refLat, refyear, refmonth, refday);
                    disp("New summer monsoon: new reference profile of density acquired.")
                end
            end

            % Rescaling function (calculates with respect to reference profile)
            Grho = nan(sizemat);
            for i = 1:sizemat(1)
                for j = 1:sizemat(2)
                    for k = 1:sizemat(3)
                        Grho(i,j,k) = rho2z4pv_new(kMaxRhoRef,rhoRef,zRef,squeeze(rhoNEM(i,j,k)));
                    end
                end
            end

            % Calculating PV from Density, U, and V
            PVNEM=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3));
            PVXNEM=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3)); 
            PVYNEM=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3)); 
            PVZNEM=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3));
            rhopv=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3));
            UPV=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3));
            VPV=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3));
            vort=nan(sizemat(1)-1,sizemat(2)-1,sizemat(3));
            for i=1:sizemat(1)-1
                for j=1:sizemat(2)-1
                    for k=depthlevel:(depthlevel+dZ-2)
                        % East and West sides PV_x
                        zetax1 = -(vNEM(i,j,k+1)-vNEM(i,j,k))/dZDiff(k+1);
                        rhox1 = 0.25*(Grho(i,j,k)+Grho(i,j+1,k)+...
                            Grho(i,j,k+1)+Grho(i,j+1,k+1));
                        zetax2 = -(vNEM(i+1,j,k+1)-vNEM(i+1,j,k))/dZDiff(k+1);
                        rhox2 = 0.25*(Grho(i+1,j,k)+Grho(i+1,j+1,k)+...
                            Grho(i+1,j,k+1)+Grho(i+1,j+1,k+1));
                        pvx = (zetax2*rhox2-zetax1*rhox1)/dXDiff(i+1,j);
                        % South and North sides PV_y
                        zetay1 = (uNEM(i,j,k+1)-uNEM(i,j,k))/dZDiff(k+1);
                        rhoy1 = 0.25*(Grho(i,j,k)+Grho(i+1,j,k)+...
                            Grho(i,j,k+1)+Grho(i+1,j,k+1));
                        zetay2 = (uNEM(i,j+1,k+1)-uNEM(i,j+1,k))/dZDiff(k+1);
                        rhoy2 = 0.25*(Grho(i,j+1,k)+Grho(i+1,j+1,k)+...
                            Grho(i,j+1,k+1)+Grho(i+1,j+1,k+1));

                        pvy = (zetay2*rhoy2-zetay1*rhoy1)/dYDiff(i,j+1);
                        % Up and Down sides PV_z
                        zetaz1=(vNEM(i+1,j,k)-vNEM(i,j,k))/dXDiff(i+1,j);
                        zetaz1=zetaz1-(uNEM(i,j+1,k)-uNEM(i,j,k))/dYDiff(i,j+1);
                        fc=0.25*(fNEMSingleLevel(i,j)+fNEMSingleLevel(i+1,j)+...
                            fNEMSingleLevel(i,j+1)+fNEMSingleLevel(i+1,j+1));
                        rhoz1= 0.25*( Grho(i,j,k)+Grho(i+1,j,k)+...
                            Grho(i,j+1,k)+Grho(i+1,j+1,k));
                        zetaz2=(vNEM(i+1,j,k+1)-vNEM(i,j,k+1))/dXDiff(i+1,j);
                        zetaz2=zetaz2-(uNEM(i,j+1,k+1)-uNEM(i,j,k+1))/dYDiff(i,j+1);
                        rhoz2= 0.25*( Grho(i,j,k+1)+Grho(i+1,j,k+1)+...
                            Grho(i,j+1,k+1)+Grho(i+1,j+1,k+1));

                        pvz = ((zetaz2+fc)*rhoz2-(zetaz1+fc)*rhoz1)/dZDiff(k+1) ;

                        if (k==1)
                            PVsurf(i,j)=(zetaz1+fc).*rhoz1; % Bretherton PV sheet
                        end

                        % Total potential vorticity
                        vort(i,j,k)=zetaz1;
                        rhoiso(i,j,k)= 0.25*(rhoNEM(i,j,k)+rhoNEM(i+1,j,k)+rhoNEM(i,j+1,k)+rhoNEM(i+1,j+1,k));
                        %
                        if (k==depthlevel+dZ-2)
                            vort(i,j,k+1)=zetaz2;
                            rhoiso(i,j,k+1)=0.25*(rhoNEM(i,j,k+1)+rhoNEM(i+1,j,k+1)+rhoNEM(i,j+1,k+1)+rhoNEM(i+1,j+1,k+1) );
                            latpv(j,1)=0.5*(latsNEM(j)+latsNEM(j+1));
                            lonpv(i,1)=0.5*(lonsNEM(i)+lonsNEM(i+1));
                            fcoriopv(i,j)=fc;
                        end
                        %
                        rhopv(i,j,k+1)= 0.125*( rhoNEM(i,j,k)+rhoNEM(i+1,j,k)+rhoNEM(i,j+1,k)+rhoNEM(i+1,j+1,k)+...
                            rhoNEM(i,j,k+1)+rhoNEM(i+1,j,k+1)+rhoNEM(i,j+1,k+1)+rhoNEM(i+1,j+1,k+1) );
                        %
                        UPV(i,j,k+1)=0.25*(uNEM(i,j,k+1)+uNEM(i,j,k)+uNEM(i,j+1,k+1)+uNEM(i,j+1,k));
                        VPV(i,j,k+1)=0.25*(vNEM(i,j,k+1)+vNEM(i,j,k)+vNEM(i+1,j,k+1)+vNEM(i+1,j,k));
                        %
                        PVNEM(i,j,k+1) = pvx+pvy+pvz; %pvz where coriolis operates
                        PVXNEM(i,j,k+1) = pvx; PVYNEM(i,j,k+1) = pvy; PVZNEM(i,j,k+1) = pvz;
                    end
                end
            end % PV calc loop
            PV = PVNEM;

            disp("PV Calculated.")
            %% LNAM CALCULATION - FROM AMEDA
            for k=depthlevel:(depthlevel+dZ-2)
                %----------------------------------------
                % Calculation of finite spatial element
                Y = latsNEM(SLnem:NLnem-1,1);
                X = lonsNEM(WLnem:ELnem-1,1);
                uNAM = uNEM(1:length(X),1:length(Y),k);
                vNAM = vNEM(1:length(X),1:length(Y),k);
                Xgrid = zeros(length(X),length(Y));
                Ygrid = zeros(length(X),length(Y));
                for p = 1:length(Y)
                    Xgrid(:,p) = X(:,1);
                end
                for p = 1:length(X)
                    Ygrid(p,:) = Y(:,1);
                end
                dx  = zeros(size(Xgrid));
                dy = zeros(size(Xgrid));
                dux = zeros(size(Xgrid));
                duy = zeros(size(Xgrid));
                dvx = zeros(size(Xgrid));
                dvy = zeros(size(Xgrid));

                %----------------------------------------
                %Spatial element in deg if grid_ll==1 or in km otherwise
                dx(2:end-1,2:end-1) = Xgrid(3:end,2:end-1) - Xgrid(1:end-2,2:end-1); %#ok<*COLND>
                dy(2:end-1,2:end-1) = Ygrid(2:end-1,3:end) - Ygrid(2:end-1,1:end-2);


                % define constants
                earth_radius = 6378.137; % km
                % kilometer (km) per degree of latitude
                R = earth_radius*pi/180; % 111.320m
                % Calcul finite spatial element in km
                dx = dx*R.*cosd(Ygrid);
                dy = dy*R;

                % in meters
                dx = dx*1000; % m
                dy = dy*1000; % m

                %Acquire mask and b variables
                mask = isnan(uNAM);
                if isempty(b)
                    b = AMEDA_fetchBfromRD("global_Rossby_Radius.mat",Xgrid,Ygrid,dx,mask);
                end

                %----------------------------------------
                % Compute speed element in m/s
                dux(2:end-1,2:end-1) = (uNAM(2:end-1,3:end) - uNAM(2:end-1,1:end-2));
                duy(2:end-1,2:end-1) = (uNAM(3:end,2:end-1) - uNAM(1:end-2,2:end-1));
                dvx(2:end-1,2:end-1) = (vNAM(2:end-1,3:end) - vNAM(2:end-1,1:end-2));
                dvy(2:end-1,2:end-1) = (vNAM(3:end,2:end-1) - vNAM(1:end-2,2:end-1));

                %----------------------------------------
                % Calculation of Okubo-Weiss criteria
                sn = (dux./dx) - (dvy./dy); % shear "cisaillement"
                ss = (dvx./dx) + (duy./dy); % strain "deformation"
                om = (dvx./dx) - (duy./dy); % vorticity "vorticité"

                okubo = sn.^2 + ss.^2 - om.^2; % in s-2

                %----------------------------------------
                % border is a parameter which prevents the constraints
                % to be applied to points too close to the domain boundaries
                % which would result in an index error
                borders = max(b(:)) + 1;

                %----------------------------------------
                % Calculation of LNAM criteria (Local Normalized Angular Momentum)
                % and LOW criteria (Local Averaged Okubo Weiss)
                %----------------------------------------

                % Initialisation
                if k == 1
                    L   = zeros([size(uNAM), maxDepthLev]);
                    LOW = nan([size(uNAM), maxDepthLev]);
                end

                %----------------------------------------
                % calculate LNAM and LOW in all domain pixels
                for i=borders:length(vNAM(:,1))-borders+1
                    for ii=borders:length(vNAM(1,:))-borders+1

                        if ~isnan(vNAM(i,ii))

                            % calculate LOW
                            OW = okubo(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii));
                            LOW(i,ii,k) = mean(OW(:));

                            % calculate LNAM
                            xlocal = Xgrid(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % Xgrid square sample of length b
                            ylocal = Ygrid(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % Ygrid square sample of length b

                            ulocal = uNAM(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % u square sample of length b
                            vlocal = vNAM(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % v square sample of length b

                            % Local Normalized Angular Momentum
                            % Use the midddle of the square as center coordinate

                            coordcentre = size(ulocal,1)-b(i,ii);

                            d_xcentre = (xlocal - xlocal(coordcentre,coordcentre));
                            d_ycentre = (ylocal - ylocal(coordcentre,coordcentre));

                            cross   = (d_xcentre.*vlocal) - (d_ycentre.*ulocal);
                            dot     = (ulocal.*d_xcentre) + (vlocal.*d_ycentre);
                            produit = sqrt(ulocal.^2 + vlocal.^2)...
                                .*sqrt(d_xcentre.^2 + d_ycentre.^2);
                            sumdp = sum(dot(:))+sum(produit(:));

                            if sumdp ~= 0
                                L(i,ii,k) = sum(cross(:)) / sumdp * sign(fNEM(i,ii,depthlevel));
                            else
                                L(i,ii,k) = 0;
                            end

                        end
                    end
                end

                thisL = L(:,:,k);

                thisL(isnan(thisL)) = 0;

                thisL(mask) = NaN;

                L(:,:,k) = thisL;
            end

            disp("LNAM computed.")
            %% EXTRACT ISOPYCNAL CHARACTERISTICS
            %set bounds of rho
            [rhoLen, ~] = size(rhoBounds);
            PVisopyc = nan(sizemat(1)-1,sizemat(2)-1,rhoLen);
            UPVisopyc=nan(sizemat(1)-1,sizemat(2)-1,rhoLen);
            VPVisopyc=nan(sizemat(1)-1,sizemat(2)-1,rhoLen);
            CTisopyc=nan(sizemat(1)-1,sizemat(2)-1,rhoLen);
            SAisopyc=nan(sizemat(1)-1,sizemat(2)-1,rhoLen);
            VORTisopyc=nan(sizemat(1)-1,sizemat(2)-1,rhoLen);
            for rhoStep = 1:rhoLen
                rhomin = rhoBounds(rhoStep,1);
                rhomax = rhoBounds(rhoStep,2);

                %loop thru i, j (X, Y)
                Zmax = nan(sizemat(1)-1,sizemat(2)-1);
                Zmin = nan(sizemat(1)-1,sizemat(2)-1);

                for i=1:sizemat(1)-1
                    for j=1:sizemat(2)-1
                        if min(isnan(rhoiso(i,j,:)))==1
                            PVisopyc(i,j,rhoStep)=NaN;
                            UPVisopyc(i,j,rhoStep)=NaN;
                            VPVisopyc(i,j,rhoStep)=NaN;
                            CTisopyc(i,j,rhoStep)=NaN;
                            SAisopyc(i,j,rhoStep)=NaN;
                            VORTisopyc(i,j,rhoStep)=NaN;
                            Zmin(i,j)=NaN;
                            Zmax(i,j)=NaN;
                        else
                            % Find z for rhomin
                            k0=min(find((rhoiso(i,j,:)-rhomin).^2==nanmin((rhoiso(i,j,:)-rhomin).^2)));
                            % SECURITY TO DEAL WITH SURFACE
                            % INVERSION OF DENSITY ...
                            if rhomin<rhoiso(i,j,1)
                                k0=1;
                            end
                            if rhoiso(i,j,k0)<rhomin
                                kmin=k0;
                                kp1=min(maxDepthLev,kmin+1);
                            else
                                if k0==1
                                    kmin=1;
                                    kp1=1;
                                else
                                    kmin=max(1,k0-1);
                                    kp1=min(maxDepthLev,kmin+1);
                                end
                            end
                            if kp1>kmin
                                zmin=dNEM(kmin)+(rhomin-rhoiso(i,j,kmin))*(dNEM(kp1)-dNEM(kmin))/...
                                    (rhoiso(i,j,kp1)-rhoiso(i,j,kmin));
                            else
                                zmin=dNEM(kmin);
                            end
                            %
                            % Find z for rhomax
                            k0=min(find((rhoiso(i,j,:)-rhomax).^2==nanmin((rhoiso(i,j,:)-rhomax).^2)));
                            if rhoiso(i,j,k0)<rhomax
                                if k0==maxDepthLev
                                    kmax=maxDepthLev;
                                    km1=kmax;
                                else
                                    kmax=min(maxDepthLev,k0+1);
                                    km1=max(1,kmax-1);
                                end
                            else
                                kmax=k0;
                                km1=max(1,kmax-1);
                            end
                            if km1<kmax
                                zmax=dNEM(kmax)+(rhomax-rhoiso(i,j,kmax))*(dNEM(km1)-dNEM(kmax))/...
                                    (rhoiso(i,j,km1)-rhoiso(i,j,kmax));
                            else
                                zmax=dNEM(kmax);
                            end
                            % layer thickness
                            Zmin(i,j)=zmin;
                            Zmax(i,j)=zmax;
                            % Integration of PV within layer between kmin/zmin and kmax/zmax
                            if kmin==kmax
                                PVisopyc(i,j,rhoStep)=NaN;
                                UPVisopyc(i,j,rhoStep)=NaN;
                                VPVisopyc(i,j,rhoStep)=NaN;
                                CTisopyc(i,j,rhoStep)=NaN;
                                SAisopyc(i,j,rhoStep)=NaN;
                                VORTisopyc(i,j,rhoStep)=NaN;
                            else
                                PVisopyc(i,j,rhoStep)=0.;
                                UPVisopyc(i,j,rhoStep)=0.;
                                VPVisopyc(i,j,rhoStep)=0.;
                                CTisopyc(i,j,rhoStep)=0.;
                                SAisopyc(i,j,rhoStep)=0.;
                                VORTisopyc(i,j,rhoStep)=0.;
                                for kloop=kmin+1:kmax
                                    PVisopyc(i,j,rhoStep)=PVisopyc(i,j,rhoStep)+PV(i,j,kloop)*dZDiff(kloop);
                                    UPVisopyc(i,j,rhoStep)=UPVisopyc(i,j,rhoStep)+UPV(i,j,kloop)*dZDiff(kloop);
                                    VPVisopyc(i,j,rhoStep)=VPVisopyc(i,j,rhoStep)+VPV(i,j,kloop)*dZDiff(kloop);
                                    CTisopyc(i,j,rhoStep)=CTisopyc(i,j,rhoStep)+ctNEM(i,j,kloop)*dZDiff(kloop);
                                    SAisopyc(i,j,rhoStep)=SAisopyc(i,j,rhoStep)+saNEM(i,j,kloop)*dZDiff(kloop);
                                    VORTisopyc(i,j,rhoStep)=VORTisopyc(i,j,rhoStep)+vort(i,j,kloop)*dZDiff(kloop);
                                end
                                % Include Dirac PV at surface for Mean
                                % isopycnal PV
                                % (Bretherton,1966/Schneider, 2003/Morel et al, 2019)
                                if kmin==1
                                    PVisopyc(i,j,rhoStep)=PVisopyc(i,j,rhoStep)+PVsurf(i,j);
                                end
                                % PV %
                                % remove part of the cell that is not within zmin - zmax
                                PVisopyc(i,j,rhoStep)=PVisopyc(i,j,rhoStep)-PV(i,j,kmin+1)*(zmin-dNEM(kmin));
                                PVisopyc(i,j,rhoStep)=PVisopyc(i,j,rhoStep)-PV(i,j,kmax)*(dNEM(kmax)-zmax);
                                % average over zmin - zmax
                                PVisopyc(i,j,rhoStep)=PVisopyc(i,j,rhoStep)/(zmax-zmin);
                                % UPV %
                                % remove part of the cell that is not within zmin - zmax
                                UPVisopyc(i,j,rhoStep)=UPVisopyc(i,j,rhoStep)-UPV(i,j,kmin+1)*(zmin-dNEM(kmin));
                                UPVisopyc(i,j,rhoStep)=UPVisopyc(i,j,rhoStep)-UPV(i,j,kmax)*(dNEM(kmax)-zmax);
                                % average over zmin - zmax
                                UPVisopyc(i,j,rhoStep)=UPVisopyc(i,j,rhoStep)/(zmax-zmin);
                                % VPV %
                                % remove part of the cell that is not within zmin - zmax
                                VPVisopyc(i,j,rhoStep)=VPVisopyc(i,j,rhoStep)-VPV(i,j,kmin+1)*(zmin-dNEM(kmin));
                                VPVisopyc(i,j,rhoStep)=VPVisopyc(i,j,rhoStep)-VPV(i,j,kmax)*(dNEM(kmax)-zmax);
                                % average over zmin - zmax
                                VPVisopyc(i,j,rhoStep)=VPVisopyc(i,j,rhoStep)/(zmax-zmin);
                                % CT %
                                % remove part of the cell that is not within zmin - zmax
                                CTisopyc(i,j,rhoStep)=CTisopyc(i,j,rhoStep)-ctNEM(i,j,kmin+1)*(zmin-dNEM(kmin));
                                CTisopyc(i,j,rhoStep)=CTisopyc(i,j,rhoStep)-ctNEM(i,j,kmax)*(dNEM(kmax)-zmax);
                                % average over zmin - zmax
                                CTisopyc(i,j,rhoStep)=CTisopyc(i,j,rhoStep)/(zmax-zmin);
                                % SA %
                                % remove part of the cell that is not within zmin - zmax
                                SAisopyc(i,j,rhoStep)=SAisopyc(i,j,rhoStep)-saNEM(i,j,kmin+1)*(zmin-dNEM(kmin));
                                SAisopyc(i,j,rhoStep)=SAisopyc(i,j,rhoStep)-saNEM(i,j,kmax)*(dNEM(kmax)-zmax);
                                % average over zmin - zmax
                                SAisopyc(i,j,rhoStep)=SAisopyc(i,j,rhoStep)/(zmax-zmin);
                                % remove part of the cell that is not within zmin - zmax
                                VORTisopyc(i,j,rhoStep)=VORTisopyc(i,j,rhoStep)-vort(i,j,kmin+1)*(zmin-dNEM(kmin));
                                VORTisopyc(i,j,rhoStep)=VORTisopyc(i,j,rhoStep)-vort(i,j,kmax)*(dNEM(kmax)-zmax);
                                % average over zmin - zmax
                                VORTisopyc(i,j,rhoStep)=VORTisopyc(i,j,rhoStep)/(zmax-zmin);
                            end
                        end
                    end
                end
            end

            disp("Isopycnal calculation completed.")

            %% ISOPYCNAL LNAM CALCULATION - COMMENT OUT AND IN SAVE IF NOT REQUIRED
            for k=1:rhoLen
                %----------------------------------------
                % Calculation of finite spatial element
                Y = latsNEM(SLnem:NLnem-1,1);
                X = lonsNEM(WLnem:ELnem-1,1);
                uNAM = UPVisopyc(1:length(X),1:length(Y),k);
                vNAM = VPVisopyc(1:length(X),1:length(Y),k);
                Xgrid = zeros(length(X),length(Y));
                Ygrid = zeros(length(X),length(Y));
                for p = 1:length(Y)
                    Xgrid(:,p) = X(:,1);
                end
                for p = 1:length(X)
                    Ygrid(p,:) = Y(:,1);
                end
                dx  = zeros(size(Xgrid));
                dy = zeros(size(Xgrid));
                dux = zeros(size(Xgrid));
                duy = zeros(size(Xgrid));
                dvx = zeros(size(Xgrid));
                dvy = zeros(size(Xgrid));

                %----------------------------------------
                %Spatial element in deg if grid_ll==1 or in km otherwise
                dx(2:end-1,2:end-1) = Xgrid(3:end,2:end-1) - Xgrid(1:end-2,2:end-1); %#ok<*COLND>
                dy(2:end-1,2:end-1) = Ygrid(2:end-1,3:end) - Ygrid(2:end-1,1:end-2);


                % define constants
                earth_radius = 6378.137; % km
                % kilometer (km) per degree of latitude
                R = earth_radius*pi/180; % 111.320m
                % Calcul finite spatial element in km
                dx = dx*R.*cosd(Ygrid);
                dy = dy*R;

                % in meters
                dx = dx*1000; % m
                dy = dy*1000; % m

                %Acquire mask and b variables
                mask = isnan(uNAM);
                if isempty(b)
                    b = AMEDA_fetchBfromRD("global_Rossby_Radius.mat",Xgrid,Ygrid,dx,mask);
                end

                %----------------------------------------
                % Compute speed element in m/s
                dux(2:end-1,2:end-1) = (uNAM(2:end-1,3:end) - uNAM(2:end-1,1:end-2));
                duy(2:end-1,2:end-1) = (uNAM(3:end,2:end-1) - uNAM(1:end-2,2:end-1));
                dvx(2:end-1,2:end-1) = (vNAM(2:end-1,3:end) - vNAM(2:end-1,1:end-2));
                dvy(2:end-1,2:end-1) = (vNAM(3:end,2:end-1) - vNAM(1:end-2,2:end-1));

                %----------------------------------------
                % Calculation of Okubo-Weiss criteria
                sn = (dux./dx) - (dvy./dy); % shear "cisaillement"
                ss = (dvx./dx) + (duy./dy); % strain "deformation"
                om = (dvx./dx) - (duy./dy); % vorticity "vorticité"

                okubo = sn.^2 + ss.^2 - om.^2; % in s-2

                %----------------------------------------
                % border is a parameter which prevents the constraints
                % to be applied to points too close to the domain boundaries
                % which would result in an index error
                borders = max(b(:)) + 1;

                %----------------------------------------
                % Calculation of LNAM criteria (Local Normalized Angular Momentum)
                % and LOW criteria (Local Averaged Okubo Weiss)
                %----------------------------------------

                % Initialisation
                if k == 1
                    Lisopyc   = zeros([size(uNAM), rhoLen]);
                    LOWisopyc = nan([size(uNAM), rhoLen]);
                end

                %----------------------------------------
                % calculate LNAM and LOW in all domain pixels
                for i=borders:length(vNAM(:,1))-borders+1
                    for ii=borders:length(vNAM(1,:))-borders+1

                        if ~isnan(vNAM(i,ii))

                            % calculate LOW
                            OW = okubo(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii));
                            LOWisopyc(i,ii,k) = mean(OW(:));

                            % calculate LNAM
                            xlocal = Xgrid(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % Xgrid square sample of length b
                            ylocal = Ygrid(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % Ygrid square sample of length b

                            ulocal = uNAM(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % u square sample of length b
                            vlocal = vNAM(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii)); % v square sample of length b

                            % Local Normalized Angular Momentum
                            % Use the midddle of the square as center coordinate

                            coordcentre = size(ulocal,1)-b(i,ii);

                            d_xcentre = (xlocal - xlocal(coordcentre,coordcentre));
                            d_ycentre = (ylocal - ylocal(coordcentre,coordcentre));

                            cross   = (d_xcentre.*vlocal) - (d_ycentre.*ulocal);
                            dot     = (ulocal.*d_xcentre) + (vlocal.*d_ycentre);
                            produit = sqrt(ulocal.^2 + vlocal.^2)...
                                .*sqrt(d_xcentre.^2 + d_ycentre.^2);
                            sumdp = sum(dot(:))+sum(produit(:));

                            if sumdp ~= 0
                                Lisopyc(i,ii,k) = sum(cross(:)) / sumdp * sign(fNEM(i,ii,depthlevel));
                            else
                                Lisopyc(i,ii,k) = 0;
                            end

                        end
                    end
                end

                thisL = Lisopyc(:,:,k);

                thisL(isnan(thisL)) = 0;

                thisL(mask) = NaN;

                Lisopyc(:,:,k) = thisL;
            end

            disp("Isopycnal LNAM computed.")

            %% Arrange to save in closed-contour algorithm format
            Metadata=struct('DATA','NEMOv3.1','Northern_Limit',NL,'Southern_Limit',SL,'Western_Limit',WL,'Eastern_Limit',EL,'Type_of_SSH','PV');
            filename_out = [output_dir + "/pv_uv_" + datecat];
            date_num=datenum(datecat,'yyyymmdd'); %MATLAB datenum format
            U = uNEM(1:length(X),1:length(Y),:);
            V = vNEM(1:length(X),1:length(Y),:);
            U_ISO = UPVisopyc(1:length(X),1:length(Y),:);
            V_ISO = VPVisopyc(1:length(X),1:length(Y),:);
            EKE = (U.^2+V.^2)./2;
            EKE_ISO = (U_ISO.^2+V_ISO.^2)./2;
            PV_ISO = PVisopyc;
            PV = PV(:,:,:);
            VORT = vort(:,:,:);
            VORT_ISO = VORTisopyc;
            RHO = rhoiso(:,:,depthlevel);
            T_ISO = CTisopyc(1:length(X),1:length(Y),:);
            S_ISO = SAisopyc(1:length(X),1:length(Y),:);
            ZOS = zosNEM(1:length(X),1:length(Y));
            OW = okubo;
            LNAM = L;
            LNAM_ISO = Lisopyc;
            LOW = LOW;
            LOW_ISO = LOWisopyc;
            F_ISO = fcoriopv;
            % Save the .mat file; go to next loop.
            save(filename_out,'Metadata','X','Y','U','V','EKE','EKE_ISO','PV','VORT','RHO','OW','LOW','LNAM',...
                'PV_ISO', 'VORT_ISO', 'T_ISO', 'S_ISO', 'U_ISO', 'V_ISO', 'ZOS', 'F_ISO', 'date_num',...
                'LNAM_ISO','LOW_ISO');
            disp(datecat + " Saved.")

            %% Calculate stats
            %Create matrices if this is at 0 (first file)
            if statsCount == 0
                VORT_CONCAT = nan(size(VORT));
                OW_CONCAT = nan(size(OW));
            end

            %Concat matrices every 10 days so memory doesn't explode but you do you
            if (mod(statsCount, 10) == 1) 
                VORT_CONCAT = [VORT_CONCAT VORT];
                OW_CONCAT = [OW_CONCAT OW];
            end

            %Increment this guy
            statsCount = statsCount + 1;

            toc
        end %END DAYS
    end %END MONTHS
end %END YEARS

%Save back stats: means
VORTSTATS(1,1) = nanmean(nanmean(VORT_CONCAT));
OWSTATS(1,1) = nanmean(nanmean(OW_CONCAT));
%Save back stats: stds
VORTSTATS(1,2) = nanstd(reshape(VORT_CONCAT, [1, numel(VORT_CONCAT)]));
OWSTATS(1,2) = nanstd(reshape(OW_CONCAT, [1, numel(OW_CONCAT)]));

%Save stats into ConfigFile_filltime
save(savedest, 'OWSTATS', 'VORTSTATS', 'statsCount', '-append');
close(h)

