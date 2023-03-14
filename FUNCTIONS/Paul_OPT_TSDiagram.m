%% Paul_OPT_TSDiagram.m (version 1.1)
%Author: Paul Ernst
%Date Created: 6/15/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared for final revision
%PE 6/15/2022 - Created
%--------------------------------------
%Purpose: Returns conservative temp and salinity back for plotting in a ts diagram
%Inputs: latlonbounds and times to grab from NEMOv3.1
%Outputs: X by 50 matrices of both CT and SA where X is the number of profiles
function [CT2Plot,SA2Plot] = Paul_OPT_TSDiagram(latlonbounds, yearsMAT, monthsMAT, daysMAT)

%Basic inputs
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);
NLnem = nearestpoint(NL, latsNEM);   SLnem = nearestpoint(SL, latsNEM);
ELnem = nearestpoint(EL, lonsNEM);   WLnem = nearestpoint(WL, lonsNEM);
dZ = 50;
depthlevel = 1;

%begin count
count = 0;
steps = length(yearsMAT)*length(monthsMAT)*length(daysMAT);

%arrange full count calculation
file2use = input_final_NEM{1,1,1};
ptemps = ncread(file2use, 'thetao',...
    [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
sals = ncread(file2use, 'so',...
    [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
[xpt, ypt, zpt] = size(ptemps);
ptempsReshape = reshape(ptemps,[xpt*ypt,zpt]);
fullCount = xpt*ypt*steps;

%allocate full result matrices
CT2Plot = zeros(fullCount,50);
SA2Plot = zeros(fullCount,50);
%% Enable file reading
h = waitbar(0,'Extraction of Variables from NEMOv3.1');
%Filter for usable files in the input directory
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

            disp("Extracting: " + datecat + " [" + file2use + "]");
            % Load relevant variables from the file (T, S, U, V, Z)
            % Reads only the specified section of file and the depth levels above/below for PV
            % gradient calculations
            ptemps = ncread(file2use, 'thetao',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
            sals = ncread(file2use, 'so',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem+1, NLnem-SLnem+1, dZ, 1]);
            dNEM = transpose(double(ncread(file2use, 'depth', depthlevel, dZ)));
            sizemat = size(ptemps);

            % Convert depth to pressure across the entire latitude band
            pNEMgrid = zeros((NLnem-SLnem+1),length(dNEM));
            for latLoop = 1:(NLnem-SLnem+1)
                pNEMgrid(latLoop,:) = gsw_p_from_z(dNEM.*-1,latsNEM(SLnem+latLoop,1));
            end

            saNEM = sals;

            % Convert temperature to conservative temperature across the entire latitude band
            ctNEM = nan(sizemat);
            for latLoop = 1:(NLnem-SLnem+1)
                ctNEM(:,latLoop,:) = gsw_CT_from_t(squeeze(saNEM(:,latLoop,:)),squeeze(ptemps(:,latLoop,:)),pNEMgrid(latLoop,:));
            end

            ctReshape = reshape(ctNEM,[xpt*ypt,zpt]);
            saReshape = reshape(saNEM,[xpt*ypt,zpt]);
            thisStep = (count-1)*xpt*ypt+1;
            nextStep = (count)*xpt*ypt;
            CT2Plot(thisStep:nextStep,:) = ctReshape;
            SA2Plot(thisStep:nextStep,:) = saReshape;
            toc
        end
    end
end
