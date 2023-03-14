%% Paul_OPT_TSDiagramAlgor.m (version 1.1)
%Author: Paul Ernst
%Date Created: 6/17/2022
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared for final revision
%PE 6/15/2022 - Created
%--------------------------------------
%Purpose: Returns conservative temp and salinity back for plotting in a ts diagram, but bins data
%           (Binned data is averaged across all profiles)
%Inputs: latlonbounds and times to grab from NEMOv3.1
%Outputs: X by 50 matrices of both CT and SA where X is the number of profiles
function [maxRSW,topminRSW,botminRSW] = Paul_OPT_TSDiagramAlgor(latlonbounds, yearsMAT, monthsMAT, daysMAT, lowbound, highbound)

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
    [WLnem,SLnem,depthlevel,1], [ELnem-WLnem, NLnem-SLnem, dZ, 1]);
[xpt, ypt, ~] = size(ptemps);

%allocate full result matrices
maxRSW = zeros(xpt, ypt); validMax = zeros(xpt, ypt);
topminRSW = zeros(xpt, ypt); validMin1 = zeros(xpt, ypt);
botminRSW = zeros(xpt, ypt); validMin2 = zeros(xpt, ypt);
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
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem, NLnem-SLnem, dZ, 1]);
            sals = ncread(file2use, 'so',...
                [WLnem,SLnem,depthlevel,1], [ELnem-WLnem, NLnem-SLnem, dZ, 1]);
            dNEM = transpose(double(ncread(file2use, 'depth', depthlevel, dZ)));
            sizemat = size(ptemps);

            % Convert depth to pressure across the entire latitude band
            pNEMgrid = zeros((NLnem-SLnem),length(dNEM));
            for latLoop = 1:(NLnem-SLnem)
                pNEMgrid(latLoop,:) = gsw_p_from_z(dNEM.*-1,latsNEM(SLnem+latLoop,1));
            end

            saNEM = sals;

            % Convert temperature to conservative temperature across the entire latitude band
            ctNEM = nan(sizemat);
            for latLoop = 1:(NLnem-SLnem)
                ctNEM(:,latLoop,:) = gsw_CT_from_t(squeeze(saNEM(:,latLoop,:)),squeeze(ptemps(:,latLoop,:)),pNEMgrid(latLoop,:));
            end

            %interpolate
            zToInterp = 575:25:5000;
            newCTNEM = interp3(X,Y,depthsNEM,ctNEM,X,Y,zToInterp);
            newSANEM = interp3(X,Y,depthsNEM,saNEM,X,Y,zToInterp);
            sizemat = size(newCTNEM);
            zpt = length(zToInterp);

            %calc rho and spic
            rhoNEM = nan(sizemat);
            spic2plot = nan(sizemat);
            for lonLoop = 1:(ELnem-WLnem)
                rhoNEM(lonLoop,:,:) = gsw_rho(squeeze(newSANEM(lonLoop,:,:)),squeeze(newCTNEM(lonLoop,:,:)),0);
                spic2plot(lonLoop,:,:) = gsw_spiciness0(squeeze(newSANEM(lonLoop,:,:)),squeeze(newCTNEM(lonLoop,:,:)));
            end

            %Calculate density
            densYes = ((rhoNEM >= lowbound) & (rhoNEM <= highbound));
            densYes(:,:,1:nearestpoint(600,zToInterp)) = 0; %depth cutoff

            % PLOT HIGH AND LOW POINTS
            %detrend
            spic2find = spic2plot;
            spic2find(~densYes) = NaN;
            rhoNEM(~densYes) = NaN;
            for xx = 1:xpt
                for yy = 1:ypt
                    x = squeeze(newSANEM(xx,yy,:));
                    %check if all nans, to just keep going
                    if (sum(~isnan(x)) == 0)
                        continue;
                    end
                    x(isnan(x)) = [];
                    y = squeeze(newCTNEM(xx,yy,:)); y(isnan(y)) = [];
                    coefficients = polyfit(x, y, 1);
                    xFit = linspace(min(x), max(x), length(x));
                    yFit = polyval(coefficients , xFit);
                    spic2detrend = gsw_spiciness0(xFit,yFit);
                    if (length(spic2detrend) < length(spic2find(xx,yy,:)))
                        spic2detrend = [spic2detrend zeros(1,length(spic2find(xx,yy,:))-length(spic2detrend))];
                    end
                    spic2find(xx,yy,:) = squeeze(spic2find(xx,yy,:)) - spic2detrend';
                end
            end
            %find max
            [~,maxSpicInd] = nanmax(spic2find,[],3);
            %find min above max
            low2find = spic2find;
            thisInd1 = zeros(xpt,ypt);
            for xx = 1:xpt
                for yy = 1:ypt
                    thisInd1(xx,yy) = min(maxSpicInd(xx,yy),zpt) + 1;
                    low2find(xx,yy,thisInd1(xx,yy):zpt) = NaN;
                end
            end
            [~,minSpicInd] = nanmin(low2find,[],3);
            %find min below max
            low2find2 = spic2find;
            thisInd2 = zeros(xpt,ypt);
            for xx = 1:xpt
                for yy = 1:ypt
                    thisInd2(xx,yy) = max(maxSpicInd(xx,yy),2) - 1;
                    low2find2(xx,yy,1:thisInd2(xx,yy)) = NaN;
                end
            end
            [~,minSpicInd2] = nanmin(low2find2,[],3);

            %Slot into return matrices
            for xx = 1:xpt
                for yy = 1:ypt
                    %Handle Max
                    toAdd = rhoNEM(xx,yy,maxSpicInd(xx,yy));
                    if isnan(toAdd)
                        toAdd = 0;
                    else
                        validMax(xx,yy) = validMax(xx,yy) + 1;
                    end
                    maxRSW(xx,yy) = maxRSW(xx,yy) + toAdd;

                    %Handle min1
                    toAdd = rhoNEM(xx,yy,minSpicInd(xx,yy));
                    if isnan(toAdd)
                        toAdd = 0;
                    else
                        validMin1(xx,yy) = validMin1(xx,yy) + 1;
                    end
                    topminRSW(xx,yy) = topminRSW(xx,yy) + toAdd;
                                      
                    %Handle min2
                    toAdd = rhoNEM(xx,yy,minSpicInd2(xx,yy));
                    if isnan(toAdd)
                        toAdd = 0;
                    else
                        validMin2(xx,yy) = validMin2(xx,yy) + 1;
                    end
                    botminRSW(xx,yy) = botminRSW(xx,yy) + toAdd;
                end
            end
            toc
        end
    end
end

% save("JustInCaseTSDiagramAlgor", "validMax", "validMin1", "validMin2", "maxRSW", "topminRSW", "botminRSW")

%NaNs
validMax(validMax==0) = NaN;
validMin1(validMin1==0) = NaN;
validMin2(validMin2==0) = NaN;
%avg and return
maxRSW = maxRSW./validMax;
topminRSW = topminRSW./validMin1;
botminRSW = botminRSW./validMin2;
