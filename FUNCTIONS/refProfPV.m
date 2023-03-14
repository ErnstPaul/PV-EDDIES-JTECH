function [kMaxRhoRef, rhoRef, zRef] = refProfPV(refLon, refLat, refyear, refmonth, refday)
%% Reference profile calculation (for rescaling aspect of rescaled PV)
%LH:  refday = 29; refmonth = 12; refyear = 3; refLon = 73.625; refLat = 8.625;
%GW:  refday = 23; refmonth = 8; refyear = 4; refLon = 53.625; refLat = 8.625;
% Reference Date + Reference Location:
% refday = 22; refmonth = 9; refyear = 4; refLon = 72; refLat = 0;
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
% Paul_OPT_LoadDefaultSettings(); %See corresponding function
disp(["Calculating reference profile for month " + num2str(refmonth)])
refLonNEM = nearestpoint(refLon, lonsNEM); refLatNEM = nearestpoint(refLat, latsNEM);
% Grab reference profile data
depthlevel = 1; dZ = 50;
file2use = input_final_NEM{refday,refmonth,refyear};
ptemps = ncread(file2use, 'thetao',...
    [refLonNEM,refLatNEM,depthlevel,1], [1, 1, dZ, 1]);
sals = ncread(file2use, 'so',...
    [refLonNEM,refLatNEM,depthlevel,1], [1, 1, dZ, 1]);
zRef = double(ncread(file2use, 'depth', depthlevel, dZ));
% Pressure from Depth
pNEMgrid = gsw_p_from_z(zRef.*-1,refLat);
% Conservative T
ctNEM = gsw_CT_from_t(squeeze(sals),squeeze(ptemps),pNEMgrid);
% Density
rhoRefOrig = gsw_rho(squeeze(sals),squeeze(ctNEM),0);
% Grab maximum
kMaxRhoRef=find(rhoRefOrig==max(rhoRefOrig));
% Slot back in
rhoRef=rhoRefOrig(1:kMaxRhoRef);
%% Plot
% plot(zRef,rhoRefOrig, 'LineWidth', 6)
% camroll(270)
% xlim([0 5000])
% xlabel("Depth")
% ylabel("Density (kg/m^3)")
% saveas(gcf, strcat("RefProf"+outputDir,num2str(depthlevel)))
end