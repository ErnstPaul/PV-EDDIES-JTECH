%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%               Filtrage Alexis             %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Originally created, unchanged from Cori Pegliasco's code + Alexis Chaigneau's Code

%% AEs
function Paul_OPT_Step6_SS(inputDir,propDir,trajDir,minlife, minamp, minrad, depthlevel)
%% récupération des vecteurs totaux L, A et R
close all
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function
seuil_lifetime = minlife; 
seuil_Amp = minamp;
seuil_Radius = minrad;
Minimum_Amplitude_Trajectories = seuil_Amp;
Minimum_Radius_Trajectories = seuil_Radius;
Minimum_Duration_Trajectories = seuil_lifetime;
input_dir_config = inputDir;
output_dir= trajDir;
output_dir= strcat(output_dir, num2str(depthlevel));
prop_dir= propDir;
prop_dir= strcat(prop_dir, num2str(depthlevel));
save(strcat(input_dir_config,'/ConfigFile_fulltime.mat'), 'Minimum_Radius_Trajectories', 'Minimum_Amplitude_Trajectories', 'Minimum_Duration_Trajectories','-append');
addpath(strcat(basepath, 'FUNCTIONS'));
input_dir = [strcat(output_dir,'/')]; %EDDY_TRAJECTORIES/'];
output_dir = input_dir; %EDDY_TRAJECTORIES/'];
load([strcat(prop_dir,'/AE_Properties.mat')],'date_num','numtraj'); %EDDY_PROPERTIES/AE_Properties.mat'],'date_num','numtraj');

list = dir([ input_dir 'AE_Traj*']);
load([input_dir list(1).name],'idtraj')

L = NaN([1 nanmax(idtraj(:))*length(list)]);
A = NaN([1 nanmax(idtraj(:))*length(list)]);
R = NaN([1 nanmax(idtraj(:))*length(list)]);

traj = struct('ID_traj',{},'Xcenter',{},'Ycenter',{},...
    'Xcentroid',{},'Ycentroid',{},'Radius',{},'Area',{},'Amplitude',{},...
    'EKE',{},'Speed',{},'Vorticity',{},'OW',{},...
    'Eddy_index',{},'Date',{},'Time_index',{},'Lifetime',{},...
    'Merging_time',{},'Merged_with_traj',{},'ID_eddies',{},'Splitting_from_traj',{},...
    'Xcontour',{},'Ycontour',{},'Splitting_time',{},'Splitted_in',{},'Eddy_found_at',{});

traj(nanmax(idtraj(:))*length(list)).Lifetime = 0;

a = 0;
c = 0;

for i=1:length(list)
    load([input_dir list(i).name])

    A(a+1:a+numtraj)=arrayfun(@(x) nanmedian(double(x.Amplitude)),AEs_traj(1:numtraj));
    R(a+1:a+numtraj)=arrayfun(@(x) nanmedian(double(x.Radius)),AEs_traj(1:numtraj));
    L(a+1:a+numtraj)=arrayfun(@(x) (x.Lifetime),AEs_traj(1:numtraj));
    
    ind = find(L(a+1:a+numtraj) >= seuil_lifetime & A(a+1:a+numtraj) >= seuil_Amp & R(a+1:a+numtraj) >= seuil_Radius);
    if isempty(ind)==0
        traj(c+1:c+length(ind)) = orderfields(AEs_traj(ind),traj);
        c = c+length(ind);
    end
    
    a = a+numtraj;
end
R(a+1:end) = [];
A(a+1:end) = [];
L(a+1:end) = [];
traj(c+1:end) = [];


%% Filtrage
id_traj = find(L >= seuil_lifetime & A >= seuil_Amp & R >= seuil_Radius); %bonne trajs
id_out = find(L < seuil_lifetime | A < seuil_Amp | R < seuil_Radius);  % mauvaises

Ntraj = 5000; % le nbre de traj dans chaque rawXX.mat
champs = [2:5 8 6 7 9 14 15 13 19 21 22 25];
ff = fields(traj);

INDEXES = NaN(3,length(id_traj));
c = 0;
for i=1:length(id_traj)
    if isempty(traj(i).Splitting_from_traj)==0 & isempty(find(id_out == traj(i).Splitting_from_traj))==0
        c = c+1;
        t = ceil(traj(i).Splitting_from_traj/Ntraj); % on cherche la bonne rawXX.mat à charger
        INDEXES(:,c) = [i t traj(i).Splitting_from_traj-(t-1)*Ntraj];
    end
end
INDEXES(:,c+1:end) = [];

IND = unique(INDEXES(2,:));
for k=1:length(IND)
    ind = find(INDEXES(2,:)==IND(k));
    load([input_dir list(IND(k)).name])
    for kk=ind
        % si la trajectoire est issue d'une trajectoire qui ne passe pas le
        % seuil, on efface le splitting from et on recopie les valeurs
        % d'avant le splitting
        i = INDEXES(1,kk);
        t = INDEXES(2,kk);
        
        temps_split = find(AEs_traj(INDEXES(3,kk)).Time_index == traj(i).Time_index(1));
        traj_modif = orderfields(AEs_traj(INDEXES(3,kk)),traj);
        
        %% association des trajectoires
        % champs en vecteurs
        for ii = 1:length(champs) % remplacer apres le temps splitt� par les valeurs de la traj issue du splitting
           traj_modif.(sprintf(ff{champs(ii)}))(temps_split+1:temps_split+traj(i).Lifetime-1) = traj(i).(sprintf(ff{champs(ii)}))(2:end);
        end
          % le reste
        traj_modif.Lifetime = length(traj_modif.Xcenter);
        traj_modif.ID_traj = traj(i).ID_traj;
        traj_modif.Splitting_from_traj = [];
        traj_modif.Splitting_time = [];
        traj_modif.Splitted_in = traj(i).Splitted_in;
        traj_modif.Merging_time = traj(i).Merging_time;
        traj_modif.Merged_with_traj = traj(i).Merged_with_traj;
        
        if traj_modif.Radius>=seuil_Radius & traj_modif.Amplitude>=seuil_Amp
        traj(i) = traj_modif;
        end
    end
end



%% Passage en cellules

clear AEs_traj
AE_traj = cell(length(id_traj),21);

Fields_traj = [  '1.ID_traj                              '; ...
    '2.Xcenter                              '; ...
    '3.Ycenter                              '; ...
    '4.Xcentroid                            '; ...
    '5.Ycentroid                            '; ...
    '6.Equivalent Radius [km]               '; ...
    '7.Eddy Area [km^2]                     '; ...
    '8.Amplitude [cm]                       '; ...
    '9.Mean EKE [(m/s)^2]                   '; ...
    '10.Mean Speed [m/s]                    '; ...
    '11.Mean Vorticity [1/s]                '; ...
    '12.Mean Okubo_Weiss parameter [1/s]    '; ...
    '13.Eddy Index                          '; ...  % 1ere position dans la cell Contours 
    '14.Date in matlab num format           '; ...
    '15.Time index                          '; ...  % 2nde position dans la cell Contours (1 = first map, 2 = second map etc)
    '16.Lifetime [days]                     '; ...
    '17.Merging Time                        '; ...  % en Time Index
    '18.Merged with Trajectory              '; ...  % ID de la trajectoire rejointe
    '19.ID eddies                           '; ...
    '20.Splitting from Trajectory           '; ...  % ID de la trajectoire mère
    '21.ID profile used  (CELL)             '];

AE_traj = squeeze(struct2cell(traj))';
AE_traj(:,22:end) = [];
AE_traj(:,21) = cell(length(id_traj),1);

save([output_dir 'AE_Filtered_Trajectories.mat'],'Metadata','AE_traj','Fields_traj','date_num')

%% CEs

%% récupération des vecteurs totaux L, A et R

load([strcat(prop_dir,'/CE_Properties.mat')],'date_num','numtraj');
list = dir([ input_dir 'CE_Traj*']);
load([input_dir list(1).name],'idtraj')

L = NaN([1 nanmax(idtraj(:))*length(list)]);
A = NaN([1 nanmax(idtraj(:))*length(list)]);
R = NaN([1 nanmax(idtraj(:))*length(list)]);

traj = struct('ID_traj',{},'Xcenter',{},'Ycenter',{},...
    'Xcentroid',{},'Ycentroid',{},'Radius',{},'Area',{},'Amplitude',{},...
    'EKE',{},'Speed',{},'Vorticity',{},'OW',{},...
    'Eddy_index',{},'Date',{},'Time_index',{},'Lifetime',{},...
    'Merging_time',{},'Merged_with_traj',{},'ID_eddies',{},'Splitting_from_traj',{},...
    'Xcontour',{},'Ycontour',{},'Splitting_time',{},'Splitted_in',{},'Eddy_found_at',{});

traj(nanmax(idtraj(:))*length(list)).Lifetime = 0;

a = 0;
c = 0;

for i=1:length(list)
    load([input_dir list(i).name])

    A(a+1:a+numtraj)=arrayfun(@(x) nanmedian(double(x.Amplitude)),CEs_traj(1:numtraj));
    R(a+1:a+numtraj)=arrayfun(@(x) nanmedian(double(x.Radius)),CEs_traj(1:numtraj));
    L(a+1:a+numtraj)=arrayfun(@(x) (x.Lifetime),CEs_traj(1:numtraj));
    
      ind = find(L(a+1:a+numtraj) >= seuil_lifetime & A(a+1:a+numtraj) >= seuil_Amp & R(a+1:a+numtraj) >= seuil_Radius);
    if isempty(ind)==0
        traj(c+1:c+length(ind)) = orderfields(CEs_traj(ind),traj);
        c = c+length(ind);
    end
    
    a = a+numtraj;
end
R(a+1:end) = [];
A(a+1:end) = [];
L(a+1:end) = [];
traj(c+1:end) = [];


%% Filtrage
id_traj = find(L >= seuil_lifetime & A >= seuil_Amp & R >= seuil_Radius); %bonne trajs
id_out = find(L < seuil_lifetime | A < seuil_Amp | R < seuil_Radius);  % mauvaises

Ntraj = 5000; % le nbre de traj dans chaque rawXX.mat
champs = [2:5 8 6 7 9 14 15 13 19 21 22 25];
ff = fields(traj);

INDEXES = NaN(3,length(id_traj));
c = 0;
for i=1:length(id_traj)
    if isempty(traj(i).Splitting_from_traj)==0 & isempty(find(id_out == traj(i).Splitting_from_traj))==0
        c = c+1;
        t = ceil(traj(i).Splitting_from_traj/Ntraj); % on cherche la bonne rawXX.mat à charger
        INDEXES(:,c) = [i t traj(i).Splitting_from_traj-(t-1)*Ntraj];
    end
end
INDEXES(:,c+1:end) = [];

IND = unique(INDEXES(2,:));
for k=1:length(IND)
    ind = find(INDEXES(2,:)==IND(k));
    load([input_dir list(IND(k)).name])
    for kk=ind
        % si la trajectoire est issue d'une trajectoire qui ne passe pas le
        % seuil, on efface le splitting from et on recopie les valeurs
        % d'avant le splitting
        i = INDEXES(1,kk);
        t = INDEXES(2,kk);
        
        temps_split = find(CEs_traj(INDEXES(3,kk)).Time_index == traj(i).Time_index(1));
        traj_modif = orderfields(CEs_traj(INDEXES(3,kk)),traj);
        
        %% association des trajectoires
        % champs en vecteurs
        for ii = 1:length(champs) % remplacer apres le temps splitt� par les valeurs de la traj issue du splitting
           traj_modif.(sprintf(ff{champs(ii)}))(temps_split+1:temps_split+traj(i).Lifetime-1) = traj(i).(sprintf(ff{champs(ii)}))(2:end);
        end
          % le reste
        traj_modif.Lifetime = length(traj_modif.Xcenter);
        traj_modif.ID_traj = traj(i).ID_traj;
        traj_modif.Splitting_from_traj = [];
        traj_modif.Splitting_time = [];
        traj_modif.Splitted_in = traj(i).Splitted_in;
        traj_modif.Merging_time = traj(i).Merging_time;
        traj_modif.Merged_with_traj = traj(i).Merged_with_traj;
        
        if traj_modif.Radius>=seuil_Radius & traj_modif.Amplitude>=seuil_Amp
        traj(i) = traj_modif;
        end
        
        
    end
end

%% Passage en cellules

clear CEs_traj
CE_traj = cell(length(id_traj),21);

Fields_traj = [  '1.ID_traj                              '; ...
    '2.Xcenter                              '; ...
    '3.Ycenter                              '; ...
    '4.Xcentroid                            '; ...
    '5.Ycentroid                            '; ...
    '6.Equivalent Radius [km]               '; ...
    '7.Eddy Area [km^2]                     '; ...
    '8.Amplitude [cm]                       '; ...
    '9.Mean EKE [(m/s)^2]                   '; ...
    '10.Mean Speed [m/s]                    '; ...
    '11.Mean Vorticity [1/s]                '; ...
    '12.Mean Okubo_Weiss parameter [1/s]    '; ...
    '13.Eddy Index                          '; ...  % 1ere position dans la cell Contours 
    '14.Date in matlab num format           '; ...
    '15.Time index                          '; ...  % 2nde position dans la cell Contours (1 = first map, 2 = second map etc)
    '16.Lifetime [days]                     '; ...
    '17.Merging Time                        '; ...  % en Time Index
    '18.Merged with Trajectory              '; ...  % ID de la trajectoire rejointe
    '19.ID eddies                           '; ...
    '20.Splitting from Trajectory           '; ...  % ID de la trajectoire mère
    '21.ID profile used  (CELL)             '];

CE_traj = squeeze(struct2cell(traj))';
CE_traj(:,22:end) = [];
CE_traj(:,21) = cell(length(id_traj),1);

save([output_dir 'CE_Filtered_Trajectories.mat'],'Metadata','CE_traj','Fields_traj','date_num')







