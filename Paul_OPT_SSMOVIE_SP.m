%% Paul_OPT_SSMOVIE_SP.m (version 2.2)
%Author: Paul Ernst
%Date Created: 11/16/2021
%Date of Last Update: 1/12/2023
%Update History:
%PE 1/12/2023 - Cleared in final revision
%PE 2/9/2022 - Updated name
%PE 11/16/2021 - Created
%--------------------------------------
%Purpose: Creates movie of VORT + IDs for AS. (Movie S2)
%Inputs: VORT, U, V, X, Y data from extracted files; Eddy tracking from Extraction
%Outputs: Full movie of Spiciness + associated IDs at a given depth level
%--------------------------------------
%% Important dynamic inputs
latlonbounds = [30, -10, 80, 40]; % [N, S, E, W] lat long boundaries
depthlevel = 2;
qs = 6;
titlestring = 'ErnstMovieS2.mp4';

%% Load variables & Inputs
%Defaults
load('DefaultData_OPT.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

%Load colormap and bounds 
load('whiteorangecmap.mat')
NL = latlonbounds(1);   SL = latlonbounds(2); WL = latlonbounds(4);   EL = latlonbounds(3);

%Input directories and configfile paths.
input_dir_config= strcat(basepath, '/EXTRACTION/INDIVIDUAL_FILES/NEM_FINAL1');
input_dir_config_traj = strcat(basepath, '/EXTRACTION/EDDY_TRAJECTORIES/NEM_FINAL2/2/');
input_dir_config_cont = strcat(basepath, '/EXTRACTION/EDDY_PROPERTIES/NEM_FINAL2/2/');
loadfile = strcat(input_dir_config, '/ConfigFile_fulltime.mat');
addpath(input_dir_config_traj) %where are the trajectories files located?
addpath(input_dir_config) %where are the files located?
load(strcat(input_dir_config_traj, '/CE_Filtered_Trajectories.mat')) %load CE trajectories
load(strcat(input_dir_config_traj, '/AE_Filtered_Trajectories.mat')) %load AE trajectories
load(strcat(input_dir_config_cont, '/Eddy_Contours.mat'))
load(loadfile);  %ConfigFile_fulltime.mat;
list = dir([ input_dir_config '/' 'pv*.mat']);
filename = [input_dir_config '/' list(1).name];
load(filename)

%Configure output
outdirinterim = ['/MOVIE/FRAMES/FINALMOVIESPIC/' num2str(depthlevel) '/'];
output_dir = strcat(basepath, outdirinterim);
[status,message,messageid] = mkdir(output_dir);

%Rename this variable
date_num_t = date_num;

%% Construct grid.
Xgrid = zeros(length(X),length(Y));
Ygrid = zeros(length(X),length(Y));
for p = 1:length(Y)
    Xgrid(:,p) = X(:,1);
end
for p = 1:length(X)
    Ygrid(p,:) = Y(:,1);
end
%% Filter for all eddies in the AS
%For CE or CE, simply ctrl+f and replace all CE with CE or vice versa
%Not advised to plot both, gets too hectic with too many numbers
%But you do you, future friend

%grab sublists of x, y, etc
idlist = cell(length(CE_traj),1);
x_a=CE_traj(:,2);
y_a=CE_traj(:,3);
%loop through every identified eddy
for i = 1:length(CE_traj)
    x_1=x_a{i,1};
    y_1=y_a{i,1};
    % is in our spatial domain at any time in its life?
    for j = 1:length(x_1)
        if ((x_1(j) > latlonbounds(4)) && (x_1(j) < latlonbounds(3)) && (y_1(j) > latlonbounds(2)) && (y_1(j) < latlonbounds(1)))
            %grab id of this guy!
            workid = CE_traj(:,1);
            idlist{i,1} = workid{i,1};
            break
        end
    end
end

%grab sublists of x, y, etc
idlist2 = cell(length(AE_traj),1);
x_a=AE_traj(:,2);
y_a=AE_traj(:,3);
%loop through every identified eddy
for i = 1:length(AE_traj)
    x_1=x_a{i,1};
    y_1=y_a{i,1};
    % is in our spatial domain at any time in its life?
    for j = 1:length(x_1)
        if ((x_1(j) > latlonbounds(4)) && (x_1(j) < latlonbounds(3)) && (y_1(j) > latlonbounds(2)) && (y_1(j) < latlonbounds(1)))
            %grab id of this guy!
            workid = AE_traj(:,1);
            idlist2{i,1} = workid{i,1};
            break
        end
    end
end

%remove empty cells real quick
idlist = idlist(~cellfun('isempty',idlist(:,1)),:);
idlist = idlist(~cellfun('isempty',idlist(:,1)),:);

%% Get dates of all Eddies: just need 4 datapoints per eddy across time:
%1: ID of Eddy
%2: Lon of occurrance (2)
%3: Lat of occurrance (3)
%4: Date of occurrance (14)

Edd = cell(4,length(idlist));
x_a=CE_traj(:,4); y_a=CE_traj(:,5); date_a=CE_traj(:,14);
%loop through all ID'd eddies
for i = 1:length(idlist)
    working = idlist{i,1};
    %extract components
    x = x_a{working,1};
    y = y_a{working,1};
    date=date_a{working,1};
    %put all the final things in their appropriate place here
    Edd{1,i} = [idlist{i}];
    Edd{2,i} = double([x]);
    Edd{3,i} = double([y]);
    Edd{4,i} = [date];
end

Edd2 = cell(4,length(idlist2));
x_a=AE_traj(:,4); y_a=AE_traj(:,5); date_a=AE_traj(:,14);
%loop through all ID'd eddies
for i = 1:length(idlist2)
    working = idlist2{i,1};
    %extract components
    x = x_a{working,1};
    y = y_a{working,1};
    date=date_a{working,1};
    %put all the final things in their appropriate place here
    Edd2{1,i} = [idlist2{i}];
    Edd2{2,i} = double([x]);
    Edd2{3,i} = double([y]);
    Edd2{4,i} = [date];
end

%% Get f matrix to subtract everything from
fNEM = zeros(length(X),length(Y));
NLnem = nearestpoint(NL, latsNEM);   SLnem = nearestpoint(SL, latsNEM);
ELnem = nearestpoint(EL, lonsNEM);   WLnem = nearestpoint(WL, lonsNEM);
for latLoop = 1:(NLnem-SLnem)
    % Calculate coriolis parameter across the entire grid
    fNEM(:,latLoop) = gsw_f(latsNEM(SLnem+latLoop,1));
end

%% Loop on files
% Loop on each extracted file (in date order)
for i=[1:length(list)]
    tic
    filename = [input_dir_config '/' list(i).name];
    load(filename, 'S_ISO', 'U_ISO','V_ISO','T_ISO', 'date_num', 'U', 'V')
    if exist('PV_ISO','var')==1
        PV_ISO = gsw_spiciness0(S_ISO(:,:,2),T_ISO(:,:,2));
        if isnan(nansum(PV_ISO(:)))
            continue
        end
    end
    U = double(U_ISO(:,:,depthlevel))*2;
    V = double(V_ISO(:,:,depthlevel))*2;
    figure('units', 'normalized', 'outerposition', [0 0 1 1], 'color', 'w', 'visible', 'off'); %fullscreen, invisible figure
    m_proj('mercator','longitude',[WL EL],'latitude',[SL NL]) %initialize map
    m_pcolor(Xgrid,Ygrid,PV_ISO) 
    hold on
    h1 = m_quiver(Xgrid(1:qs:end,1:qs:end),Ygrid(1:qs:end,1:qs:end),U(1:qs:end,1:qs:end),V(1:qs:end,1:qs:end), 'k');
    set(h1,'AutoScale','off')
    set(h1,'LineWidth',1.2)
    shading interp
    hold on
    m_coast('patch',[.5 .5 .5]);
    m_grid('box', 'fancy','fontsize',36);
    hold on

    %Reference Vector
    diffX = EL-WL; diffY = NL-SL;
    [hpv5, htv5] = m_vec(50, EL-diffY*0.15, NL+diffX*-0.1, 100, 0, 'w', 'key', '1 m s^{-1}',...
        'shaftwidth', 4, 'headwidth', 20, 'headlength', 24);
    set(htv5,'FontSize',28);

    %ID NUMBERS
    for pls = 1:length(idlist)
        if ismember(date_num,Edd{4,pls})
            for n = 1:length(Edd{4,pls})
                if (date_num == Edd{4,pls}(n))
                    indexmat = n;
                    break;
                end
            end
            %plot a text annotation to the right of the eddy center if
            %it's on the correct date and in the correct area
            if ((Edd{2,pls}(indexmat) > WL) && (Edd{2,pls}(indexmat) < (EL-2))...
                    && (Edd{3,pls}(indexmat) > (SL+.5)) && (Edd{3,pls}(indexmat) < NL))
                m_text((Edd{2,pls}(indexmat)),(Edd{3,pls}(indexmat)),int2str(Edd{1,pls}), 'color', 'k', 'fontsize', 12)
            end
        end
    end

    for pls = 1:length(idlist2)
        if ismember(date_num,Edd2{4,pls})
            for n = 1:length(Edd2{4,pls})
                if (date_num == Edd2{4,pls}(n))
                    indexmat = n;
                    break;
                end
            end
            %plot a text annotation to the right of the eddy center if
            %it's on the correct date and in the correct area
            if ((Edd2{2,pls}(indexmat) > WL) && (Edd2{2,pls}(indexmat) < (EL-2))...
                    && (Edd2{3,pls}(indexmat) > (SL+.5)) && (Edd2{3,pls}(indexmat) < NL))
                m_text((Edd2{2,pls}(indexmat)),(Edd2{3,pls}(indexmat)),int2str(Edd2{1,pls}), 'color', 'k', 'fontsize', 12)
            end
        end
    end

    %CONTOURS
    contoursForLoop = CEs(:,i,:);
    numContours = length(contoursForLoop);
    for contCount = 1:numContours
        m_line(contoursForLoop{contCount,1,1},contoursForLoop{contCount,1,2}, 'color', 'k', 'linewi', 2)
    end
    %CONTOURS
    contoursForLoop = AEs(:,i,:);
    numContours = length(contoursForLoop);
    for contCount = 1:numContours
        m_line(contoursForLoop{contCount,1,1},contoursForLoop{contCount,1,2}, 'color', 'k','marker','.', 'linewi', 2)
    end
   

    %this is a section dedicated to the date marking in the top right
    %corner of the plot
    namestart = filename(106:113);
    %namestart = filename(77:84);
    name2use = [namestart(5:6) '/' namestart(7:8) '/' namestart(1:4)];

    %NOTE NOTE NOTE NOTE NOTE NOTE: 70, 28 is the (lon, lat) that this
    %date annotation will be, and this WILL vary on your geographic
    %area. There is no way to detect easily where is the best place for
    %this label-- please adjust it to where you need it manually!
    m_text(68,28,name2use,'color','w','fontsize',40);
    set(gca,'fontsize',36);
    colorbar(gca,'Position',[0.732031250000002 0.113526570048309 0.0194901019588132 0.811594202898551]); colormap(jet);
    %I don't know how this colorbar position varies with screen resolution. You
    %may need to fiddle with this.
    ylabel('Latitude');
    xlabel('Longitude');
    title(['PV_I_S_O (s^-^1) and all IDs At Depth'], 'fontsize', 36);
    caxis([0 2]);
    %This color axis is generally pretty good for Spiciness. But, you can
    %certainly figure out if you have a better one.
    %colormap(CustomColormap) %colormap for the figure);
    set(gcf,'PaperPositionMode','auto');
    %save the figure to a filename in the output directory
    disp(['Frame:' ' ' name2use]);
    filenamestring = ([output_dir namestart '.png']);
    filename = char(filenamestring);
    %we use export_fig because it automatically crops and is nice
    export_fig(filename);
    %we have to delete the figure or else memory go kaboom because
    %MATLAB is a fun and engaging memory management simulator
    delete(gcf)
    disp("Frame " + num2str(i) + " is done in " + toc + " seconds.");
end

%% Make Movie
filedir=dir(output_dir);
addpath(output_dir);
filedir = filedir(~[filedir.isdir]);
[~,idx] = sort([filedir.datenum]);
filedir = filedir(idx);
video=VideoWriter(titlestring,'MPEG-4'); %create the video object, using mp4 by default
video.FrameRate=15;%set video frame rate. number of frames per second. default is 30 frames per second
open(video); %open the file for writing
for dday=1:length(filedir) %where N is the number of images; we use 3 because the only other thing in this directory
    %should be "." and ".." which are the present
    %and current directories. If you have other
    %random stuff in here, then just note that you
    %may need to change "3" to something else. Check
    %the "filedir" variable.
    disp(['processing video:' ' ' filedir(dday).name]);
    I=imread([filedir(dday).name]); %read the next image
    writeVideo(video,I);
end
close(video); %close the file
disp('finished');