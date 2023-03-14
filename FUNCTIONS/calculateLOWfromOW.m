load('DefaultData_SS.mat') %See Paul_OPT_CreateDefaultSettings.m
Paul_OPT_LoadDefaultSettings(); %See corresponding function

depthlevel = 1;
outputDir = 'NEM_PGW_RSW'; %name of the directory we're going to be using
output_dir_interim = strcat(basepath, ['/EXTRACTION/INDIVIDUAL_FILES/' outputDir]);
input_dir = strcat(output_dir_interim, [num2str(depthlevel) + "/"]);
inputDir = char(input_dir);

%% Enable file reading
%Filter for usable files in the input directory
list = dir([ inputDir '/' 'pv*.mat']);
configFile = dir([ inputDir '/' 'Config*.mat']);
load([configFile.folder '/' configFile.name]);
b = [];

for fileNum = 1:length(list)
    tic
    disp(['Extracting LOW for:    ' list(fileNum).name ''])
    filename = [inputDir  list(fileNum).name];
    %Load and reaffirm presence of all necessary variables
    load(filename, 'U_ISO', 'V_ISO', 'X', 'Y')

    %loop on iso layers
    for k = 1:2
            uNAM = U_ISO(1:length(X),1:length(Y),k);
            vNAM = V_ISO(1:length(X),1:length(Y),k);
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
                LOWisopyc = nan([size(uNAM), 2]);
            end

            %----------------------------------------
            % calculate LNAM and LOW in all domain pixels
            for i=borders:length(vNAM(:,1))-borders+1
                for ii=borders:length(vNAM(1,:))-borders+1

                    if ~isnan(vNAM(i,ii))

                        % calculate LOW
                        OW = okubo(i-b(i,ii):i+b(i,ii),ii-b(i,ii):ii+b(i,ii));
                        LOWisopyc(i,ii,k) = mean(OW(:));
                    end
                end
            end
    end
    LOW_ISO = LOWisopyc;
    save(filename, 'LOW_ISO', '-append')
    disp([num2str(fileNum) + " of " + num2str(length(list))])
    toc
end