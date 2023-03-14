%%
function [Anticyclonic_Cell,Cyclonic_Cell,Nanti,Ncyclo] = extraction_eddies(X,Y,SSH,minDist,min_amp,min_R,Extraction_Type)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% EXTRACTION DES TOURBILLONS PAR LA METHODE DES CONTOURS CONCENTRIQUES     %
%                             FERMES.                                      %
%                                                                          %
% 1) on identifie les centres possibles en cherchant les minima / maxima   %
%    locaux;                                                               %
% 2) on cherche pour chaque minimum si il ya des contours fermés tous les  %
%    1 mm; le contour le plus externe correspond au bord du tourbillon     %
% 3) si l'amplitude du tourbillon est superieure a un certain critere      %
%    (min_amp), on le garde                                                %
% 4) ce code n'est applicable qu'en approximation geostrophique, pour des  %
%    latitudes > 5 degree                                                  %
%                                                                          %
%  inputs:  X,Y = Vecteurs de coordonnees                                  %
%           SSH = sea level anomalies en m                                 %
%           minDist = Distance minimum (en points de grille) ou l'on       %
%                     identifie les centres possibles (choisir ~150 km)    %
%           min_amp = amplitude minimum que l'on conserve pour les         %
%                     tourbillons                                          %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create matrices of coordinate
[yy,xx] = meshgrid(Y,X);

% Look for non NAN local maximum (Anticylone) and minimum (Cyclone)
% Store them in a sorted vectors f without NaN

f1 = localMaximum(SSH,minDist);  f1(isnan(SSH(f1)))=[];
f2 = localMaximum(-SSH,minDist);  f2(isnan(SSH(f2)))=[];
f=sort([f1(:);f2(:)]);

% Initialize number of eddies at 0
Ncyclo = 0;
Nanti = 0;

thresh = 0.001;

indbad = 1;
%Time step on Local maximum
count = 1:length(f);
while isempty(indbad)==0
    indbad = [];
    
    for i=count %1:length(f)
        indbad1(1:2) = 0;
        % hd is the local SSH maximum considered
        hd = SSH(f(i));
        
        % Find indices of longitudes and lattitude around maximum +/- 5
        % degrees. For vectors (needed for contours function)
        % and in the matrices
        indlon = find(X>=xx(f(i))-5 & X<=xx(f(i))+5);
        indlat = find(Y<=yy(f(i))+5 & Y>=yy(f(i))-5);
        indf =  find(yy(f)<=yy(f(i))+5 & yy(f)>=yy(f(i))-5 & xx(f)>=xx(f(i))-5 & xx(f)<=xx(f(i))+5);
        
        ssh = SSH(indlon,indlat);
        xxx = xx(indlon,indlat);
        yyy = yy(indlon,indlat);
        
        
        % Dichotomy algorithm to find contour of eddies
        
        % Anticyclones
        
        % Set initial amplitude differences with center of 2 m (fin-deb)
        % finish2=fin at beginnig. IsGood is defined=0 as no contour has be
        % identified yet.
        deb = 0;
        fin = 1 + (2*rand(1)-1)*thresh;
        finish2 = fin;
        IsGood = 0;
        
        % iind will be used to store the anplitude and index of identified
        % contour corresponding to criterion in the area studied
        iind = [];
        
        % While loop on contours espaced by 1 mm
        while fin-deb>thresh
            
            % Store in c all contour with values hd-fin in the perimeter
            % defined presviously indicated by indlat & indlon
            c=contourc(X(indlon),Y(indlat),SSH(indlon,indlat)',[hd-fin hd-fin]);
            
            % If there is at least one contour detected
            if isempty(c)==0
                
                % Split contours in a structure array of contours
                [cstruct c] = splitcontours(c);
                
                % Count the number of contours stored in the structure
                NbreC = length(cstruct);
                
                % Loop on these contours
                for j=1:NbreC
                    
                    % Store contours X and Y postion in vectors xd and yd
                    xd = cstruct(j).x;
                    yd = cstruct(j).y;
                    
                    % If it is a closed contour
                    if xd(end) == xd(1) && yd(end) ==yd(1)
                        
                        % If the local maximum lied in this contour and if
                        % there is only on loca maximum using mex inpoly.c
                        if inpoly([xx(f(i));yy(f(i))],[xd(:) yd(:)]')==1 && length(find(inpoly([xx(f(indf))';yy(f(indf))'],[xd(:) yd(:)]'))==1)==1
                            
                            iiind = find(inpoly([xxx(:) yyy(:)]',[xd(:) yd(:)]')==1);
                            if all(isfinite(ssh(iiind)))==1 & length(isfinite(ssh(iiind)))>=4
                                % Give flag of 1 to the variable is good that
                                % saying that a contour corresponds to defined
                                % criterions
                                IsGood = 1;
                                iind = [fin j];
                            end
                            break
                        end
                    end
                end
            end
            
            % Condition on the detection of contour corresponding to criterions
            if IsGood == 1 % One contour was detected
                
                % Divide and shift the area of research toward exterior
                % (toward fin). ie Divide distance between contour and
                % exterior by 2. Reset IsGood to 0.
                deb = fin;
                fin = deb + (finish2-deb)/2;
                IsGood = 0;
                
            elseif IsGood == 0 % No contour was detected
                
                % Divide and shift the area of research toward interior
                % (toward center). ie Divide distance between contour and
                % center by 2.
                finish2 = fin;
                fin = deb + (finish2-deb)/2;
            end
        end
        
        % If a contour stored in ind was detected with amplitude superior
        % than the minimum defined when the function was called
        if isempty(iind)==0 && iind(1)>=min_amp
            
            % Find contours delimiting the closed area around one center
            % where the difference of amplitude with the local maximum is
            % maximum and that correspond to previous criterions
            c=contourc(X(indlon),Y(indlat),SSH(indlon,indlat)',[hd-iind(1) hd-iind(1)]);
            
            % Split contours and take the contours values identified by the
            % iind(2) which store index of good contour (j)
            [cstruct ~] = splitcontours(c);
            xd = cstruct(iind(2)).x;
            yd = cstruct(iind(2)).y;
            
            
        else
            indbad1(1) = 1;
            
        end
        
        % Cyclones
        
        % Set initial amplitude differences with center of 200 cm (fin-deb)
        % finish2=fin at beginnig. IsGood is defined=0 as no contour has be
        % identified yet.
        deb = 0;
        fin = 1 + (2*rand(1)-1)*thresh;
        finish2 = fin;
        IsGood = 0;
        
        % iind will be used to store the anplitude and index of identified
        % contour corresponding to criterion in the area studied
        iind = [];
        
        % While loop on contours espaced by 1 mm
        while fin-deb>thresh
            
            % Store in c all contour with values hd-fin in the perimeter
            % defined presviously indicated by indlat & indlon
            c=contourc(double(X(indlon)),double(Y(indlat)),SSH(indlon,indlat)',[hd+fin hd+fin]);
            
            % If there is at least one contour detected
            if isempty(c)==0
                
                % Split contours in a structure array of contours
                [cstruct ~] = splitcontours(c);
                
                % Count the number of contours stored in the structure
                NbreC = length(cstruct);
                
                % Loop on these contours
                for j=1:NbreC
                    
                    % Store contours X and Y postion in vectors xd and yd
                    xd = cstruct(j).x;
                    yd = cstruct(j).y;
                    
                    % If it is a closed contour
                    if xd(end) == xd(1) && yd(end) ==yd(1)
                        
                        % If the local maximum lied in this contour and if
                        % there is only on loca maximum (Using mex
                        if inpoly([xx(f(i));yy(f(i))],[xd(:) yd(:)]')==1 && length(find(inpoly([xx(f(indf))';yy(f(indf))'],[xd(:) yd(:)]'))==1)==1
                            iiind = find(inpoly([xxx(:) yyy(:)]',[xd(:) yd(:)]')==1);
                            if all(isfinite(ssh(iiind)))==1 & length(isfinite(ssh(iiind)))>=4
                                
                                % Give flag of 1 to the variable is good that
                                % saying that a contour corresponds to defined
                                % criterions
                                IsGood = 1;
                                iind = [fin j];
                                break
                            end
                        end
                    end
                end
            end
            
            % Condition on the detection of contour corresponding to criterions
            if IsGood == 1 % One contour was detected
                
                % Divide and shift the area of research toward exterior
                % (toward fin). ie Divide distance between contour and
                % exterior by 2. Reset IsGood to 0.
                deb = fin;
                fin = deb + (finish2-deb)/2;
                IsGood = 0;
                
            elseif IsGood == 0 % No contour was detected
                
                % Divide and shift the area of research toward interior
                % (toward center). ie Divide distance between contour and
                % center by 2.
                finish2 = fin;
                fin = deb + (finish2-deb)/2;
            end
        end
        
        % If a contour stored in ind was detected with amplitude superior
        % than the minimum defined when the function was called
        if isempty(iind)==0 && iind(1)>=min_amp
            
            % Find contours delimiting the closed area around one center
            % where the difference of amplitude with the local maximum is
            % maximum and that correspond to previous criterions
            c=contourc(double(X(indlon)),double(Y(indlat)),SSH(indlon,indlat)',[hd+iind(1) hd+iind(1)]);
            
            % Split contours and take the contours values identified by the
            % iind(2) which store index of good contour (j)
            [cstruct ~] = splitcontours(c);
            xd = cstruct(iind(2)).x;
            yd = cstruct(iind(2)).y;
            
        else
            indbad1(2) = 1;
            
        end
        
        if indbad1(1)==1 & indbad1(2)==1
            indbad = [indbad i];
            break
        end
    end
    
    f(indbad) = [];
count = i:length(f);
end







% Create Cells for Cyclone and Anticylone. Row dimension large but will be
% reduced when cells will be filled
if Extraction_Type==1
Cyclonic_Cell=cell(length(f),10);
Anticyclonic_Cell=cell(length(f),10);
else
Cyclonic_Cell=cell(length(f),17);
Anticyclonic_Cell=cell(length(f),17);
end


%Time step on Local maximum
for i=1:length(f)
    
    % hd is the local SSH maximum considered
    hd = SSH(f(i));
    
    % Find indices of longitudes and lattitude around maximum +/- 4
    % degrees. For vectors (needed for contours function)
    % and in the matrices
    indlon = find(X>=xx(f(i))-5 & X<=xx(f(i))+5);
    indlat = find(Y<=yy(f(i))+5 & Y>=yy(f(i))-5);
    indf =  find(yy(f)<=yy(f(i))+5 & yy(f)>=yy(f(i))-5 & xx(f)>=xx(f(i))-5 & xx(f)<=xx(f(i))+5);
    
    ssh = SSH(indlon,indlat);
    xxx = xx(indlon,indlat);
    yyy = yy(indlon,indlat);
    % Binary search algorithm (Dichotomie) to find contour of eddies
    
    % Anticyclones
    
    % Set initial amplitude differences with center of 2 m (fin-deb)
    % finish2=fin at beginnig. IsGood is defined=0 as no contour has be
    % identified yet.
    deb = 0;
    fin = 1+(2*rand(1)-1)*thresh;
    finish2 = fin;
    IsGood = 0;
    
    % iind will be used to store the anplitude and index of identified
    % contour corresponding to criterion in the area studied
    iind = [];
    
    % While loop on contours espaced by 1 mm
    while fin-deb>thresh
        
        % Store in c all contour with values hd-fin in the perimeter
        % defined presviously indicated by indlat & indlon
        c=contourc(X(indlon),Y(indlat),SSH(indlon,indlat)',[hd-fin hd-fin]);
        
        % If there is at least one contour detected
        if isempty(c)==0
            
            % Split contours in a structure array of contours
            [cstruct c] = splitcontours(c);
            
            % Count the number of contours stored in the structure
            NbreC = length(cstruct);
            
            % Loop on these contours
            for j=1:NbreC
                
                % Store contours X and Y postion in vectors xd and yd
                xd = cstruct(j).x;
                yd = cstruct(j).y;
                
                % If it is a closed contour
                if xd(end) == xd(1) && yd(end) ==yd(1)
                    
                    % If the local maximum lied in this contour and if
                    % there is only on loca maximum using mex inpoly.c
                    if inpoly([xx(f(i));yy(f(i))],[xd(:) yd(:)]')==1 && length(find(inpoly([xx(f(indf))';yy(f(indf))'],[xd(:) yd(:)]'))==1)==1
                        iiind = find(inpoly([xxx(:) yyy(:)]',[xd(:) yd(:)]')==1);
                        if all(isfinite(ssh(iiind)))==1 & length(isfinite(ssh(iiind)))>=4
                            
                            % Give flag of 1 to the variable is good that
                            % saying that a contour corresponds to defined
                            % criterions
                            IsGood = 1;
                            iind = [fin j];
                            break
                        end
                    end
                end
            end
        end
        
        % Condition on the detection of contour corresponding to criterions
        if IsGood == 1 % One contour was detected
            
            % Divide and shift the area of research toward exterior
            % (toward fin). ie Divide distance between contour and
            % exterior by 2. Reset IsGood to 0.
            deb = fin;
            fin = deb + (finish2-deb)/2;
            IsGood = 0;
            
        elseif IsGood == 0 % No contour was detected
            
            % Divide and shift the area of research toward interior
            % (toward center). ie Divide distance between contour and
            % center by 2.
            finish2 = fin;
            fin = deb + (finish2-deb)/2;
        end
    end
    
    % If a contour stored in ind was detected with amplitude superior
    % than the minimum defined when the function was called
    if isempty(iind)==0 && iind(1)>=min_amp
        
        % Find contours delimiting the closed area around one center
        % where the difference of amplitude with the local maximum is
        % maximum and that correspond to previous criterions
        c=contourc(X(indlon),Y(indlat),SSH(indlon,indlat)',[hd-iind(1) hd-iind(1)]);
        
        % Split contours and take the contours values identified by the
        % iind(2) which store index of good contour (j)
        [cstruct ~] = splitcontours(c);
        xd = cstruct(iind(2)).x;
        yd = cstruct(iind(2)).y;
        
        [b,a] = calcul_aire_rayon(xd,yd);
        
        %         b=areaint(yd,xd,elips);
        %         a = sqrt(b/pi);
        %
              if a>=min_R
        % Indicates that one more anticylones was detected
        Nanti = Nanti+1;
        
        Anticyclonic_Cell{Nanti,8}=b;
        Anticyclonic_Cell{Nanti,7}=a;
        
        % Fill values of X and Y center
        %Anticyclonic_Cell{Nanti,1}=xx(f(i));
        %Anticyclonic_Cell{Nanti,2}=yy(f(i));
        
        % Compute centroid of selected contour
        geom = polygeom(xd,yd);
        
        % Fill values of X and Y  geoïd
        Anticyclonic_Cell{Nanti,3}=geom(2);
        Anticyclonic_Cell{Nanti,4}=geom(3);
        
        % Fill tables of X and Y contours
        Anticyclonic_Cell{Nanti,5}=xd;
        Anticyclonic_Cell{Nanti,6}=yd;
        
        % Fill values of Amplitude as fin = h(center)-h(contour)
        Anticyclonic_Cell{Nanti,9}=iind(1);
              end
        
    end
    
    % Cyclones
    
    % Set initial amplitude differences with center of 200 cm (fin-deb)
    % finish2=fin at beginnig. IsGood is defined=0 as no contour has be
    % identified yet.
    deb = 0;
    fin = 1 + (2*rand(1)-1)*thresh;
    finish2 = fin;
    IsGood = 0;
    
    % iind will be used to store the anplitude and index of identified
    % contour corresponding to criterion in the area studied
    iind = [];
    
    % While loop on contours espaced by 1 mm
    while fin-deb>thresh
        
        % Store in c all contour with values hd-fin in the perimeter
        % defined presviously indicated by indlat & indlon
        c=contourc(double(X(indlon)),double(Y(indlat)),SSH(indlon,indlat)',[hd+fin hd+fin]);
        
        % If there is at least one contour detected
        if isempty(c)==0
            
            % Split contours in a structure array of contours
            [cstruct ~] = splitcontours(c);
            
            % Count the number of contours stored in the structure
            NbreC = length(cstruct);
            
            % Loop on these contours
            for j=1:NbreC
                
                % Store contours X and Y postion in vectors xd and yd
                xd = cstruct(j).x;
                yd = cstruct(j).y;
                
                % If it is a closed contour
                if xd(end) == xd(1) && yd(end) ==yd(1)
                    
                    % If the local maximum lied in this contour and if
                    % there is only on loca maximum (Using mex
                    if inpoly([xx(f(i));yy(f(i))],[xd(:) yd(:)]')==1 && length(find(inpoly([xx(f(indf))';yy(f(indf))'],[xd(:) yd(:)]'))==1)==1
                        iiind = find(inpoly([xxx(:) yyy(:)]',[xd(:) yd(:)]')==1);
                        if all(isfinite(ssh(iiind)))==1 & length(isfinite(ssh(iiind)))>=4
                            
                            % Give flag of 1 to the variable is good that
                            % saying that a contour corresponds to defined
                            % criterions
                            IsGood = 1;
                            iind = [fin j];
                            break
                        end
                    end
                end
            end
        end
        
        % Condition on the detection of contour corresponding to criterions
        if IsGood == 1 % One contour was detected
            
            % Divide and shift the area of research toward exterior
            % (toward fin). ie Divide distance between contour and
            % exterior by 2. Reset IsGood to 0.
            deb = fin;
            fin = deb + (finish2-deb)/2;
            IsGood = 0;
            
        elseif IsGood == 0 % No contour was detected
            
            % Divide and shift the area of research toward interior
            % (toward center). ie Divide distance between contour and
            % center by 2.
            finish2 = fin;
            fin = deb + (finish2-deb)/2;
        end
    end
    
    % If a contour stored in ind was detected with amplitude superior
    % than the minimum defined when the function was called
    if isempty(iind)==0 && iind(1)>=min_amp
        
        % Find contours delimiting the closed area around one center
        % where the difference of amplitude with the local maximum is
        % maximum and that correspond to previous criterions
        c=contourc(double(X(indlon)),double(Y(indlat)),SSH(indlon,indlat)',[hd+iind(1) hd+iind(1)]);
        
        % Split contours and take the contours values identified by the
        % iind(2) which store index of good contour (j)
        [cstruct ~] = splitcontours(c);
        xd = cstruct(iind(2)).x;
        yd = cstruct(iind(2)).y;
        
        [b,a] = calcul_aire_rayon(xd,yd);
        %b=areaint(yd,xd,elips);
        %a = sqrt(b/pi);
        
              if a>=min_R
        
        % Indicates that one more anticylones was detected
        Ncyclo = Ncyclo+1;
        
        Cyclonic_Cell{Ncyclo,8}=b;
        Cyclonic_Cell{Ncyclo,7}=a;
        
        % Compute centroid of selected contour
        geom = polygeom(xd,yd);
        
        % Fill values of X and Y center
        % Cyclonic_Cell{Ncyclo,1}=xx(f(i));
        % Cyclonic_Cell{Ncyclo,2}=yy(f(i));
        
        % Fill values of X and Y  geoïd
        Cyclonic_Cell{Ncyclo,3}=geom(2);
        Cyclonic_Cell{Ncyclo,4}=geom(3);
        
        % Fill tables of X and Y contours
        Cyclonic_Cell{Ncyclo,5}=xd;
        Cyclonic_Cell{Ncyclo,6}=yd;
        
        % Fill values of Amplitude as fin = h(center)-h(contour)
        Cyclonic_Cell{Ncyclo,9}=iind(1);
              end
    end
    
    %Reduce cells to their non-NaN extend
    Anticyclonic_Cell(Nanti+1:end,:) = [];
    Cyclonic_Cell(Ncyclo+1:end,:) = [];
    
end


