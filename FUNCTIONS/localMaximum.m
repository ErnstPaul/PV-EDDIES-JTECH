function f = localMaximum(x,minDist)
% function f = localMaximum(x,minDist)
% This function returns the indexes of local maximum in the data x.
% x can be a vector or a matrix
% minDist is the minimum distance between two peaks (local maximas)
% minDist should be a vector in which each argument corresponds to it's
% relevent dimension
% Example:
% x = randn(100,30,10);
% minDist = [10 3 5];
% peak = (x,minDist);

ind=find(isnan(x));
x(ind)=0;
if nargin < 2
    minDist = size(x)/10;
end

dimX = length ( size(x) );
if length(minDist) ~= dimX
    % In case minimum distance isn't defined for all of x dimensions
    % I use the first value as the default for all of the dimensions
    minDist = minDist( ones(dimX,1) );
end

% validity checks
minDist = ceil(minDist);
minDist = max( [minDist(:)' ; ones(1,length(minDist))] );
minDist = min( [minDist ; size(x)] );

%--------
% this section comes to solve the problem of a plato
% without this code, points with the same hight will be recognized as peaks
y = sort(x(:));
dY = diff(y);
% finding the minimum step in the data
minimumDiff = min( dY(dY ~= 0) );
%adding noise which won't affect the peaks
x = x + rand(size(x))*minimumDiff;
%--------

se = ones(minDist);
X = imdilate(x,se);
f = find(x == X);

% tic
% toP = X;
% load('/Volumes/Lacie-SAN/SAN2/Paul_Eddies/Subsurface_Optimization/EXTRACTION/INDIVIDUAL_FILES/NEM_FULL1/pv_uv_20160101.mat', 'X')
% load('/Volumes/Lacie-SAN/SAN2/Paul_Eddies/Subsurface_Optimization/EXTRACTION/INDIVIDUAL_FILES/NEM_FULL1/pv_uv_20160101.mat', 'Y')
% Xgrid = zeros(length(X),length(Y));
% Ygrid = zeros(length(X),length(Y));
% for p = 1:length(Y)
%     Xgrid(:,p) = X(:,1);
% end
% for p = 1:length(X)
%     Ygrid(p,:) = Y(:,1);
% end
% figure; pcolor(Xgrid,Ygrid,ZOS); hold on; colormap(jet); shading flat; scatter(Xgrid(f),Ygrid(f),'mx')
% xlim([40 80])
% ylim([-10 30])