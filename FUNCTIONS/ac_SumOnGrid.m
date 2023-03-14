function S = ac_SumOnGrid(Xd,Yd,Data,xg,yg)

% Compute the sum of a parameter on a grid centered at xg yg
%  input: - Xd,Yd = positions of the irregularly distributed data
%         - Data = The Data to be summed
%         - xg,yg = vectors of the grid 
%                 !! All the data (Xd,Yd) MUST be bounded by (xg,yg)
%
% output:  - S = sum of the data on the grid
%----------------------------------------------------------------

    LonBinIndexes = FindBinIndexes(Xd,xg);
    LatBinIndexes = FindBinIndexes(Yd,yg);
    
    S = full(sparse(LonBinIndexes,LatBinIndexes,Data,length(xg),length(yg)));
    
    
    function BinIndexes = FindBinIndexes(X,CenterBins)   
 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% For each X, this functions finds the index of the bin it belongs to in
% CenterBins:  BinIndexes = yt_FindBinIndexes(X,CenterBins) ;
%
%     inputs:  X = data irregularly distributed
%              CenterBins = a SORTED !vector with the centers of the bins
%                            Sorting is not checked for !
%     outputs: BinIndexes = vector the size of X, containing for each X
%     value, the index of the "CenterBins" it belongs to.
%
%   Yann Tremblay (May 2014), based on function factor_scale_for_sparse by
%   Alexis Chaigneau (2006)
%
%       Double precision bug from earlier function solved
%       Can now accomodate irregular grids spacing
%       For use with data in geographic grid space :
%               - cannot accomodate grids with dateline crossing
%               - run this function 2 times : for Lat and for Long
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% No error checking for speed...

% Get half bin intervals
% produce column vector no matter input vector direction
CenterBins=CenterBins(:);
HalfBinSize=diff(CenterBins)/2 ;

% Vector of bins edges (not centers)
Bins=sort( [ CenterBins(1)-HalfBinSize(1) ; ...
             CenterBins(1:end-1)+HalfBinSize ; ...
             CenterBins(end)+HalfBinSize(end) ] ) ;   
         
% Calculate output
% for i=1:length(CenterBins)
%     ind = find(X>=Bins(i) & X<=Bins(i+1) );
%     if ~isempty(ind)
%         BinIndexes(ind) = i;
%     end
% end

[~,BinIndexes] = histc(X,Bins);


   


