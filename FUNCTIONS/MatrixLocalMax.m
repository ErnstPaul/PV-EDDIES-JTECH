function Locations = MatrixLocalMax(A,varargin)
%MATRIXLOCALMAX Find local maxima in a matrix by returning a
% new matrix that equal one at the maxima and zero elsewhere.
%
% Locations = MatrixLocalMax(A) finds all isolated local maxima
% of the matrix A.
%
% Locations = MatrixLocalMax(A,Reach) finds local maxima of the
% matrxi A that are separated by at least Reach points. Reach
% should be an integer larger than or equal to one. Reach is
% included to filter maxima in matrices with noise.
%
% Note: If two points within Reach distance of one another are
% identical, MatrixLocalMax only returns one value
%
% Example:
% A = peaks;
% Locations = MatrixLocalMax(A);
%
% See also MAX


%
% Author: David Eyre
% University of Utah
% Copyright 2000
%
    if nargin==1
        Reach = 1;
    elseif nargin==2
        Reach = max(round(varargin{1}),1);
    else
        error('Incorrect number of input arguments')
    end


        [m,n] = size(A);
    P = (1:m)+Reach;
    Q = (1:n)+Reach;
   
%
% Pad A so boundary conditions become transparent
%
    APad = (min(min(A))-1)*ones(m+(2*Reach),n+(2*Reach));
    APad(P,Q) = A;
   
%
% Line up arrays to determine local maxima in boxes centered
% around every point with radius = Reach. This uses a lot of
% memory but is more efficient that non-vector approaches
%
    Count = 0;
    Sheets = (2*Reach+1)^2;
    Work = zeros(m,n,Sheets);
        for j=-Reach:Reach
        for k=-Reach:Reach
            Count = Count+1;
            Work(:,:,Count) = APad(P+j,Q+k);
        end
    end
    [WorkMax,IMax] = max(Work,[],3);
   
%
% Maxima correspond to points mapped to themselves
%
    Invariant = ((2*Reach+1)^2+1)/2;
    Locations = 1-sign(abs(IMax-Invariant));