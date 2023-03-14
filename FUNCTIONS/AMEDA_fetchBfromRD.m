
function [b] = AMEDA_fetchBfromRD(mat_Rd,x,y,Dx,mask)

%Load Rd file for global ocean
load(mat_Rd)
Rd_Lim = 12;
Rd_Chelton(Rd_Chelton<Rd_Lim/3) = Rd_Lim/3; 
Rd = interp2(lon_Rd,lat_Rd,Rd_Chelton,x,y,'*linear');
Rd(mask)=nan;

% Resolution parameters:
%----------------------------------------------
% mask for resolution parameters
IND = Rd<200 & ~isnan(Rd);

% gama is resolution coefficient which is the number of pixels per Rd.
% After test gama>2.4 is required to get the max number of eddies.
gama = nan(size(Rd));
gama(IND) = Rd(IND) ./ Dx(IND); % [0.1-1.5] for AVISO 1/8 (0.8 in average)
b = nanmax(1,round((1.2*gama)/2)); % always 1 for AVISO 1/8
b(isinf(b)) = 1;
