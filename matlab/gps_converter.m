nargin=length(VarName1);
zone=30;
northp=true;
lat=zeros(1,nargin); 
lon=zeros(1,nargin);
lat_d=zeros(1,nargin);
lon_d=zeros(1,nargin);
lat_m=zeros(1,nargin); 
lon_m=zeros(1,nargin);  
alt=zeros(1,nargin); 
time=zeros(1,nargin);
for j=1:nargin,
   
  lat_d(j)=floor(VarName4(j,1)/100);
  lon_d(j)=floor(VarName6(j,1)/100);
  
  lat_m(j)=(VarName4(j,1)-lat_d(j)*100)/60;
  lon_m(j)=(VarName6(j,1)-lon_d(j)*100)/60;
  
  lat(j)=lat_d(j)+lat_m(j);
  lon(j)=-lon_d(j)-lon_m(j);
  
  alt(j)=VarName8(j,1);
  
   time(j)=VarName3(j,1);
     
end

[x, y, gam, k] = utm_fwd(zone, northp, lat, lon);

% R_1=cov(x);
% 
% R_2=cov(y);
% 
% R_3=cov(alt);
% 
% R=[R_1 0 ;0 R_2 ];


