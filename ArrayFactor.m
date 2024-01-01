% Date : 12.06.2021
% axis settings
% theta = 0:0.01:2*pi;
theta = pi/2; %theta kesit 
phi= 0:0.01:2*pi;
%phi=pi/2; %phi kesit
axis=phi; 
%--------------------------------------- 
f=3e9;
c=3e8;
lambda = c/f;
L=lambda/2;

theta0 = 90; % yönlenme açıları
phi0 = 7; 
A=1; %genlik
N=9; % eleman sayısı
d=lambda/2; % elemanlar arası mesafe
%---------------------------------------
ele_pattern=((cos(cos(theta)*pi*L/lambda)-cos(pi*L/lambda))./(sin(theta))).^2;% normalize ışıma siddeti
ele_pattern= ones(size(phi)).*F; %phi de sürekli hale getirdik

AF=A;
for i=1:floor(N/2)
   dn= i*d;
   B=-(sin(theta0*pi/180)*sin(phi0*pi/180)*dn*2*pi/lambda);
   AF=AF+ 2*A*cos((sin(theta).*sin(phi)*dn*2*pi/lambda)+B);  
end

figure,
polarplot(axis,ele_pattern)
ax=gca;
ax.ThetaZeroLocation= 'top';
title('Element pattern(phi)(theta=90)')

figure,
polarplot(axis,AF)
ax=gca;
ax.ThetaZeroLocation= 'top';
title('Array factor(theta=90)')

figure,
polarplot(axis,AF.*ele_pattern)
ax=gca;
ax.ThetaZeroLocation= 'top';
title('Total(theta=90)')

