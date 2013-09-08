function transit_time_analytical(L, hL, W, K, ne)
% analytical solution of groundwater transit time using equation (17) on
% R. Chesnaux et al.
% An Analytical Solution for Ground Water Transit Time through Unconfined
% Aquifers. Ground Water, 2005

Lp   = L;                   % unit: m, aqifer length
hLp  = hL;                     % unit: m, fixed head at the downgradient boundary
hmax = 0.05*Lp;               % unit: m, upgradient boundary fixed head, by
                              % Dupuit-Forchheimer limits, hmax/Lp < 0.1
alpha = Lp^2 + K*hLp^2/W;
dx = 5;
x = (0:dx:L)';
t1 = ne*sqrt(alpha/(K*W))*(Lp*sqrt(1/Lp^2-1/alpha)-x.*sqrt(1./x.^2-1/alpha) + ...
    log((sqrt(alpha)./x+sqrt(alpha./x.^2-1))/(sqrt(alpha)/Lp+sqrt(alpha/Lp^2-1))));
t2 = hLp*ne/W*log(Lp./x);
plot(x,t1/(60*60*24*365),'b-')
hold on
plot(x,t2/(60*60*24*365),'r--')
xlabel('Distance (m)')
ylabel('Travel Time (years)')
title('Transit time from water table to exit boundary')
