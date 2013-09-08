Lp   = 500;                   % unit: m, aqifer length
hLp  = 2;                     % unit: m, fixed head at the downgradient boundary
W    = 0.5/(365*24*60*60);    % unit: m/year -> m/s, surface recharge
K    = 1e-4;                  % unit: m/s, hydraulic conductivity
ne   = 0.35;                  % effective porosity
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
ylim([0 20])
xlabel('Distance (m)')
ylabel('Travel Time (years)')
title('Transit time from water table to exit boundary')
