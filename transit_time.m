function transit_time(x,y,Ksi,ha,L,dKsi)
% calculate the transit time using equation (20) on R. Chesnaux et al.
% An Analytical Solution for Ground Water Transit Time through Unconfined
% Aquifers. Ground Water, 2005
% delta Ksi = 8e-8

v = dKsi:dKsi:max(max(Ksi));
figure
[C, ~] = contour(x,y,Ksi,v);

sind = zeros(length(v),1);
eind = zeros(length(v),1);

i = 1;
e = 0;  % ending point
while e < size(C, 2)
    s = e + 2;      % starting point
    sind(i) = s;
    e = e + C(2, s-1) + 1;
    eind(i) = e;
    i = i + 1;
end

nx   = 100;
icx1 = zeros(1, nx);
icy1 = zeros(1, nx);
dp   = zeros(1, nx);
ds   = zeros(1, nx);
t    = zeros(1, length(v));
xt   = zeros(1, length(v));
ne   = 0.35;
for n = 1 : length(v)-1
    cx1 = C(1,sind(n):eind(n));
    cy1 = C(2,sind(n):eind(n));
    hh = interp1(x,ha,cx1);
    cy1 = cy1.*hh/cy1(1);
    plot(cx1, cy1);
    hold on
    
    cx2 = C(1,sind(n+1):eind(n+1));
    cy2 = C(2,sind(n+1):eind(n+1));   
    hh = interp1(x,ha,cx2);
    cy2 = cy2.*hh/cy2(1);
    plot(cx2, cy2);
    ylim([0 10])
    drawnow
    
    icx2 = linspace(min(cx2),max(cx2),nx);
    icy2 = interp1(cx2,cy2,icx2);
    
    icx = linspace(min(cx1),max(cx1),10000);
    icy = interp1(cx1,cy1,icx);    
    for i = 1 : nx
        d = sqrt((icx2(i)-icx).^2+(icy2(i)-icy).^2);
        [dmin, ind] = min(d);
        icx1(i) = icx(ind);
        icy1(i) = icy(ind);
        dp(i)   = dmin;
    end
%     icx1 = linspace(min(cx1),max(cx1),nx);
%     icy1 = interp1(cx1,cy1,icx1);
%     dp = icy2 - icy1;
    xm = (icx1 + icx2) / 2;
    ym = (icy1 + icy2) / 2;
    ds(2:end) = sqrt((xm(2:end)-xm(1:end-1)).^2+(ym(2:end)-ym(1:end-1)).^2);
    s = cumsum(ds);
    t(n) = ne/dKsi*trapz(s, dp)/365/24/60/60;
    xt(n) = icx1(1);
end
xt(end) = L;
plot(x,ha,'b--','LineWidth',2)
ylim([0 10])
figure
plot(xt,t,'ko','MarkerFaceColor','k','MarkerSize',3); hold on