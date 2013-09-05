L  = 500;        % unit: m, aquifer length
D  = 10;         % unit: m, aquifer depth
dx = 5;          % unit: m, x direction step size 
dy = 1;          % unit: m, y direction step size

KH = 1.0e-4;    % unit: m/s, KXX, horizontal hydraulic head
KV = KH/1000;   % unit: m/s, KYY, vertical   hydraulic head

alpha = 1/dx^2;
beta  = 1/dy^2;
alphaKsi = 1/dx^2;
betaKsi  = 1/dy^2;

x = (0:dx:L)';
y = (0:dy:D)';

W    = 0.5/(365*24*60*60);    % unit: m/year -> m/s, surface recharge
K    = 1e-4;                  % unit: m/s, hydraulic conductivity
hL = 2;                     % unit: m, fixed head at the downgradient boundary
% unit: m, upgradient boundary fixed head, by 
% Dupuit-Forchheimer limits, hmax/Lp < 0.1
h0 = sqrt(W/K*L^2+hL^2);
ha       = sqrt(W/K*(L^2-x.^2)+hL^2);

h        = zeros(length(y),length(x));
h(1,:)   = ha;
h(:,end) = hL;

Ksi      = zeros(length(y),length(x));
Ksi(1,:) = W*x;

I = (2:size(h,1)-1)';
J = (2:size(h,2)-1)';

error = 1e6;
nstep = 0;
nplot = 100;
while error > 1e-3
    nstep = nstep + 1;
    h_old = h;
    for j = 2 : size(h,2)-1
        for i = 2 : size(h,1)-1
            h(i,j) = (alpha*(h(i,j+1)+h(i,j-1)) + beta*(h(i+1,j)+h(i-1,j))) / ...
                (2*(alpha+beta));
        end
    end
%     h(I,J) = (alpha*(h(I,J+1)+h(I,J-1)) + beta*(h(I+1,J)+h(I-1,J))) / ...
%         (2*(alpha+beta));
    h(2:end,1) = 4/3*h(2:end,2)-1/3*h(2:end,3);
    h(end,1:end-1) = 4/3*h(end-1,1:end-1)-1/3*h(end-2,1:end-1);
    error = norm(h-h_old);
    if mod(nstep, nplot) == 0 || error - 1e-3 < eps
        [~,handle] = contour(x,y,flipud(h),20);
        set(handle,'ShowText','on')
        %pcolor(h)
        %shading interp
        %colorbar
        hold on
        plot(x,ha,'k','LineWidth',3)
        %ylim([0 10])
        %plot(x,h(:,:)); hold on
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||h||_{2} = ', num2str(error)])
        drawnow
        hold off
    end
end

qx = zeros(size(h));
qy = zeros(size(h));
qx(:,1)   = -K*(h(:,2)-h(:,1))/dx;
qx(:,J)   = -K*(h(:,J+1)-h(:,J))/dx;
qx(:,end) = -K*(h(:,end)-h(:,end-1))/dx;
qy(I,:) = -K*(h(I,:)-h(I-1,:))/dy;
ne   = 0.35;                  % effective porosity
vx = qx/ne;
vy = qy/ne;