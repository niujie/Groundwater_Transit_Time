L  = 500;        % unit: m, aquifer length
B  = 10;         % unit: m, aquifer depth
dx = 5;          % unit: m, x direction step size 
dy = 1;          % unit: m, y direction step size

%KH = 1.0e-4;    % unit: m/s, KXX, horizontal hydraulic head
%KV = KH/1000;   % unit: m/s, KYY, vertical   hydraulic head

alpha = 1/dx^2;
beta  = 1/dy^2;

x = (0:dx:L)';
y = (0:dy:B)';
[X,Y]=meshgrid(x,y);

W    = 0.5/(365*24*60*60);    % unit: m/year -> m/s, surface recharge
K    = 1e-4;                  % unit: m/s, hydraulic conductivity
hL   = 2;                     % unit: m, fixed head at the downgradient boundary
% unit: m, upgradient boundary fixed head, by 
% Dupuit-Forchheimer limits, hmax/Lp < 0.1
h0 = sqrt(W/K*L^2+hL^2);
ha       = sqrt(W/K*(L^2-x.^2)+hL^2);
h        = zeros(length(y),length(x));
h(1,:)   = ha;
h(:,end) = hL;

I = (2:size(h,1)-1)';
J = (2:size(h,2)-1)';

xp         = (0:length(x)-1)'*dx^2/sqrt(dx^2+dy^2);
qx         = zeros(size(h));
qx(:,end)  = -K*(3*h(:,end)-4*h(:,end-1)+h(:,end-2))/(2*dx);
Ksi        = zeros(size(h));
Ksi(1,1)   = -(-3*h(1,1)+4*h(1,2)-h(1,3))/(2*dx);
Ksi(1,J)   = -(h(1,J+1)-h(1,J-1))/(2*dx);
Ksi(1,end)   = -(3*h(1,end)-4*h(1,end-1)+h(1,end-2))/(2*dx);
%Ksi(1,:)   = W*xp;

error = 1e6;
nstep = 0;
nplot = 1;
while error > 1e-3
    nstep = nstep + 1;
    h_old = h;
    Ksi_old = Ksi;
    for j = 2 : size(h,2)-1
        for i = 2 : size(h,1)-1
            h(i,j) = (alpha*(h(i,j+1)+h(i,j-1)) + beta*(h(i+1,j)+h(i-1,j))) / ...
                (2*(alpha+beta));
            Ksi(i,j) = (alpha*(Ksi(i,j+1)+Ksi(i,j-1)) + beta*(Ksi(i+1,j)+Ksi(i-1,j))) / ...
                (2*(alpha+beta));            
        end
    end
%     h(I,J) = (alpha*(h(I,J+1)+h(I,J-1)) + beta*(h(I+1,J)+h(I-1,J))) / ...
%         (2*(alpha+beta));
    h(2:end,1) = 4/3*h(2:end,2)-1/3*h(2:end,3); % left
    h(end,1:end-1) = 4/3*h(end-1,1:end-1)-1/3*h(end-2,1:end-1); % bottom
    qx(:,end)  = -K*(3*h(:,end)-4*h(:,end-1)+h(:,end-2))/(2*dx);
    Ksi(:,end) = Ksi(:,end-1) + qx(:,end);
    error = norm(h-h_old)+norm(Ksi-Ksi_old);
    if mod(nstep, nplot) == 0 || error - 1e-3 < eps
        [~,handle1] = contour(x,y,h,20);
        %set(handle1,'ShowText','on')
        hold on
        [~,handle2] = contour(x,y,Ksi,50);
        set(gca,'YDir','reverse');
        %plot(x,ha,'k','LineWidth',3)
        %ylim([0 10])
        qx = zeros(size(h));
        qy = zeros(size(h));
        qx(:,1)   = -K*(h(:,2)-h(:,1))/dx;
        qx(:,J)   = -K*(h(:,J+1)-h(:,J))/dx;
        qx(:,end) = -K*(h(:,end)-h(:,end-1))/dx;
        qy(1,:) = -K*(h(2,:)-h(1,:))/dy;
        qy(I,:) = -K*(h(I,:)-h(I-1,:))/dy;
        qy(end,:) = -K*(h(end,:)-h(end-1,:))/dy;
        ne   = 0.35;                  % effective porosity
%         vx = qx/ne;
%         vy = qy/ne;
        vx = zeros(size(h));
        vy = zeros(size(h));
        vx(1,:) = -(Ksi(2,:)-Ksi(1,:))/dy;
        vx(I,:) = -(Ksi(I+1,:)-Ksi(I-1,:))/(2*dy);
        vx(end,:) = -(Ksi(end,:)-Ksi(end-1,:))/dy;
        vy(:,1) = (Ksi(:,2)-Ksi(:,1))/dx;
        vy(:,J) = (Ksi(:,J+1)-Ksi(:,J-1))/(2*dx);
        vy(:,end) = (Ksi(:,end)-Ksi(:,end-1))/dx;
        hold on;
        handle3 = quiver(X,Y,vx,vy);
        adjust_quiver_arrowhead_size(handle3, 0.1);
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end
end