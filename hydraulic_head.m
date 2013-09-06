L  = 500;        % unit: m, aquifer length
B  = 10;         % unit: m, aquifer depth
dx = 5;          % unit: m, x direction step size 
dy = 1;          % unit: m, y direction step size

KH = 1.0e-4;    % unit: m/s, KXX, horizontal hydraulic head
KV = KH/1000;   % unit: m/s, KYY, vertical   hydraulic head

alpha = KH/dx^2;
beta  = KV/dy^2;
alphaKsi = 1/dx^2/KH;
betaKsi  = 1/dy^2/KV;

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
h(1,:)   = ha(1);
index    = find(y>=B-hL);
h(index,end) = hL;

I = (2:size(h,1)-1)';
J = (2:size(h,2)-1)';

qx = zeros(size(h));
qy = zeros(size(h));

error = 1e6;
nstep = 0;
nplot = 100;
omega = 1.0;
while error > 1e-6
    nstep = nstep + 1;
    h_old = h;
    
    index = find(y<B-hL);
    h(index,end) = 4/3*h(index,end-1)-1/3*h(index,end-2);  % right
    h(2:end,1) = 4/3*h(2:end,2)-1/3*h(2:end,3); % left
    h(end,1:end-1) = 4/3*h(end-1,1:end-1)-1/3*h(end-2,1:end-1); % bottom    
    % Gaussian iteration
    for j = 2 : size(h,2)-1
        for i = 2 : size(h,1)-1
            h(i,j) = (alpha*(h(i,j+1)+h(i,j-1)) + beta*(h(i+1,j)+h(i-1,j))) / ...
                (2*(alpha+beta));           
        end
    end
    % vector form is just Jocoby iteration
%     h(I,J) = (alpha*(h(I,J+1)+h(I,J-1)) + beta*(h(I+1,J)+h(I-1,J))) / ...
%         (2*(alpha+beta));
    % SOR iteration
    h = (1-omega)*h_old + omega*h;
    error = norm(h-h_old);
    if mod(nstep, nplot) == 0 || error - 1e-6 < eps
        [~,handle1] = contour(x,y,h,20);
        hold on
        set(gca,'YDir','reverse');
        %plot(x,ha,'k','LineWidth',3)
        %ylim([0 10])
        qx(:,1)   = -K*(-3*h(:,1)+4*h(:,2)-h(:,3))/(2*dx);
        qx(:,J)   = -K*(h(:,J+1)-h(:,J-1))/(2*dx);
        qy(1,:) = -K*(-3*h(1,:)+4*h(2,:)-h(3,:))/(2*dy);
        qy(I,:) = -K*(h(I+1,:)-h(I-1,:))/(2*dy);
        qy(end,:) = -K*(3*h(end,:)-4*h(end-1,:)+h(end-2,:))/(2*dy);
        ne   = 0.35;                  % effective porosity
        vx = qx/ne;
        vy = qy/ne;
        vx = vx./max(vx(:))*86400;
        vy = vy./max(vy(:))*86400;

        hold on;
        handle2 = quiver(X,Y,vx,vy);
        %adjust_quiver_arrowhead_size(handle2, 0.1);
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end
end

Ksi        = zeros(size(h));
Ksi(1,:)   = W*x;
qx(:,end)  = -K*(3*h(:,end)-4*h(:,end-1)+h(:,end-2))/(2*dx);
Ksi(:,end) = qx(:,end) .* y(end:-1:1);

error = 1e6;
nstep = 0;
nplot = 1;
omega = 1.5;
while error > 1e-6
    nstep = nstep + 1;
    Ksi_old = Ksi;
    % Gaussian iteration
    for j = 2 : size(h,2)-1
        for i = 2 : size(h,1)-1
            Ksi(i,j) = (alphaKsi*(Ksi(i,j+1)+Ksi(i,j-1)) + betaKsi*(Ksi(i+1,j)+Ksi(i-1,j))) / ...
                (2*(alphaKsi+betaKsi));
        end
    end
    % SOR iteration
    Ksi = (1-omega)*Ksi_old + omega*Ksi;
    error = norm(Ksi-Ksi_old);   
    if mod(nstep, nplot) == 0 || error - 1e-6 < eps
        [~,handle1] = contour(x,y,h,20);
        hold on
        [~,handle2] = contour(x,y,Ksi,20);
        set(gca,'YDir','reverse');
        vx(1,:) = -(-3*Ksi(1,:)+4*Ksi(2,:)-Ksi(3,:))/(2*dy);
        vx(I,:) = -(Ksi(I+1,:)-Ksi(I-1,:))/(2*dy);
        vx(end,:) = -(3*Ksi(end,:)-4*Ksi(end-1,:)+Ksi(end-2,:))/(2*dy);
        vy(:,1) = (-3*Ksi(:,1)+4*Ksi(:,2)-Ksi(:,3))/(2*dx);
        vy(:,J) = (Ksi(:,J+1)-Ksi(:,J-1))/(2*dx);
        vy(:,end) = (3*Ksi(:,end)-4*Ksi(:,end-1)+Ksi(:,end-2))/(2*dx);
        vx = vx./max(vx(:))*86400;
        vy = vy./max(vy(:))*86400;
        hold on;
        handle3 = quiver(X,Y,vx,vy);
        adjust_quiver_arrowhead_size(handle3, 0.1);
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end    
end
set(handle1,'ShowText','on')
set(handle2,'ShowText','on')