L  = 500;        % unit: m, aquifer length
B  = 10;         % unit: m, aquifer depth
dx = 5;          % unit: m, x direction step size 
dy = 1;        % unit: m, y direction step size
ne = 0.35;                  % effective porosity

x = (0:dx:L)';
y = (0:dy:B)';
[X,Y] = meshgrid(x,y);

W  = 0.5/(365*24*60*60);    % unit: m/year -> m/s, surface recharge
K  = 1e-4;                  % unit: m/s, hydraulic conductivity
KH = K;                     % unit: m/s, KXX, horizontal hydraulic head
KV = K/1000;                % unit: m/s, KYY, vertical   hydraulic head
alpha = KH/dx^2;
beta  = KV/dy^2;
alphaKsi = 1/dx^2/KH;
betaKsi  = 1/dy^2/KV;

hL = 2;                     % unit: m, fixed head at the downgradient boundary
% unit: m, upgradient boundary fixed head, by 
% Dupuit-Forchheimer limits, hmax/Lp < 0.1
h0 = sqrt(W/K*L^2+hL^2);
ha       = sqrt(W/K*(L^2-x.^2)+hL^2);
h        = zeros(length(y),length(x))+B;
h(:,end) = hL;
h(end,:) = ha;

I = (2:size(h,1)-1)';
J = (2:size(h,2)-1)';

qx = zeros(size(h));
qy = zeros(size(h));
vx = zeros(size(h));
vy = zeros(size(h));

nplot = 10000;
omega = 1.5;
error = 1e6;
nstep = 0;
while error > 1e-8
    nstep = nstep + 1;
    h_old = h;

    %h(end,:) = 2*dy*W/KV/3 + 4/3*h(end-1,:) - 1/3*h(end-2,:); % up
    h(2:end-1,1) = 4/3*h(2:end-1,2)-1/3*h(2:end-1,3); % left
    h(1,1:end-1) = 4/3*h(2,1:end-1)-1/3*h(3,1:end-1); % bottom    
    % Gaussian iteration
%     for j = 2 : size(h,2)-1
%         for i = 2 : size(h,1) - 1
%             h(i,j) = (alpha*(h(i,j+1)+h(i,j-1)) + beta*(h(i+1,j)+h(i-1,j))) / ...
%                 (2*(alpha+beta));           
%         end
%     end
    % vector form is just Jocoby iteration
    h(I,J) = (alpha*(h(I,J+1)+h(I,J-1)) + beta*(h(I+1,J)+h(I-1,J))) / ...
        (2*(alpha+beta));
    % SOR iteration
    % h = (1-omega)*h_old + omega*h;
    error = norm(h-h_old);
    if mod(nstep, nplot) == 0 || error - 1e-8 < eps
        [~,handle1] = contour(x,y,h,20);
        hold on
%         plot(x,ha,'k','LineWidth',3)
%         ylim([0 10])
%         qx(:,1)   = -K*(-3*h(:,1)+4*h(:,2)-h(:,3))/(2*dx);
%         qx(:,J)   = -K*(h(:,J+1)-h(:,J-1))/(2*dx);
%         qx(:,end) = -K*(3*h(:,end)-4*h(:,end-1)+h(:,end-2))/(2*dx);
%         qy(1,:) = -K*(-3*h(1,:)+4*h(2,:)-h(3,:))/(2*dy);
%         qy(I,:) = -K*(h(I+1,:)-h(I-1,:))/(2*dy);
%         qy(end,:) = -K*(3*h(end,:)-4*h(end-1,:)+h(end-2,:))/(2*dy);
%         vx = qx/ne;
%         vy = qy/ne;
%         vx = vx./max(vx(:))*86400;
%         vy = vy./max(vy(:))*86400;
% 
%         hold on;
%         handle2 = quiver(X,Y,vx,vy);
%         %adjust_quiver_arrowhead_size(handle2, 0.1);
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end
end

Ksi        = zeros(size(h));
Ksi(end,:) = W*x;

error = 1e6;
nstep = 0;
nplot = 100;
while error > 1e-12
    nstep = nstep + 1;
    Ksi_old = Ksi;
    % Gaussian iteration
    for j = 2 : size(h,2)-1
        for i = 2 : size(h,1)-1
            Ksi(i,j) = (alphaKsi*(Ksi(i,j+1)+Ksi(i,j-1)) + betaKsi*(Ksi(i+1,j)+Ksi(i-1,j))) / ...
                (2*(alphaKsi+betaKsi));
        end
    end
    % top boundary
%     for j = 2 : size(h,2)-1
%         Ksi(end,j) = dy/3/dx*(h(end,j-1)-h(end,j+1)) - h(end,j) + 4/3*h(end-1,j) - h(end-2,j)/3 + ...
%             4/3*Ksi(end-1,j) - Ksi(end-2,j)/3 - 2*dy*W/3/KV;
%     end    
    Ksi(:,end) = 4/3*Ksi(:,end-1) - 1/3*Ksi(:,end-2);   % right
    % SOR iteration
    Ksi = (1-omega)*Ksi_old + omega*Ksi;
    error = norm(Ksi-Ksi_old);   
    if mod(nstep, nplot) == 0 || error - 1e-12 < eps
        [~,handle1] = contour(x,y,h,20);
        hold on
        [C,handle2] = contour(x,y,Ksi,20);
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

transit_time(x,y,Ksi,ha,L,5e-7);
transit_time(x,y,Ksi,ha,L,8e-8)
transit_time_analytical;