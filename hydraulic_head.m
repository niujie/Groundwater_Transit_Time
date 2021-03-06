L  = 500;        % unit: m, aquifer length
B  = 10;         % unit: m, aquifer depth
dx = 5;          % unit: m, x direction step size 
dy = 1;          % unit: m, y direction step size
ne = 0.35;       % effective porosity

x = (0:dx:L)';
y = (0:dy:B)';
[X,Y] = meshgrid(x,y);

W  = 0.5/(365*24*60*60);    % unit: m/year -> m/s, surface recharge
K  = 1e-4;                  % unit: m/s, hydraulic conductivity
KH = K;                     % unit: m/s, KXX, horizontal hydraulic head
KV = K;                % unit: m/s, KYY, vertical   hydraulic head
alpha = KH/dx^2;
beta  = KV/dy^2;
alphaKsi = 1/dx^2/KV;
betaKsi  = 1/dy^2/KH;

hL = 2;                     % unit: m, fixed head at the downgradient boundary
% unit: m, upgradient boundary fixed head, by 
% Dupuit-Forchheimer limits, hmax/Lp < 0.1
h0 = sqrt(W/K*L^2+hL^2);
ha       = sqrt(W/K*(L^2-x.^2)+hL^2);
h        = zeros(length(y),length(x))+B;
h(:,end) = hL;
h(end,:) = ha;
for i = 1 : size(x)
    Y(:,i) = Y(:,i) .* ha(i) / max(Y(:,i));
end

I = (2:size(h,1)-1)';
J = (2:size(h,2)-1)';

qx = zeros(size(h));
qy = zeros(size(h));

nplot = 10000;
omega = 1.5;
error = 1e6;
nstep = 0;
while error > 1e-8
    nstep = nstep + 1;
    h_old = h;

    % up boundary condition
    % this one gives the right velocity field
    % but head can't converge to analytical water table
    % and converges very slowly
    % h(end,:) = 2*dy*W/KV/3 + 4/3*h(end-1,:) - 1/3*h(end-2,:);
    % h(end,:) = dy*W/KV + h(end-1, :);  % first order
    
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
%     h = (1-omega)*h_old + omega*h;
    error = norm(h-h_old);
    if mod(nstep, nplot) == 0 || error - 1e-8 < eps
        [~,handle1] = contour(X,Y,h,20);
        hold on
        plot(x,ha,'b--','LineWidth',2)
        ylim([0 10])

        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end
end
% when doing the particle tracking to estimate traveling time
% velocities calculated from the derivative of head are correct
qx(:,1)   = -KH*(-3*h(:,1)+4*h(:,2)-h(:,3))/(2*dx);
qx(:,J)   = -KH*(h(:,J+1)-h(:,J-1))/(2*dx);
qx(:,end) = -KH*(3*h(:,end)-4*h(:,end-1)+h(:,end-2))/(2*dx);
qy(1,:) = -KV*(-3*h(1,:)+4*h(2,:)-h(3,:))/(2*dy);
qy(I,:) = -KV*(h(I+1,:)-h(I-1,:))/(2*dy);
qy(end,:) = -KV*(3*h(end,:)-4*h(end-1,:)+h(end-2,:))/(2*dy);
vx = qx/ne*86400;
vy = qy/ne*86400;

hold on;
handle2 = quiver(X,Y,vx,vy,3.0);
adjust_quiver_arrowhead_size(handle2, 0.1);

Ksi        = zeros(size(h));
% top boundary, use this condition doesn't mathe the velocities obtained
% from qx = -K dh/dx and qy = -K dh/dy
% but the transit time result matches
% Ksi(:,end) = qx(:,end) .* Y(:,end);
Ksi(end,:) = W*x;

error = 1e6;
nstep = 0;
nplot = 100;
while error > 1e-13
    nstep = nstep + 1;
    Ksi_old = Ksi;
    % Gaussian iteration
    for j = 2 : size(h,2)-1
        for i = 2 : size(h,1)-1
            Ksi(i,j) = (alphaKsi*(Ksi(i,j+1)+Ksi(i,j-1)) + betaKsi*(Ksi(i+1,j)+Ksi(i-1,j))) / ...
                (2*(alphaKsi+betaKsi));
        end
    end
    % top boundary, use this condition mathes the velocities obtained from
    % qx = -K dh/dx and qy = -K dh/dy
    % but the transit time result doesn't match
%     for j = 2 : size(h,2)-1
%         Ksi(end,j) = KH*dy/3/dx*(h(end,j-1)-h(end,j+1)) + 4/3*Ksi(end-1,j) - Ksi(end-2,j)/3;
%     end
   Ksi(:,end) = 4/3*Ksi(:,end-1) - 1/3*Ksi(:,end-2);   % right
    % SOR iteration
    Ksi = (1-omega)*Ksi_old + omega*Ksi;
    error = norm(Ksi-Ksi_old);   
    if mod(nstep, nplot) == 0 || error - 1e-13 < eps
        [~,handle1] = contour(X,Y,h,20);
        hold on
        plot(x,ha,'b--','LineWidth',2)
        [C,handle2] = contour(X,Y,Ksi,20);

        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        ylim([0 10])
        drawnow
        hold off
    end    
end
% when doing the particle tracking to estimate traveling time
% velocities calculated from the derivative of stream function are incorrect
% use the derivative of head above instead
% qx(1,:) = (-3*Ksi(1,:)+4*Ksi(2,:)-Ksi(3,:))/(2*dy);
% qx(I,:) = (Ksi(I+1,:)-Ksi(I-1,:))/(2*dy);
% qx(end,:) = (3*Ksi(end,:)-4*Ksi(end-1,:)+Ksi(end-2,:))/(2*dy);
% qy(:,1) = -(-3*Ksi(:,1)+4*Ksi(:,2)-Ksi(:,3))/(2*dx);
% qy(:,J) = -(Ksi(:,J+1)-Ksi(:,J-1))/(2*dx);
% qy(:,end) = -(3*Ksi(:,end)-4*Ksi(:,end-1)+Ksi(:,end-2))/(2*dx);
% vx = qx/ne*86400;
% vy = qy/ne*86400;
hold on;
plot(x,ha,'b--','LineWidth',2)
handle3 = quiver(X,Y,vx,vy,3.0);
adjust_quiver_arrowhead_size(handle3, 0.1);
        
set(handle1,'ShowText','on')
set(handle2,'ShowText','on')
ylim([0 10])

transit_time(x,y,Ksi,ha,L,5e-7);
transit_time_analytical(L, hL, W, K, ne);
transit_time(x,y,Ksi,ha,L,8e-8);
transit_time_analytical(L, hL, W, K, ne);
particle_tracking;