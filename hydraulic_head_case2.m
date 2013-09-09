% example from E.O. Friend et al.
% The Dual Formulation of Flow for Contaminant Transort Modeling 1.
% WRR 1985
L  = 1000;       % unit: m, aquifer length
B  = 60;         % unit: m, aquifer depth
dx = 50;         % unit: m, x direction step size 
dy = 10;         % unit: m, y direction step size
ne = 0.35;     % effective porosity

x = (0:dx:L)';
y = (0:dy:B)';
[X,Y] = meshgrid(x,y);

K  = 1e-3;                  % unit: m/s, hydraulic conductivity
alpha = K/dx^2;
beta  = K/dy^2;
alphaKsi = 1/dx^2/K;
betaKsi  = 1/dy^2/K;

H = 4e-4 * L;
h        = zeros(length(y),length(x))+B;
h(end,:) = B + H/2*(1+cos(pi*x/L));

I = (2:size(h,1)-1)';
J = (2:size(h,2)-1)';

qx = zeros(size(h));
qy = zeros(size(h));
vx = zeros(size(h));
vy = zeros(size(h));

nplot = 100;
omega = 1.5;
error = 1e6;
nstep = 0;
while error > 1e-8
    nstep = nstep + 1;
    h_old = h;

    h(2:end-1,1) = 4/3*h(2:end-1,2)-1/3*h(2:end-1,3); % left
    h(2:end-1,end) = 4/3*h(2:end-1,end-1) - 1/3*h(2:end-1,end-2);   % right
    h(1,:) = 4/3*h(2,:)-1/3*h(3,:); % bottom    
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
        qx(:,1)   = -K*(-3*h(:,1)+4*h(:,2)-h(:,3))/(2*dx);
        qx(:,J)   = -K*(h(:,J+1)-h(:,J-1))/(2*dx);
        qx(:,end) = -K*(3*h(:,end)-4*h(:,end-1)+h(:,end-2))/(2*dx);
        qy(1,:) = -K*(-3*h(1,:)+4*h(2,:)-h(3,:))/(2*dy);
        qy(I,:) = -K*(h(I+1,:)-h(I-1,:))/(2*dy);
        qy(end,:) = -K*(3*h(end,:)-4*h(end-1,:)+h(end-2,:))/(2*dy);
        vx = qx/ne;
        vy = qy/ne;

        hold on;
        handle2 = quiver(X,Y,vx,vy);
        adjust_quiver_arrowhead_size(handle2, 0.3);        
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end
end

Ksi        = zeros(size(h));

error = 1e6;
nstep = 0;
nplot = 10;
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
    for j = 2 : size(h,2)-1
        Ksi(end,j) = K*dy/3/dx*(h(end,j-1)-h(end,j+1)) + 4/3*Ksi(end-1,j) - Ksi(end-2,j)/3;
    end
    % SOR iteration
    Ksi = (1-omega)*Ksi_old + omega*Ksi;
    error = norm(Ksi-Ksi_old);   
    if mod(nstep, nplot) == 0 || error - 1e-12 < eps
        [~,handle1] = contour(x,y,h,10);
        hold on
        [~,handle2] = contour(x,y,Ksi,20);
        qx(1,:) = (-3*Ksi(1,:)+4*Ksi(2,:)-Ksi(3,:))/(2*dy);
        qx(I,:) = (Ksi(I+1,:)-Ksi(I-1,:))/(2*dy);
        qx(end,:) = (3*Ksi(end,:)-4*Ksi(end-1,:)+Ksi(end-2,:))/(2*dy);
        qy(:,1) = -(-3*Ksi(:,1)+4*Ksi(:,2)-Ksi(:,3))/(2*dx);
        qy(:,J) = -(Ksi(:,J+1)-Ksi(:,J-1))/(2*dx);
        qy(:,end) = -(3*Ksi(:,end)-4*Ksi(:,end-1)+Ksi(:,end-2))/(2*dx);
        
        vx = qx/ne;
        vy = qy/ne;        
        hold on;
        handle3 = quiver(X,Y,vx,vy);
        adjust_quiver_arrowhead_size(handle3, 0.3);
        title(['Iteration steps: ', num2str(nstep), ...
            ',    ||error||_{2} = ', num2str(error)])
        drawnow
        hold off
    end    
end
set(handle1,'ShowText','on')
set(handle2,'ShowText','on')