% particle tracking to evaluate the traveling time from water table to the
% exit boundary in groundwater, using 4th order Runge Kutta scheme

close all

% vx = vx / 86400;
% vy = vy / 86400;
[X,Y] = meshgrid(x,y);
% XI and YI are skewed mesh, just for plotting
XI = X;
YI = Y;
for i = 1 : size(x)
    YI(:,i) = Y(:,i) .* ha(i) / max(Y(:,i));
end
dt = 1;
nplot = 100;

T = zeros(size(x));
for i = 1 : length(x)
    xp = x(i);
    yp = B;

    nstep = 0;
    xs = [];
    ys = [];
    while xp > 0 && xp < L && yp > 0
        nstep = nstep + 1;

        % interpolate the water table position at xp
        hh = interp1(x, ha, xp);
        % re-scale the yp position based on the skewed water table position
        % for plotting
        ypi = yp.*hh/B;    
        xs = [xs, xp];
        ys = [ys, ypi];
        if mod(nstep, nplot) == 0
            plot(x,ha,'b--','LineWidth',2)
            hold on;    
            handle3 = quiver(XI,YI,vx,vy,3.0);
            adjust_quiver_arrowhead_size(handle3, 0.1);
            % note ypi here
            plot(xp, ypi, 'ro', 'MarkerSize', 3, 'MarkerFaceColor', 'r')
            plot(xs, ys, 'k--')
            hold off
            title(['number of steps: ', num2str(nstep), ' total time: ', ...
                num2str(nstep*dt)])
            xlim([0 L])
            ylim([0 B])
            drawnow
        end
        [xp,yp] = RK4_Particle_Tracking_2D(xp,yp,X,Y,vx,vy,dt);
    end
    
    T(i) = dt * nstep;
end