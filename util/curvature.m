% Plot Gaussian curvature against C(theta)

theta = 0:0.001:pi;

r = 20/(2*pi);  % Minor radius
n = 4;        % Ratio n=R/r

a = r*sqrt(n^2 - 1);
eta = atanh(a/(n*r));

curvature = cos(theta)./((r^2)*(n+cos(theta)));

ttheta = acos(R/r - (a^2)./(r*(R+r*cos(theta))));
ttheta2 = acos(R/r - (a^2 * curvature)./(cos(theta)));

C = 10*(cosh(eta) - cos(ttheta2)).^2/(a^2);

figure('Name','Curvature','position', [0, 0, 300, 250]);
hold all
plot(theta, curvature);
plot([0 pi],[0 0])
ylabel('Gaussian curvature');
xlabel('\theta');
xlim([0 pi]);
ax = gca;
ax.XTick = [0 pi/2 pi];
ax.XTickLabels = {'0','\pi/2','\pi'};
legend('Torus','Flat')

figure('Name','Scaling','position', [10, 10, 300, 250]);
hold all
plot(theta, C);
plot([0 pi],[1 1])
ylabel('Diffusion scaling');
xlabel('\theta');
xlim([0 pi]);
ax = gca;
ax.XTick = [0 pi/2 pi];
ax.XTickLabels = {'0','\pi/2','\pi'};
legend('Torus','Flat')

figure('Name','Both','position', [20, 20, 300, 250]);
hold all
plot(curvature, C)
xlabel('Gaussian curvature')
ylabel('Diffusion scaling')
xlim([-0.1 0.033])
ax = gca;
ax.XTick = [-0.1 -0.05 0 0.03];
ax.XTickLabels = {'-0.1','-0.05','0','0.03'};

