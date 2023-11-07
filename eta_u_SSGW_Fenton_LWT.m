clear

addpath Fenton/ 
addpath SSGW/


H = 0.8 ; % m
h = 2.0 ; % m
T = 10.0 ; % s
g = 9.806 ; % m/s^2
N = 2^12 ;%1024 ; %1024 ; 
tol = 1e-10 ;

%% 
omg = 2 * pi / T ;
[kh] = dispersionLZ(h, T) ; 
k = kh / h ; 
L = 2*pi / k  ;
kH2 = k * H / 2 ; 
kH = k * H ; 

%% LWT
theta_LWT = linspace (0, 2*pi, 100) ; 
eta_LWT = H /2 * cos (theta_LWT) ;
zele_LWT = linspace (-h, max(eta_LWT)) ;
u_max_LWT = ( omg * H / 2 / sinh(kh) ) .* cosh (k * (zele_LWT + h)) ; 

%% SSGW
if 1.0-tanh(kh) < tol
    scaleVel = sqrt(g/k);
    scaleLength = 1/k;
    fprintf('deep water. scale is 1/k.\n')
else
    scaleVel = sqrt(g*h);
    scaleLength = h;
    fprintf('finite water depth. scale is h.\n')
end
 
% get x and y for free surface 
[zs, ws, PP, ~] = SSGW_computeW(kh,kH2,N,tol, []) ;

% get w at free surface
xs = real(zs) ;
ys = imag(zs) ; 
u_surface = scaleVel* ( real(ws) + PP(4) ) ; 
dal = pi / PP(2) / N  ;

% u vertical profile 
ZZ = [xs(1)+1i*linspace(-h, (ys(1)-dal)*scaleLength, 10)/scaleLength] ;  % umax 
[~,~,~, W] = SSGW_computeW(kh,kH2,N,tol, ZZ) ;


%% Fenton's
tic
% % [eta, uvel, theta, zele, uvel_surface] = StreamFunction_surface_u (H, h, T, plottype, phase_theta, zlocation)
% [~, uvel, theta, zele, uvel_fenton_surface] = StreamFunction_surface_u (H, h, T, 'temporal', 0, 0) ; 
% toc

% plot(theta, eta, ':r', 'linewidth', 2)
[~, uvel, ~, zele] = StreamFunction_surface_u (H, h, T, 'vertical', 0, 0) ; 
[eta, ~, theta, ~] = StreamFunction_surface_u (H, h, T, 'temporal', 0, -h) ; 


figure; hold on; box on
plot(u_max_LWT, zele_LWT, '-k', 'LineWidth', 3)
plot(scaleVel*real(W), scaleLength*imag(ZZ), '-b', 'LineWidth', 3)
plot(uvel, zele, '-r', 'LineWidth', 3)
legend('LWT', 'SSGW', 'Fenton', 'orientation', 'vertical')
xlim([0, 2])
ylim([-inf, 1.2])
xlabel('u_{max} (m/s)')
ylabel('z (m)')
set(gcf, 'Position', [ 2409         164         338         474])
set (gca, 'fontsize', 22)
 
figure; hold on; box on
plot(theta_LWT, eta_LWT, '-k', 'LineWidth', 3)
plot(k*scaleLength*xs, scaleLength*ys, '-b', 'LineWidth', 3)
plot(theta, eta, '-r', 'LineWidth', 3)
legend('SSGW', 'Fenton', 'LWT', 'orientation', 'vertical')
xlim([0, 2*pi])
xlabel('kx - \omegat')
ylabel('\eta (m)')
set (gca, 'fontsize', 22)

rmpath Fenton/ 
rmpath SSGW/
