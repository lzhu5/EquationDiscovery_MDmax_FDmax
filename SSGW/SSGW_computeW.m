function [zs,ws,PP, W] = SSGW_computeW(kd,kH2,N,tol, Z)
%% created by Clamond. 
%% modified by Ling Zhu, July 12, 2023

% SSGW: Steady Surface Gravity Waves.
%       Computation of irrotational 2D periodic surface pure gravity waves 
%       of arbitrary length in arbitrary depth. 
%
% MANDATORY INPUT PARAMETERS:
% kd  = k*d   : relative depth (wavenumber "k" times mean water depth "d").
% kH2 = k*H/2 : steepness (half the total wave height "H" times the wavenumber "k").
%
% OPTIONAL INPUT PARAMETERS:
% N   : number of positive Fourier modes (default, N=2048).
% tol : tolerance (default, tol=1e-14).
% Z  : Complex abscissa where fields are desired inside the fluid (default Z = []).
%      Z should be strictly below the surface, i.e., -1 <= imag(Z) < eta(real(Z))
%      y = eta(x) being the equation of the free surface.

%
% OUTPUT PARAMETERS:
% zs  = complex abscissas at the free surface (at the computational nodes).
% ws  = complex velocity at the free surface (at the computational nodes) in 
%       the frame of reference moving with the wave.
% PP  = Physical Parameters: PP(1)=depth, PP(2)=wavenumber, PP(3)=waveheight, 
%       PP(4)=celerity c_e, PP(5)=celerity c_s, PP(6)=Bernoulli constant, 
%       PP(7)=crest height, PP(8)=trough height, PP(9)=impulse, 
%       PP(10)=potential energy, PP(11)=kinetic energy, PP(12)=radiation stress,
%       PP(13)=momentum flux, PP(14)=energy flux, PP(15)=group velocity.
% W    : Complex velocity in the bulk at abscissas Z.
%
% NOTE 1: This program computes waves of arbitrary length for all heights 
% up to about 99% of the maximum one. It is not designed to compute (near) 
% limiting waves. 
%
% NOTE 2: The output quantities are dimensionless with the following scaling:
%   In deep water:   rho = g = k = 1.
%   In finite depth: rho = g = d = 1.
%
% EXAMPLE 1. To compute a wave of steepness kH2=0.3 in infinite depth:
% [zs,ws,PP]=SSGW(inf,0.3);
%
% EXAMPLE 2. To compute a cnoidal wave with height-over-depth=0.5 and 
% length-over-depth=100:
% Hd=0.5; Ld=100; kd=2*pi/Ld; kH2=pi*Hd/Ld; [zs,ws,PP]=SSGW(kd,kH2);
% 
% EXAMPLE 3. For steep and long waves, the default number of Fourier modes
% must be increased. For instance, in order to compute a cnoidal wave with 
% height-over-depth=0.7 and length-over-depth=10000:
% Hd=0.7; Ld=10000; kd=2*pi/Ld; kH2=pi*Hd/Ld; [zs,ws,PP]=SSGW(kd,kH2,2^19);
%
% Edit the m-file for more details.

% For details of the algorithm and the notations, read:
% Clamond, D. & Dutykh, D. 2018. Accurate fast computation of steady 
% two-dimensional surface gravity waves in arbitrary depth.  
% J. Fluid Mech. 844, pp. 491-518.
%
% This m-file was written with the purpose of clarity. The notations closely 
% match those of the paper above.

% Authors: D. Clamond & D. Dutykh.
% Version: 2017-02-08.
% Revision: 2018-05-31.

%--------------------------------------------------------------------------

% Check input parameters.
if nargin<2                                             
   error('Two dimensionless scalar parameters must be provided.');
end
if kd<0 || imag(kd)~=0 || kH2<0 || imag(kH2)~=0
   error('Input scalar parameters kd and kH2 must be real and positive.');
end
if nargin<3
   N = 2048;
end
if nargin<4
   tol =1e-14;
end
if nargin<5
   Z =[];
end

% Determine depth and choose parameters.
if 1-tanh(kd) < tol                                                        % Deep water case.
   d   = inf;                                                                % Depth.
   k   = 1;                                                                  % Wavenumber.
   g   = 1;                                                                  % Acceleration due to gravity.
   lam = 1/k;                                                                % Characteristic wavelength lambda.
else                                                                       % Finite depth case.
   d   = 1;                                                                  % Depth.
   k   = kd/d;                                                               % Wavenumber.
   g   = 1;                                                                  % Acceleration due to gravity.
   lam = tanh(kd)/k;                                                         % Characteristic wavelength lambda.   
end
c02 = g*lam;                                                               % Linear phase velocity squared.
H   = 2*kH2/k;                                                             % Total wave height.
L   = pi/k;                                                                % Half-length of the computational domain (with c_r=c_e).
dal = L/N;                                                                 % Delta alpha.
dk  = pi/L;                                                                % Delta k.

% Vectors.
va  = (0:2*N-1)'*dal;                                                      % Vector of abscissas in the conformal space.
vk  = [ 0:N-1 -N:-1 ]'*dk;                                                 % Vector of wavenumbers.

% Initial guess for the solution:
Ups = (H/2)*(1+cos(k*va));                                                 % Airy solution for Upsilon.
sig = 1;                                                                   % Parameter sigma.

% Commence Petviashvili's iterations.
err  = inf;                                                                % Enforce loop entry.
iter = 0;                                                                  % Iterations counter.
tic;                                                                       % Start clocking.
while (err > tol)

  % Compute sigma and delta.      
  mUps = mean(Ups);                                                        % << Upsilon >>.
  Ys   = Ups - mUps;                                                       % Y_s.
  if d==inf                                                                % Deep water.
     sig = 1;                                                                % sigma.
     CYs = real(ifft(abs(vk).*fft(Ys)));                                     % C{ Y_s }.
     mys = -Ys'*CYs/N/2;                                                     % << y_s >>.
  else                                                                     % Finite depth.
     C_hat  =  vk.*coth((sig*d)*vk);      C_hat(1)  = 1/(sig*d);             % Operator C in Fourier space.
     S2_hat = (vk.*csch((sig*d)*vk)).^2;  S2_hat(1) = 1/(sig*d)^2;           % Operator S^2 in Fourier space.
     Ys_hat  = fft(Ys);
     E  = mean(Ys.*real(ifft(C_hat.*Ys_hat))) + (sig-1)*d;                   % Equation for sigma.
     dE = d - d*mean(Ys.*real(ifft(S2_hat.*Ys_hat)));                        % Its derivative.
     sig = sig - E/dE;                                                       % Newton new sigma.
     mys = (sig-1)*d;                                                        % << y_s >>.
  end
  del = mys - mUps;                                                        % Parameter delta.
  C_hat  =  vk.*coth((sig*d)*vk);  C_hat(1) = 1/(sig*d);                   % Updated operator C in Fourier space.
 
  % Compute Bernoulli constant B.
  Ups2  = Ups.*Ups;                                                        % Upsilon^2.
  mUps2 = mean(Ups2);                                                      % << Upsilon^2 >>.
  CUps  = real(ifft(C_hat.*fft(Ups)));                                     % C{ Upsilon }.
  CUps2 = real(ifft(C_hat.*fft(Ups2)));                                    % C{ Upsilon^2 }.
  DCU   = CUps(N+1) -  CUps(1);                                            % C{ Upsilon }_trough - C{ Upsilon }_crest.
  DCU2  = CUps2(N+1) - CUps2(1);                                           % C{ Upsilon^2 }_trough - C{ Upsilon^2 }_crest.
  Bg    = 2*del - H/sig*(1+del/d+sig*CUps(1))/DCU + DCU2/DCU/2;            % B/g.
  
  % Define linear operators in Fourier space.
  Cinf_hat = abs(vk);  Cinf_hat(1) = 0;                                    % Operator C_inf.  
  CIC_hat  = tanh((sig*d)*abs(vk));                                        % Operator C_inf o C^{-1}.      
  if d==inf, CIC_hat(1) = 1; end                                           % Regularisation.
  L_hat    = (Bg-2*del)*Cinf_hat - ((1+del/d)/sig)*CIC_hat;                % Operator L.
  IL_hat   = 1./L_hat;  IL_hat(1) = 1;                                     % Operator L^-1.
 
  % Petviashvili's iteration.
  Ups_hat = fft(Ups);                                                      % Fourier transform of Upsilon.
  CUps_hat = C_hat.*Ups_hat;
  LUps = real(ifft(L_hat.*Ups_hat));                                       % L{Upsilon}.
  Ups2_hat = fft(Ups.*Ups);                                                % Fourier transform of Upsilon^2.
  NUps_hat = CIC_hat.*fft(Ups.*real(ifft(CUps_hat)));
  NUps_hat = NUps_hat + Cinf_hat.*Ups2_hat/2;                              % Nonlinear term in Fourier space.
  NUps = real(ifft(NUps_hat));                                             % N{ Upsilon }.
  S = (Ups'*LUps)/(Ups'*NUps);                                             % Weight.
  U = S*S*real(ifft(NUps_hat.*IL_hat));                                    % New Upsilon.
  U = H * ( U - U(N+1) ) / ( U(1) - U(N+1) );                              % Enforce mean value.
  
  % Update values.
  err = norm(U-Ups,inf);                                                   % Error measured by the L_inf norm.
  Ups = U;                                                                 % New Upsilon.
  iter = iter+1;
  
end
time = toc;


% Post processing.
IH_hat = -1i*coth(sig*d*vk);  IH_hat(1) = 0;                               % Inverse Hilbert transform.
Ys  = Ups - mean(Ups);
Ys_hat  = fft(Ys);
CYs = real(ifft(C_hat.*Ys_hat));
Xs  = real(ifft(IH_hat.*Ys_hat));
mys = -Ys'*CYs/N/2;
Zs  = Xs + 1i*Ys;
dZs = ifft(1i*vk.*fft(Zs));
zs  = va + 1i*mys + Zs;
dzs = 1 + dZs;
B   = g*Bg;
ce  = sum( (1+CYs)./abs(dzs).^2 )/2/N;
ce  = sqrt(B/ce);
cs  = sig*ce;
ws  = -ce./dzs;
a   = max(imag(zs));
b   = -min(imag(zs));

xs = [ real(zs(N+1:end))-2*pi/k ; real(zs(1:N)) ];
ys = [ imag(zs(N+1:end)) ; imag(zs(1:N)) ];

if d==inf 
    Bce2d = 0;
    IC = 1./abs(vk); IC(1) = 0;
else
    Bce2d = (B-ce^2)*d;
    IC = tanh(vk*sig*d)./vk; IC(1) = sig*d;                                % Inverse C-operator.
end
ydx  = real(dzs).*imag(zs);
intI = -ce*mys;                                                            % Impulse.
intV = mean(ydx.*imag(zs))*g/2;                                            % Potential energy.
intK = intI*ce/2;                                                          % Kinetic energy.
intSxx = 2*ce*intI - 2*intV + Bce2d;                                       % Radiation stress.
intS = intSxx - intV + g*d^2/2;                                            % Momentum flux.
intF = Bce2d*ce/2 + (B+ce^2)*intI/2 + (intK-2*intV)*ce;                    % Energy flux.
cg   = intF/(intK+intV);                                                   % Group velocity.
K1   = imag(zs)'*(imag(zs)/2-Bg)/2/N;
ICydx  = real(ifft(IC.*fft(ydx))); 
errfun = norm( imag(zs) - Bg + sqrt(Bg^2+2*K1-2*ICydx) ,inf);              % Residual.

% Output data.
PP(1)=d; PP(2)=k; PP(3)=H; PP(4)=ce; PP(5)=cs; PP(6)=B; PP(7)=a; PP(8)=b; 
PP(9)=intI; PP(10)=intV; pp(11)=intK; PP(12)=intSxx; PP(13)=intS;
PP(14)=intF; PP(15)=cg;

%% added by Ling Zhu, July 12, 2023
% Output at desired locations in the bulk.
if length(Z)>1
    Ceta  = real(ifft(C_hat.*fft(imag(zs))));               % C(eta)
    dexi   = (ifft(1i*vk.*fft(imag(zs))));              % d eta / d xi
    
    if isempty(Z)==0
       Z = Z(:);
       LZ = length(Z);
       W = zeros(LZ,1); 
       dzsm1 = ce/cs - 1 + Ceta + 1i*dexi ;
       factor = 0.25*1i * k/pi; 
       
       for n=1:LZ
            zzz = dzsm1./tan(0.5*k*zs-0.5*k*Z(n)) - conj(dzsm1)./tan(0.5*conj(k*zs)-1i*kd-0.5*k*Z(n))  ;
            W(n)  = factor*dal*ce*sum(zzz) ;
        end
    end
else
    W = [];
end

