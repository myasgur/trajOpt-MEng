function Zdot = varEqsTrajOpt(~,Z,params)

%%%%%%%%%%%%%%% Problem Setup %%%%%%%%%%%%%%%%%%%%%%%%%
N = 14; % number of states + costates
Y = Z(1:N); % state vector
PHI = Z(N+1:(N^2+N)); 
phi = reshape(PHI,N,N); % state transition matrix

mu = params.mu; % sun-earth mass parameter
rho = params.rho; % continuation parameter
epsilon = params.epsilon; % homotopic performance index
Tmax_SI = params.Tmax; % maximum thrust (N)
c_SI = params.c; % exhaust velocity (km/s)

MUnit = params.MUnit; % nondim mass unit (-->kg)
MscUnit = params.MscUnit; % nondim S/C mass unit
TUnit = params.TUnit; % nondim time unit (-->seconds)
DUnit = params.DUnit; % nondim distance unit (-->km)
VUnit = params.VUnit; % nondim speed unit (-->km/s)
units = [MUnit,MscUnit,TUnit,DUnit,VUnit];

%%%%%%%%%%%%%%% Nondimensionalize %%%%%%%%%%%%%%%%%%%%%
[Tmax,c]=SI2nondim(units,Tmax_SI,c_SI);

%%%%%%%%%%%%%%%%%% Unpack State %%%%%%%%%%%%%%%%%%%%%%%%

r = Y(1:3); v = Y(4:6); m = Y(7); 
lr = Y(8:10); lv = Y(11:13); lm = Y(14);

x = r(1); y = r(2); z = r(3);
vx = v(1); vy = v(2); vz = v(3);
lvx = lv(1); lvy = lv(2); lvz = lv(3); lv_mag = norm(lv); 

%%%%%%% Switching Function ( Throttle Control ) %%%%%%%%
S = 1-c.*lv_mag./m - lm;
% u = 0.5*(1-tanh(S./rho)); %hyperbolic tangent smoothing
    
if S > epsilon
    u = 0;
elseif S < -epsilon
    u = 1;
else
    u = 0.5*(1-tanh(S./rho)); %hyperbolic tangent smoothing
end

%%%%%%%%%%%%%%%% Common Distances %%%%%%%%%%%%%%%%%%%
r1 = sqrt((( mu+x  ).^2+y.^2+z.^2));
r2 = sqrt((  mu+x-1).^2+y.^2+z.^2);
r13 = r1.^3;
r23 = r2.^3;
r15 = r1.^5;
r25 = r2.^5;

%%%%%%%%%%%%%% g(r) and h(v) functions %%%%%%%%%%%%%%%%
g = [(x - (1-mu).*(x+mu)./r13 - mu.*(x+mu-1)./r23);
     (y - (1-mu).*y./r13 - mu.*y./r23);
     (-(1-mu).*z./r13 - mu.*z./r23)];
h = [2*vy; -2*vx; 0];

%%%%%%%%%%%%%%%% G and H Matrices %%%%%%%%%%%%%%%%%%%%%%
G = zeros(3); dGlvdr = zeros(3);

G(1,1) = 1-(1-mu)./r13+3.*(1-mu).*(x+mu).^2./r15-mu./r23+3.*mu*(x+mu-1).^2./r25;
G(2,2) = 1-(1-mu)./r13+3.*(1-mu).*y.^2./r15-mu./r23+3.*mu*y.^2./r25;
G(3,3) =  -(1-mu)./r13+3.*(1-mu).*z.^2./r15-mu./r23+3.*mu*z.^2./r25;
G(1,2) = 3.*(1-mu).*(x+mu).*y./r15+3.*mu*(x+mu-1).*y./r25;
G(2,1) = G(1,2);
G(1,3) = 3.*(1-mu).*(x+mu).*z./r15+3.*mu*(x+mu-1).*z./r25;
G(3,1) = G(1,3);
G(2,3) = 3.*(1-mu).*y.*z./r15+3.*mu.*y.*z./r25;
G(3,2) = G(2,3);

H = [0 2 0;
    -2 0 0;
     0 0 0];

Gxx = (3.*(mu+x).*(mu-1)./r15 - 3.*mu.*(mu+x-1)./r25);
Gyy = (3.*y.*(mu-1)./r15 - 3.*mu.*y./r25);
Gzz = (3.*z.*(mu-1)./r15 - 3.*mu.*z./r25);

dGlvdr(1,:) = [2.*lvx.*Gxx + lvy.*Gyy + lvz.*Gzz;
                          lvy.*Gxx;
                          lvz.*Gxx];

dGlvdr(2,:) = [           lvx.*Gyy;
               lvx.*Gxx - 2.*lvy*Gyy - lvz.*Gzz;
                          lvz.*Gyy];

dGlvdr(3,:) = [           lvx.*Gzz;
                          lvy.*Gzz;
               lvx.*Gxx - lvy*Gyy - 2.*lvz*Gzz];

%%%%%%%%%%%%%%%% Equations of Motion %%%%%%%%%%%%%%%%%%%
Ydot = [v;
        g+h-(lv./lv_mag).*u.*Tmax./m;
        -u.*Tmax./c;
        -G.'*lv;
        -lr-H.'*lv;
        -lv_mag.*u.*Tmax./m.^2];

%%%%%%%%%%%%%%%%%% Jacobian Matrix %%%%%%%%%%%%%%%%%%%%
Df = zeros(N);
Df(1:3,4:6) = eye(3);
Df(4:6,1:3) = G;
Df(4:6,4:6) = H;
Df(4:6,7) = (lv./lv_mag).*u.*Tmax./m.^2;
Df(4:6,11:13) = -(u.*Tmax./m).*((eye(3)./lv_mag) - (lv*lv.'/lv_mag^3));
Df(8:10,1:3) = -dGlvdr;
Df(8:10,11:13) = -G.';
Df(11:13,8:10) = -eye(3);
Df(11:13,11:13) = -H.';
Df(14,7) = 2.*lv_mag.*u.*Tmax./m.^3;
Df(14,11:13) = -(lv.'./lv_mag).*u.*Tmax./m.^2;

%%%%%%%%%%%%%%% Variational Equations %%%%%%%%%%%%%%%%%%

phidot = Df*phi;

%%%%%%%%%%%%%%%%  Reshape and Collect %%%%%%%%%%%%%%%%%%

Zdot = [Ydot(:);phidot(:)];



end