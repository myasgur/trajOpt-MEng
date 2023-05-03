function Ydot = cr3bp_EOM_time(t,Y)

global mu T_max c

%%%%%%%%%%%%%%%%%% Unpack State %%%%%%%%%%%%%%%%%%%%%%%%

r = Y(1:3); v = Y(4:6); m = Y(7); 
lr = Y(8:10); lv = Y(11:13); lm = Y(14);

x = r(1); y = r(2); z = r(3);
vx = v(1); vy = v(2); vz = v(3);
lvx = lv(1); lvy = lv(2); lvz = lv(3); lv_mag = norm(lv); 

%%%%%%% Switching Function ( Throttle Control ) %%%%%%%%
St = -c.*lv_mag./m - lm;
% u = 0.5*(1-tanh(S./rho)); %hyperbolic tangent smoothing
    
if St > 0
    u = 0;
elseif St < 0
    u = 1;
else
    u = 0.5; 
end

%%%%%%%%%%%%%%%% Common Distances %%%%%%%%%%%%%%%%%%%
r1 = sqrt(( mu+x  ).^2+y.^2+z.^2);
r2 = sqrt(( mu+x-1).^2+y.^2+z.^2);
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
G = zeros(3); %dGlvdr = zeros(3);

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

%%%%%%%%%%%%%%%% Equations of Motion %%%%%%%%%%%%%%%%%%%
Ydot = [v;
        g+h-(lv./lv_mag).*u.*T_max./m;
        -u.*T_max./c;
        -G'*lv;
        -lr-H'*lv;
        -lv_mag.*u.*T_max./m.^2];
end