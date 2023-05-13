function [g,h]=ghFunctions(r,v)

global mu

x = r(1); y = r(2); z = r(3);
vx = v(1); vy = v(2); vz = v(3);

%%%%%%%%%%%%%%%% Common Distances %%%%%%%%%%%%%%%%%%%
r1 = sqrt(( mu+x  ).^2+y.^2+z.^2);
r2 = sqrt(( mu+x-1).^2+y.^2+z.^2);
r13 = r1.^3;
r23 = r2.^3;

%%%%%%%%%%%%%% g(r) and h(v) functions %%%%%%%%%%%%%%%%
g = [(x - (1-mu).*(x+mu)./r13 - mu.*(x+mu-1)./r23);
     (y - (1-mu).*y./r13 - mu.*y./r23);
     (-(1-mu).*z./r13 - mu.*z./r23)];
h = [2*vy; -2*vx; 0];

end