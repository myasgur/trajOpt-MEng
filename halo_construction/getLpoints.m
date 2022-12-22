function Ls = getLpoints(mu)
% INPUT: 
%           mu = mass parameter of 3-body system
% OUTPUT: 
%           Ls = x-positions of collinear Lagrange points, L1-3 [3x1]

syms gam
f = @(x) x - (1 -mu)*(x+mu)./abs(x+mu).^3 - mu*(x - 1+mu)./abs(x - 1 + mu).^3;
g1 = gam^5 - (3-mu)*gam^4 + (3-2*mu)*gam^3 - mu*gam^2 + 2*mu*gam - mu==0;
g2 = gam^5 + (3-mu)*gam^4 + (3-2*mu)*gam^3 - mu*gam^2 - 2*mu*gam - mu==0;

gamma1 = vpasolve(g1,0.1-mu);
gamma2 = vpasolve(g2,0.1-mu);

L1 = (1-mu)-double(gamma1(1));
L2 = (1-mu)+double(gamma2(1));
L3 = fzero(f,[-1-mu-0.02,-1-mu+0.02]);
Ls = [L1,L2,L3];
end