function lams = psoSearchVars(x)
global l0
% convert search variable to angles (Eq (32))
beta = zeros(size(x)); %intialize
beta(1:3) = 0.5*pi*x(1:3);
beta(4:5) = pi*(x(4:5)-0.5);
beta(6:7) = 2*pi*x(6:7);

% convert angles to costates (Eq (33))
l0 = sin(beta(1));

lr0 = cos(beta(1))*cos(beta(2))*cos(beta(3))*...
    [cos(beta(4))*cos(beta(6));
     cos(beta(4))*sin(beta(6)); 
     sin(beta(4))];

lv0 = cos(beta(1))*cos(beta(2))*sin(beta(3))*...
    [cos(beta(5))*cos(beta(7));
     cos(beta(5))*sin(beta(7)); 
     sin(beta(5))];

lm0 = cos(beta(1))*sin(beta(2));

lams = [lr0;lv0;lm0];

end