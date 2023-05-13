function [value,isterminal,dir] = cr3bp_event(~,z)
    global DU mu

    r1sq = sqrt((mu+z(1)).^2+z(2).^2+z(3).^2);
    r2sq = sqrt((mu+z(1)-1).^2+z(2).^2+z(3).^2);
    value = [r1sq - 6378.137/(DU);
             r2sq - 1738.1/(DU)]; 
    isterminal = [1;1]; % terminate
    dir = [0;0]; % from above or below
end