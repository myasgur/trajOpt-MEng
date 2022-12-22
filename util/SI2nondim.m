function [Tmax_nondim,c_nondim]=SI2nondim(units,Tmax_SI,c_SI)
% Nondimensional units parameters for Earth-Sun-Spacecraft 3-body system
MUnit = units(1); % nondim mass unit (-->kg)
MscUnit = units(2); % nondim S/C mass unit
TUnit = units(3); % nondim time unit (-->seconds)
DUnit = units(4); % nondim distance unit (-->km)
VUnit = units(5); % nondim speed unit (-->km/s)
FUnit = 1000*MscUnit*DUnit/TUnit^2; %nondim force (-->N)

% m_nondim = m_SI./MscUnit;
Tmax_nondim=Tmax_SI./FUnit;
c_nondim=c_SI./VUnit;
    
end