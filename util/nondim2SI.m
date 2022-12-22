function [Tmax_SI,c_SI]=nondim2SI(units,Tmax_nondim,c_nondim)
% Nondimensional units parameters for Earth-Sun-Spacecraft 3-body system
MUnit = units(1); % nondim mass unit (-->kg)
MscUnit = units(2); % nondim S/C mass unit
TUnit = units(3); % nondim time unit (-->seconds)
DUnit = units(4); % nondim distance unit (-->km)
VUnit = units(5); % nondim speed unit (-->km/s)
FUnit = 1000*MscUnit*DUnit/TUnit^2; %nondim force (-->N)

% m_SI = m_nondim.*MscUnit;
Tmax_SI=Tmax_nondim.*FUnit;
c_SI=c_nondim.*VUnit;
    
end