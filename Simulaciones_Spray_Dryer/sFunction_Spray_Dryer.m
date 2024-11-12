%% Creado por Colomba Dazzarola el 4/4/2024

function [sys, x0,str,ts] = sFunction_Spray_Dryer(t,x,u,flag)
%EL código siguiente no se modifica
% global qt qc; 

switch flag
   case 0
   [sys,x0,str,ts] = mdlInitializeSizes;
   case 1 
   sys = mdlDerivatives(t,x,u);
   case 3
   sys = mdlOutputs(t,x,u);
   case {2,4,9}
   sys=[];
   otherwise
   error(['Unhandled flag = ',num2str(flag)]);
end
end

function [sys,x0,str,ts] = mdlInitializeSizes
% Declaration of number of differential equations to integrate, of varying
% input variables that enter the macro and outputs that it will have.

sizes=simsizes;
sizes.NumContStates  = 3;        % Number of differential equations to integrate
sizes.NumDiscStates  = 0;  
sizes.NumOutputs     = 4;        % Number of output variables that the code have
sizes.NumInputs      = 7;        % Number of input variables to the code
sizes.DirFeedthrough = 1;        % If input information is required in the calculation of outputs
sizes.NumSampleTimes = 1; 
sys = simsizes(sizes);
str = [];
ts = [0 0];
x0 = [1.52e-2   3   100+273.15]; % Initial conditions

end 

function sys = mdlDerivatives(~,x,u)

% State variables 
Ho      =   x(1);       % Air absolute humidity outlet air (kg vapor/kg dry air)
Xo      =   x(2);       % Water content outlet solid (kg water/kg solids)
T4      =   x(3);       % Temperature air outlet (K)

% Input variables
G       =   u(1);       % Dry air flow-rate (kg/s)
T_i     =   u(2);       % Temperature inlet (K)
Hi      =   u(3);       % Air absolute humidity air inlet (kg vapor/kg dry air)

F       =   u(4);       % Solids flow-rate (kg/s)
T_F     =   u(5);       % Temperature liquid feed (K)
Xi      =   u(6);       % Water content inlet (kg water/kg solids)
ri      =   u(7);

% Model constants: values used in Farias, 1995
a       =   1000;       % Empiral constant for semi-skimmed milk (kg/m3)
b       =   290;        % Empiral constant for semi-skimmed milk (kg/m3)

% Model constants: values used in Zaror and Perez-Correa, 1991
Nu      =   2;            % Nusselt number
m       =   1.5;          % Solid moisture equilibrium constant
M_A     =   7;            % Air hold-up in the drier chamber (kg)
T_sat   =   43 + 273.15;  % Temperature at saturation (K)    
Xoc     =   0.25;         % Water content critical (kg water/kg solids)

% Algebraic equation: mass of solids in droplet
ms      =   ri^3 * (4*pi * (a*(1+Xi) + b))/(3 * (1 + Xi)^2); % Mass of solids in a droplet (kg solids)

% Algebraic Equations: relations used in Zaror and Perez-Correa, 1991
De      =   (1.38*10^(-11)*T4 - 1.55*10^(-9));       % Effective diffusivity (m2/s)
kA      =   (7.4*10^(-8)*T4 + 4.19*10^(-6));         % Air termal conductivity (kW/m/K)
lambda_sat = 3180.14 - 2.508 * T_sat;                % Latent heat of vaporisation (kJ/kg)

% Algebraic Equations: Mass fraction
%xa_1    =   Xo / (1 - Xo);      % Mass fraction of water in current 1 (kg water/kg droplet)
%xs_1    =   1 - xa_1;           % Mass fraction of solids in current 1 (kg solids/kg droplet)
yv_4    =   Ho / (1 - Ho);      % Mass fraction of water vapour in current 4 (kg vapour/kg humid air)
%yA_4    =   1 - yv_4;           % Mass fraction of dry air in current 4 (kg dry air/kg humid air)
yv_3    =   Hi / (1 - Hi);      % Mass fraction of water vapour in current 3 (kg vapour/kg humid air)

% Algebraic Equations: Mass flowrates 
F_1     =   F*(1 + Xi);         % Mass flowrate of current 1 (kg/s)
%F_2     =   F*(1 + Xo);         % Mass flowrate of current 2 (kg/s)
F_3     =   G*(1 + Hi);         % Mass flowrate of current 3 (kg/s)
F_4     =   G*(1 + Ho);         % Mass flowrate of current 4 (kg/s)

% Algebraic Equations: Specific heat capacity
% Regression presented in Munir (2016). Engineering toolbox aprox.
% 3.98-4.02
Cp_1    =   (3774.48 + 1.15 * (T_i - 273.15) + 3.93*10^(-3) * (T_i - 273.15)^2)/1000; 
% Regression presented in Geankoplis (1993)
Cp_3    =   1.005 + 1.88 * Hi; % Calor humedo o capacidad calorifica del aire humedo (kJ/kg aire seco/K)
Cp_4    =   1.005 + 1.88 * Ho; % Calor humedo o capacidad calorifica del aire humedo (kJ/kg aire seco/K)
% Assume regression
Cv_4    =   0.718 + 1.4108 * Ho;

% Algebraic Equations: Mass hold-up gas (dry air + water vapour) (kg gas)
Mg      =   M_A * (1 + Ho);

% Algebraic Equations: used in Farias 1995
if Xo >= Xoc
    rd = ((3 * (1 + Xo)^2 * ms)/(4 * pi * (a*(1 + Xo) + b)))^(1/3);     % droplet radio (m)
    adrop = 4 * pi * rd^2;                                              % droplet superficial area (m^2)
    h = Nu * kA / (2*rd);                                               % heat transfer coefficent (kW/m2/K)
else 
    rd = ((3 * (1 + Xoc)^2 * ms)/(4 * pi * (a*(1 + Xoc) + b)))^(1/3);   % droplet radio (m)
    adrop = 4 * pi * rd^2;                                              % droplet superficial area (m^2)
    h = Nu * kA / (2*rd);                                               % heat transfer coefficent (kW/m2/K)
end

% Differential equations: mass balances Ho, Xo, T_4
sys(1)  = G/M_A * (Hi - Ho) + F/M_A * (Xi - Xo);

if Xo >= Xoc
    sys(2)  = (h / lambda_sat * adrop / ms * (T_sat - T4)); 
else
    sys(2)  = (-4*pi^2/(2*rd)^2) * De * (Xo - m*Ho);
end

sys(3) = (F_1 * Cp_1 * (T_F - T_sat) + F_3 * Cp_3 * (T_i - T_sat) - F_4 * Cp_4 * (T4 - T_sat) + lambda_sat * (F_3 * yv_3 - F_4 * yv_4))/(Mg * Cv_4);

end

function sys = mdlOutputs(~,x,u)

% State variables 
Ho      =   x(1);       % Air absolute humidity outlet air (kg vapor/kg dry air)
Xo      =   x(2);       % Water content outlet solid (kg water/kg solids)
T4      =   x(3);       % Temperature air outlet (K)

% Input variables required
Xi      =   u(6);       % Water content inlet (kg water/kg solids)
ri      =   u(7);

% Algebraic Equations
a       =   1000;       % Empiral constant for semi-skimmed milk (kg/m3)
b       =   290;        % Empiral constant for semi-skimmed milk (kg/m3)
Xoc     =   0.25;         % Water content critical (kg water/kg solids)
ms      =   ri^3 * (4*pi * (a*(1+Xi) + b))/(3 * (1 + Xi)^2);    % Mass of solids in a droplet (kg solids)

if Xo >= Xoc
    rd = ((3 * (1 + Xo)^2 * ms)/(4 * pi * (a*(1 + Xo) + b)))^(1/3);     % droplet radio (m)
else 
    rd = ((3 * (1 + Xoc)^2 * ms)/(4 * pi * (a*(1 + Xoc) + b)))^(1/3);   % droplet radio (m)
end

% Outputs:
sys = [T4 Ho Xo rd];

end
