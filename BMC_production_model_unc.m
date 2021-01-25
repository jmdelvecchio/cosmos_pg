function [P_total_10Be, P_total_26Al] = BMC_production_model(MassDepth)

% muon production for bear meadows core
Latitude = 40.7481;     % bear meadows core latitude
Longitude = -77.7457;   % bear meadows core longitude
Elevation = 562;        % bear meadows core elevation from lidar

Pressure = ERA40atm(Latitude,Longitude,Elevation);
RigidityCutoff = lat_to_Rc(Latitude);
SpallationAttenuationLength = interpLsp(Pressure,RigidityCutoff);
SpallationScalingFactor = stone2000(Latitude,Pressure,1);

Pslhl_10Be = 4.01;  % from Borchers 2016 for Stone scaling scheme
Pslhl_26Al = 27.93; % from Borchers 2016 for Stone scaling scheme

P_spallation_10Be = SpallationScalingFactor.*Pslhl_10Be.*exp(-MassDepth./SpallationAttenuationLength);
P_spallation_26Al = SpallationScalingFactor.*Pslhl_26Al.*exp(-MassDepth./SpallationAttenuationLength);

consts10Be.Natoms = 2.006e22;               % atoms O per gram of quartz
consts10Be.k_neg = 0.00191.*0.704.*0.1828;  % from Balco 2017 Beacon Heights calibration, incorporates Fc and Fd from Heisinger 2002
consts10Be.sigma0 = 0.280e-30;              % cm2, from Balco 2017 Beacon Heights calibration

consts26Al.Natoms = 1.003e22;               % atoms Si per gram of quartz
consts26Al.k_neg = 0.0133.*0.296.*0.6559;   % from Balco 2017 Beacon Heights calibration, incorporates Fc and Fd from Heisinger 2002
consts26Al.sigma0 = 3.89e-30;               % cm2, from Balco 2017 Beacon Heights calibration


P_muon_10Be0 = P_mu_total_alpha1(MassDepth, Pressure, consts10Be,'no')';
P_muon_26Al0 = P_mu_total_alpha1(MassDepth, Pressure, consts26Al,'no')';

P_muon_10Be = P_muon_10Be0 + (P_muon_10Be0 * 0.1 * randn(1,1)); %10% uncertainty
P_muon_26Al = P_muon_26Al0 + (P_muon_26Al0 * 0.1 * randn(1,1)); %10% uncertainty


P_total_10Be = P_spallation_10Be + P_muon_10Be;
P_total_26Al = P_spallation_26Al + P_muon_26Al;