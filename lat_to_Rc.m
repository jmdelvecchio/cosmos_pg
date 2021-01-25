function out = lat_to_Rc(lat)

% This applies Lifton polynomial to convert geographic lat to Rc. 

latr = lat.*pi./180; % convert to radians

% Define polynomial in cos(lat) from LSD2014, Equation 2
pRc = [-448.004 1189.18 -1152.15 522.061 -103.241 6.89901 0];
% Apply polynomial formula
out = polyval(pRc,cos(latr));