function [ Biorefinery_array] = makeBiorefineryArray(  )
%MAKEBIOREFINERYARRAY Summary of this function goes here
%   Detailed explanation goes here

%% Set n_refineries here, and specifics for each refinery below.
n_refineries = 3;

%% Preallocation
Biorefinery_array(1:n_refineries) = Biorefinery;

%% Set biorefinery specifications here
Biorefinery_array(1) = Biorefinery('Present day biorefinery. Commercial-scale cellulosic biorefineries', 2020, 0.433, 0.4035, 0, 0.0296, 0.263, 0.737, 0, 0.00025, 0.614, 0.386, 'Laser et al. 2009a, Field et al. 2020');
Biorefinery_array(2) = Biorefinery('Future biorefinery', 2020, 0.712, 0.541, 0.158, 0.013, 0.461, 0.522, 0.017, -999, 0.614, 0.386, 'Laser et al. 2009b, Field et al. (2020)');
Biorefinery_array(3) = Biorefinery('Future BECCS', 2020, 0.676, 0.541, 0.158, -0.023, 0.461, 0.041, 0.498, -999, 0.614, 0.386, 'Field et al. 2020');




end

