function plot_potentials_as_function_of_emission_factor(FEM, Nor_mask)

binary_pixels_to_consider = FEM.pe > 0;

x_vec_emission_factor = [-400:5:400];

prefix_output = 'High_res_';
output_folder = 'Output/';

prefix_output = [output_folder prefix_output];

time = FEM.time_ac;
idx_cohort = zeros(1,length(time));

time_cohort = [2010 2000 1990];

for i = 1:length(time)
    for j = 1:length(time_cohort)
        if time(i) >= time_cohort(j)
            idx_cohort(i) = j;
            break
        end
    end

end


fe_elec_potential = zeros(length(time_cohort),length(x_vec_emission_factor));
fe_ft_potential = zeros(length(time_cohort),length(x_vec_emission_factor));
fe_elec_ccs_potential = zeros(length(time_cohort),length(x_vec_emission_factor));
fe_ft_ccs_potential = zeros(length(time_cohort),length(x_vec_emission_factor));

fe_elec_potential_Nor = zeros(length(time_cohort),length(x_vec_emission_factor));
fe_ft_potential_Nor = zeros(length(time_cohort),length(x_vec_emission_factor));
fe_elec_ccs_potential_Nor = zeros(length(time_cohort),length(x_vec_emission_factor));
fe_ft_ccs_potential_Nor = zeros(length(time_cohort),length(x_vec_emission_factor));



for i = 1:length(x_vec_emission_factor)
    fprintf([num2str(x_vec_emission_factor(i)) ', '])
    for t = 1:length(FEM.time_ac)

        fe_elec_this  = FEM.fe_elec_after_aban_year(:,:,t);
        fe_elec_ccs_this = FEM.fe_elec_ccs_after_aban_year(:,:,t);
        fe_ft_this = FEM.fe_ft_after_aban_year(:,:,t);
        fe_ft_ccs_this = FEM.fe_ft_ccs_after_aban_year(:,:,t);

        binary_elec_this = FEM.emission_factor_electricity(:,:,t) <= x_vec_emission_factor(i);
        binary_elec_ccs_this = FEM.emission_factor_electricity_CCS(:,:,t) <= x_vec_emission_factor(i);
        binary_ft_this = FEM.emission_factor_FT(:,:,t) <= x_vec_emission_factor(i);
        binary_ft_ccs_this = FEM.emission_factor_FT_CCS(:,:,t) <= x_vec_emission_factor(i);
        
        fe_elec_potential(idx_cohort(t),i) = fe_elec_potential(idx_cohort(t),i) + sum(sum(fe_elec_this(binary_elec_this)));
        fe_elec_ccs_potential(idx_cohort(t),i) = fe_elec_ccs_potential(idx_cohort(t),i) + sum(sum(fe_elec_ccs_this(binary_elec_ccs_this)));

        fe_ft_potential(idx_cohort(t),i) = fe_ft_potential(idx_cohort(t),i) + sum(sum(fe_ft_this(binary_ft_this)));
        fe_ft_ccs_potential(idx_cohort(t),i) = fe_ft_ccs_potential(idx_cohort(t),i) + sum(sum(fe_ft_ccs_this(binary_ft_ccs_this)));
        
        %Nor
        fe_elec_potential_Nor(idx_cohort(t),i) = fe_elec_potential_Nor(idx_cohort(t),i) + sum(sum(fe_elec_this(binary_elec_this & Nor_mask)));
        fe_elec_ccs_potential_Nor(idx_cohort(t),i) = fe_elec_ccs_potential_Nor(idx_cohort(t),i) + sum(sum(fe_elec_ccs_this(binary_elec_ccs_this & Nor_mask)));

        fe_ft_potential_Nor(idx_cohort(t),i) = fe_ft_potential_Nor(idx_cohort(t),i) + sum(sum(fe_ft_this(binary_ft_this & Nor_mask)));
        fe_ft_ccs_potential_Nor(idx_cohort(t),i) = fe_ft_ccs_potential_Nor(idx_cohort(t),i) + sum(sum(fe_ft_ccs_this(binary_ft_ccs_this & Nor_mask)));


    end

end

%% PLOT

transparancy = 0.5;
line_width = 1;

% Emission factors, alternative techs
emf_solar_wind_elec = [2 16];
emf_fossil_with_ccs_elec = [44 73];
emf_natural_gas_elec = [136 146];
emf_coal_elec = [220 259];
emf_petrol_fuel = 92.4;
emf_diesel_fuel = 93.9;

% Colors
color_solar_wind = [1 1 0];
color_coal =  [0 0 0];
color_ng = [0.4940 0.1840 0.5560];
color_fossil_with_ccs = [0.8500 0.3250 0.0980];
color_diesel = [0.6350 0.0780 0.1840];

% REF Scarlat et al. (2022), Applied energy
emf_Norway_elec_prod = 7.8;
emf_EU27_elec_prod = 86.1;

legend_text={'2011 - 2020', '2001 - 2010', '1993 - 2000'};

figure 

%% Plot bioelectricity
area(x_vec_emission_factor, 10^-9*fe_elec_potential')
hold on
legend(legend_text, 'location', 'northwest')
xline(0);

%Solar/wind
x_points = [emf_solar_wind_elec(1) emf_solar_wind_elec(1) emf_solar_wind_elec(2) emf_solar_wind_elec(2)];
y_points = [0 12 12 0];

a=fill(x_points, y_points, color_solar_wind);
a.FaceAlpha = transparancy;

%fossil elec with ccs
x_points = [emf_fossil_with_ccs_elec(1) emf_fossil_with_ccs_elec(1) emf_fossil_with_ccs_elec(2) emf_fossil_with_ccs_elec(2)];

a=fill(x_points, y_points, color_fossil_with_ccs);
a.FaceAlpha = transparancy;

%natural gas
x_points = [emf_natural_gas_elec(1) emf_natural_gas_elec(1) emf_natural_gas_elec(2) emf_natural_gas_elec(2)];

a=fill(x_points, y_points, color_ng);
a.FaceAlpha = transparancy;

%coal
x_points = [emf_coal_elec(1) emf_coal_elec(1) emf_coal_elec(2) emf_coal_elec(2)];

a=fill(x_points, y_points, color_coal);
a.FaceAlpha = transparancy;

ylabel('EJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')

title('Bioelectricity')

filename = [prefix_output 'bioelectricity.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_bioelectricity.mat', 'fe_elec_potential', 'x_vec_emission_factor');

%% Plot bioelectricity with ccs

figure 
area(x_vec_emission_factor, 10^-9*fe_elec_ccs_potential')
hold on
ylabel('EJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text, 'location', 'northwest')
title('BECCS - Bioelectricity')

xline(0);

%Solar/wind
x_points = [emf_solar_wind_elec(1) emf_solar_wind_elec(1) emf_solar_wind_elec(2) emf_solar_wind_elec(2)];
y_points = [0 8 8 0];
color = [1 1 0];
a=fill(x_points, y_points, color);
a.FaceAlpha = transparancy;

%fossil elec with ccs
x_points = [emf_fossil_with_ccs_elec(1) emf_fossil_with_ccs_elec(1) emf_fossil_with_ccs_elec(2) emf_fossil_with_ccs_elec(2)];
color = [0.8500 0.3250 0.0980];
a=fill(x_points, y_points, color);
a.FaceAlpha = transparancy;

%natural gas
x_points = [emf_natural_gas_elec(1) emf_natural_gas_elec(1) emf_natural_gas_elec(2) emf_natural_gas_elec(2)];
color = [0.4940 0.1840 0.5560];
a=fill(x_points, y_points, color);
a.FaceAlpha = transparancy;

%coal
x_points = [emf_coal_elec(1) emf_coal_elec(1) emf_coal_elec(2) emf_coal_elec(2)];
color = [0 0 0];
a=fill(x_points, y_points, color);
a.FaceAlpha = transparancy;


filename = [prefix_output 'bioelectricity_ccs.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_bioelectricity_ccs.mat', 'fe_elec_ccs_potential', 'x_vec_emission_factor');

%% Plot FT diesel



figure 
area(x_vec_emission_factor, 10^-9*fe_ft_potential')
hold on
ylabel('EJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text, 'location', 'northwest')
title('Fischer-Tropsch diesel')

xline(0)

%xline(emf_petrol_fuel, 'color', color_diesel,'LineWidth', line_width)
xline(emf_diesel_fuel, 'color', color_diesel,'LineWidth', line_width)

filename = [prefix_output 'ft_diesel.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_ft_diesel.mat', 'fe_ft_potential', 'x_vec_emission_factor');

%% Plot FT diesel ccs

figure 
area(x_vec_emission_factor, 10^-9*fe_ft_ccs_potential')
hold on
ylabel('EJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text,'location', 'northwest')
title('BECCS - Fischer-Tropsch diesel')

xline(0);
%xline(emf_petrol_fuel, 'color', color_diesel,'LineWidth', line_width)
xline(emf_diesel_fuel,'color', color_diesel,'LineWidth', line_width);


filename = [prefix_output 'ft_diesel_ccs.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_ft_diesel_ccs.mat', 'fe_ft_ccs_potential', 'x_vec_emission_factor');

%% NORWAY

%% Bioelectricity NOR
figure 

area(x_vec_emission_factor, 10^-6*fe_elec_potential_Nor')
hold on
ylabel('PJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text, 'location', 'northwest')
title('Bioelectricity')

xline(0)

xline(emf_Norway_elec_prod, 'b','LineWidth', line_width)
xline(emf_EU27_elec_prod, 'C','LineWidth', line_width)
%xline(emf_natural_gas_elec, 'color','#7E2F8E')

%natural gas
x_points = [emf_natural_gas_elec(1) emf_natural_gas_elec(1) emf_natural_gas_elec(2) emf_natural_gas_elec(2)];
y_points= ([0 2.5 2.5 0]);

a=fill(x_points, y_points, color_ng);
a.FaceAlpha = transparancy;

filename = [prefix_output 'Norway_bioelectricity.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_bioelectricity_Norway.mat', 'fe_elec_potential_Nor', 'x_vec_emission_factor');

%% Bioelectricity with CCS NOR

figure 
area(x_vec_emission_factor, 10^-6*fe_elec_ccs_potential_Nor')
hold on
ylabel('PJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text, 'location', 'northwest')
title('BECCS - Bioelectricity')

xline(0)

xline(emf_Norway_elec_prod, 'b','LineWidth', line_width)
xline(emf_EU27_elec_prod, 'C','LineWidth', line_width)
%xline(emf_natural_gas_elec, 'color','#7E2F8E')

%natural gas
x_points = [emf_natural_gas_elec(1) emf_natural_gas_elec(1) emf_natural_gas_elec(2) emf_natural_gas_elec(2)];
y_points= ([0 1.6 1.6 0]);

a=fill(x_points, y_points, color_ng);
a.FaceAlpha = transparancy;

filename = [prefix_output 'Norway_bioelectricity_ccs.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_bioelectricity_ccs_Norway.mat', 'fe_elec_ccs_potential_Nor', 'x_vec_emission_factor');

%% FT diesel NOR
% FT diesel
figure 
area(x_vec_emission_factor, 10^-6*fe_ft_potential_Nor')
hold on
ylabel('PJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text, 'location', 'northwest')
title('Fischer-Tropsch diesel')

xline(0)

%xline(emf_petrol_fuel, 'color',color_diesel,'LineWidth', line_width)
xline(emf_diesel_fuel, 'color',color_diesel,'LineWidth', line_width)

filename = [prefix_output 'Norway_ft_diesel.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_ft_diesel_Norway.mat', 'fe_ft_potential_Nor', 'x_vec_emission_factor');

%% FT diesel with CCS NOR
% FT diesel
figure 
area(x_vec_emission_factor, 10^-6*fe_ft_ccs_potential_Nor')
hold on
ylabel('PJ year^{-1}')
xlabel('kg CO_{2}eq GJ^{-1}')
legend(legend_text,'location', 'northwest')
title('BECCS - Fischer-Tropsch diesel')


xline(0)
%xline(emf_petrol_fuel, 'color', color_diesel,'LineWidth', line_width)
xline(emf_diesel_fuel,'color', color_diesel,'LineWidth', line_width)

filename = [prefix_output 'Norway_ft_diesel_ccs.pdf'];
print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_ft_diesel_ccs_Norway.mat', 'fe_ft_ccs_potential_Nor', 'x_vec_emission_factor');


