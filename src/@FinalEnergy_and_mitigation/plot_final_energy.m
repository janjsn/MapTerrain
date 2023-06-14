function plot_final_energy(obj)
%PLOT_FINAL_ENERGY Summary of this function goes here
%   Detailed explanation goes here

filename = 'Output/final_energy_global.pdf';


if exist(filename, 'file')
    delete(filename)
end

figure
h = barh([10^-9*obj.fe_ft_ccs_tot 10^-9*obj.fe_elec_ccs_tot] );
xlabel('EJ yr^{-1}');

yticklabels({'BIO-FT','BIO-EL'});

print('-vector','-dpdf', '-r1000', filename)
%save('Output/src_data_carbon_fluxes_v2.mat', 'plotMatrix');
end

