function plot_carbon_fluxes(obj)


filename = 'Output/carbon_fluxes.pdf';
if exist(filename, 'file')
    delete(filename)
end

fe_elec_ccs_tot = obj.fe_elec_ccs_tot;
fe_ft_ccs_tot = obj.fe_ft_ccs_tot;


cf_continued_regrowth = obj.carbon_flux_30y_avg_continued_natural_regr_Mg_CO2eq_pr_yr_tot;
cf_elec_ccs = obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_tot;
cf_ft_ccs = obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_tot;

emf_solar_wind_elec = [2 16];
emf_fossil_with_ccs_elec = [44 73];
emf_natural_gas_elec = [136 146];
emf_coal_elec = [220 259];
emf_petrol_fuel = 92.4;
emf_diesel_fuel = 93.9;

avoided_em_solar_wind_elec = -10^-3*mean(emf_solar_wind_elec)*fe_elec_ccs_tot;
avoided_em_fossil_ccs_elec = -10^-3*mean(emf_fossil_with_ccs_elec)*fe_elec_ccs_tot;
avoided_em_natural_gas_elec = -10^-3*mean(emf_natural_gas_elec)*fe_elec_ccs_tot;
avoided_em_coal_elec = -10^-3*mean(emf_coal_elec)*fe_elec_ccs_tot;

avoided_em_diesel = -10^-3*emf_diesel_fuel*fe_ft_ccs_tot;


plotMatrix = zeros(8,8);
plotMatrix(1,1) = cf_continued_regrowth;
plotMatrix(2,2) = cf_elec_ccs;
plotMatrix(3,3) = cf_ft_ccs;

plotMatrix(4,2) = cf_elec_ccs;
plotMatrix(4,4) = avoided_em_solar_wind_elec;

plotMatrix(5,2) = cf_elec_ccs;
plotMatrix(5,5) = avoided_em_fossil_ccs_elec;

plotMatrix(6,2) = cf_elec_ccs;
plotMatrix(6,6) = avoided_em_natural_gas_elec;

plotMatrix(7,2) = cf_elec_ccs;
plotMatrix(7,7) = avoided_em_coal_elec;

plotMatrix(8,3) = cf_ft_ccs;
plotMatrix(8,8) = avoided_em_diesel;

plotMatrix = plotMatrix*10^-9;

figure
barh(plotMatrix, 'stacked');
hold on
xlabel('GtCO2eq yr^{-1}');
yticklabels({'NR', 'BIO-EL', 'BIO-FT', 'BIO-EL AV-REN', 'BIO-EL AV-FOSwCCS', 'BIO-EL AV-NG', 'BIO-EL AV-COAL', 'BIO-FT AV-D'});
legend_labels = {'Natural regrowth', 'Bioelectricity', 'Biofuel FT diesel', 'Avoided emissions PV/wind', 'Avoided emissions fossil w/CCS',...
    'Avoided emissions natural gas', 'Avoided emissions coal', 'Avoided emissions fossil diesel'};
legend(legend_labels, 'Location','southoutside');

print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_carbon_fluxes.mat', 'plotMatrix');

%% SECOND DESIGN
plotMatrix_2 = zeros(8,4);
plotMatrix_2(1,1) = cf_continued_regrowth;


plotMatrix_2(2,3) = cf_ft_ccs;

plotMatrix_2(3,3) = cf_ft_ccs;
plotMatrix_2(3,4) = avoided_em_diesel;

plotMatrix_2(4,2) = cf_elec_ccs;

plotMatrix_2(5,2) = cf_elec_ccs;
plotMatrix_2(5,4) = avoided_em_solar_wind_elec;

plotMatrix_2(6,2) = cf_elec_ccs;
plotMatrix_2(6,4) = avoided_em_fossil_ccs_elec;

plotMatrix_2(7,2) = cf_elec_ccs;
plotMatrix_2(7,4) = avoided_em_natural_gas_elec;

plotMatrix_2(8,2) = cf_elec_ccs;
plotMatrix_2(8,4) = avoided_em_coal_elec;

plotMatrix_2 = plotMatrix_2*10^-9;

filename = 'Output/carbon_fluxes_v2.pdf';
if exist(filename, 'file')
    delete(filename)
end

figure
h = barh(plotMatrix_2, 'stacked');
xlabel('GtCO2eq yr^{-1}');
yticklabels({'NATURAL REGROWTH', 'BIO-FT', 'BIO-FT AV-D','BIO-EL', 'BIO-EL AV-REN', 'BIO-EL AV-FOSwCCS', 'BIO-EL AV-NG', 'BIO-EL AV-COAL' });
legend_labels = {'Natural regrowth', 'Bioelectricity w/CCS', 'Biofuel FT diesel w/CCS', 'Avoided emissions'};
legend(legend_labels, 'Location','southwest');

h(4).LineStyle = '--';
h(4).LineWidth = 1;

h(1).FaceColor = [.2 .6 .5];
h(3).FaceColor = 'b';
h(2).FaceColor = 'y';
h(4).FaceColor = 'm';

print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_carbon_fluxes_v2.mat', 'plotMatrix');
end

