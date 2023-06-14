function export_Norway(obj, Norway_mask)


filename = 'Output/Norway_wTrondelag.nc';

if exist(filename,'file')
    delete(filename)
end

step = 360/4320;
lon = [-180+step/2:step:180-step/2];
lat = [90-step/2:-step:-90+step/2];

dummy_matrix = zeros(length(lon),length(lat));
[~, ~, ~, lat_bnds,  ~]  = aggregateMatrix2givenDimensions( dummy_matrix,lon,lat, length(lon), length(lat) );
[ area_after_lat ] = get_cell_area_per_latitude( lat_bnds, abs(lon(1)-lon(2)));


abandoned_cropland_fractions = zeros(length(lon),length(lat));
for lats = 1:length(lat)
    abandoned_cropland_fractions(:,lats) = obj.land_availability(:,lats)/area_after_lat(lats);
end

abandoned_cropland_fractions(~Norway_mask) = -999;
abandoned_cropland_hectare = obj.land_availability;
abandoned_cropland_hectare(~Norway_mask) = -999;

pe_Norway = obj.pe;
fprintf(['Norway primary energy potential: ' num2str(10^-6*sum(sum(pe_Norway(Norway_mask)))) '\n']);
pe_tot_Norway  = 10^-6*sum(sum(pe_Norway(Norway_mask)));
pe_Norway(~Norway_mask) = -999;

fe_el_Norway = obj.fe_bioelectricity;
fprintf(['Norway bioelectricity potential: ' num2str(10^-6*sum(sum(fe_el_Norway(Norway_mask)))) '\n']);
fe_el_tot_Norway =10^-6*sum(sum(fe_el_Norway(Norway_mask)));
fe_el_Norway(~Norway_mask) = -999;

fe_el_ccs_Norway = obj.fe_bioelectricity_with_css_GJ_per_year;
fprintf(['Norway bioelectricity potential with CCS: ' num2str(10^-6*sum(sum(fe_el_ccs_Norway(Norway_mask)))) '\n']);
fe_el_ccs_tot_Norway = 10^-6*sum(sum(fe_el_ccs_Norway(Norway_mask)));
fe_el_ccs_Norway(~Norway_mask) = -999;

fe_ft_Norway = obj.fe_FT_diesel;
fprintf(['Norway FT potential: ' num2str(10^-6*sum(sum(fe_ft_Norway(Norway_mask)))) '\n']);
fe_ft_tot_Norway = 10^-6*sum(sum(fe_ft_Norway(Norway_mask)));
fe_ft_Norway(~Norway_mask) = -999;

fe_ft_ccs_Norway = obj.fe_FT_diesel_with_ccs_GJ_per_year;
fprintf(['Norway FT potential: ' num2str(10^-6*sum(sum(fe_ft_ccs_Norway(Norway_mask)))) '\n']);
fe_ft_ccs_tot_Norway = 10^-6*sum(sum(fe_ft_ccs_Norway(Norway_mask)));
fe_ft_ccs_Norway(~Norway_mask) = -999;

figure
barh([ fe_ft_ccs_tot_Norway fe_el_ccs_tot_Norway]);

fossil_fuel_road_transport_demand_Norway = (25.9+6.7)*3.6;
current_biofuel_use_Norway = 4.2*3.6;
ambition_gap_target_2030 = ((4.5+1)*3.6)-current_biofuel_use_Norway;


%REF SSB https://www.ssb.no/transport-og-reiseliv/landtransport/artikler/utfordringer-med-fornybart-drivstoff

hold on
%plot([2.5 4.5], [fossil_fuel_road_transport_demand_Norway*0.05 fossil_fuel_road_transport_demand_Norway*0.05], '--', 'LineWidth', 2);
%plot([2.5 4.5], [0.25*current_biofuel_use_Norway 0.25*current_biofuel_use_Norway], '--', 'LineWidth', 2);
%plot([2.5 4.5], [ambition_gap_target_2030 ambition_gap_target_2030], '--', 'LineWidth', 2);

% scatter(3, fossil_fuel_road_transport_demand_Norway*0.05);
% scatter(3,0.25*current_biofuel_use_Norway);
% scatter(3,ambition_gap_target_2030)


yticklabels({ 'BIO-FT','BIO-EL',});
xlabel('PJ year^{-1}');
%legend({'Final energy', '5% of fossil diesel supply', '25% of current biofuel supply', 'Ambition gap 2030'}, 'Location', 'Northwest');
%legend({'Final energy', '25% of current biofuel supply', 'Ambition gap 2030'}, 'Location', 'Northwest');

filename_plot = [ 'Norway_fe_potential_v2.pdf'];
if exist(filename_plot, 'file')
    delete(filename_plot)
end

print('-vector','-dpdf', '-r1000', filename_plot)

%% PLOT CARBON FLUXES
% REF Scarlat et al. (2022), Applied energy
emf_Norway_elec_prod = 7.8;
emf_EU27_elec_prod = 86.1;
% REF HANSSEN ET AL. (2020)
emf_natural_gas_elec = [136 146];
emf_diesel_fuel = 93.9;

%Carbon fluxes
cf_el_ccs = sum(sum(obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_2d(Norway_mask)));
cf_ft_ccs = sum(sum(obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_2d(Norway_mask)));
cf_natural_regrowth = -sum(sum(obj.carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_yr_tot_2d(Norway_mask)));

%Avoided emissions
avoided_em_Norway_el = -10^3*sum(sum(fe_el_ccs_tot_Norway*emf_Norway_elec_prod));
avoided_em_EU27_el = -10^3*sum(sum(fe_el_ccs_tot_Norway*emf_EU27_elec_prod));
avoided_em_natural_gas_el = -10^3*sum(sum(fe_el_ccs_tot_Norway*mean(emf_natural_gas_elec)));
avoided_em_diesel_fuel = -10^3*sum(sum(fe_ft_ccs_tot_Norway*emf_diesel_fuel));

plotMatrix = zeros(7,4);
plotMatrix(1,1) = cf_natural_regrowth;

plotMatrix(2,3) = cf_ft_ccs;

plotMatrix(3,3) = cf_ft_ccs;
plotMatrix(3,4) = avoided_em_diesel_fuel;

plotMatrix(4,2) = cf_el_ccs;

plotMatrix(5,2) = cf_el_ccs;
plotMatrix(5,4) = avoided_em_Norway_el;

plotMatrix(6,2) = cf_el_ccs;
plotMatrix(6,4) = avoided_em_EU27_el;

plotMatrix(7,2) = cf_el_ccs;
plotMatrix(7,4) = avoided_em_natural_gas_el;

plotMatrix = 10^-6*plotMatrix;

filename = 'Output/carbon_fluxes_NOR_v2.pdf';
if exist(filename, 'file')
    delete(filename)
end

figure
h = barh(plotMatrix, 'stacked');
xlabel('MtCO2eq yr^{-1}');
yticklabels({'NATURAL REGROWTH', 'BIO-FT', 'BIO-FT AV-D','BIO-EL', 'BIO-EL AV-NOR', 'BIO-EL AV-EU27', 'BIO-EL AV-NG'});
legend_labels = {'Natural regrowth', 'Bioelectricity w/CCS', 'Biofuel FT diesel w/CCS', 'Avoided emissions'};
legend(legend_labels, 'Location','southwest');

h(4).LineStyle = '--';
h(4).LineWidth = 1;

h(1).FaceColor = [.2 .6 .5];
h(3).FaceColor = 'b';
h(2).FaceColor = 'y';
h(4).FaceColor = 'm';
xlim([-1 0]);

print('-vector','-dpdf', '-r1000', filename)
save('Output/src_data_carbon_fluxes_NOR_v2.mat', 'plotMatrix');


%% EXPORT


nccreate(filename,'lat','Dimensions',{'lat' length(lat)});
ncwriteatt(filename, 'lat', 'standard_name', 'latitude');
ncwriteatt(filename, 'lat', 'long_name', 'latitude');
ncwriteatt(filename, 'lat', 'units', 'degrees_north');
ncwriteatt(filename, 'lat', '_CoordinateAxisType', 'Lat');
ncwriteatt(filename, 'lat', 'axis', 'Y');

nccreate(filename,'lon','Dimensions',{'lon' length(lon)});
ncwriteatt(filename, 'lon', 'standard_name', 'longitude');
ncwriteatt(filename, 'lon', 'long_name', 'longitude');
ncwriteatt(filename, 'lon', 'units', 'degrees_east');
ncwriteatt(filename, 'lon', '_CoordinateAxisType', 'Lon');
ncwriteatt(filename, 'lon', 'axis', 'X');

nccreate(filename, 'abandoned_cropland_hectare', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_hectare', 'standard_name', 'abandoned_cropland_hectare');
ncwriteatt(filename, 'abandoned_cropland_hectare', 'long_name', 'abandoned_cropland_hectare');
ncwriteatt(filename, 'abandoned_cropland_hectare', 'units', 'ha');
ncwriteatt(filename, 'abandoned_cropland_hectare', 'missing_value', '-999');

nccreate(filename, 'abandoned_cropland_fractions', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_fractions', 'standard_name', 'abandoned_cropland_fractions');
ncwriteatt(filename, 'abandoned_cropland_fractions', 'long_name', 'abandoned_cropland_fractions');
ncwriteatt(filename, 'abandoned_cropland_fractions', 'units', 'fractions');
ncwriteatt(filename, 'abandoned_cropland_fractions', 'missing_value', '-999');

nccreate(filename, 'primary_energy', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy', 'standard_name', 'primary_energy');
ncwriteatt(filename, 'primary_energy', 'long_name', 'primary_energy');
ncwriteatt(filename, 'primary_energy', 'units', 'TJ year-1');
ncwriteatt(filename, 'primary_energy', 'missing_value', '-999');

nccreate(filename, 'final_energy_bioelectricity', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_bioelectricity', 'standard_name', 'final_energy_bioelectricity');
ncwriteatt(filename, 'final_energy_bioelectricity', 'long_name', 'final_energy_bioelectricity');
ncwriteatt(filename, 'final_energy_bioelectricity', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_bioelectricity', 'missing_value', '-999');

nccreate(filename, 'final_energy_bioelectricity_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'standard_name', 'final_energy_bioelectricity_ccs');
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'long_name', 'final_energy_bioelectricity_ccs');
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'missing_value', '-999');

nccreate(filename, 'final_energy_fischer_tropsch_diesel', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'standard_name', 'final_energy_fischer_tropsch_diesel');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'long_name', 'final_energy_fischer_tropsch_diesel');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'missing_value', '-999');

nccreate(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'standard_name', 'final_energy_fischer_tropsch_diesel_ccs');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'long_name', 'final_energy_fischer_tropsch_diesel_ccs');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'missing_value', '-999');

% nccreate(filename, 'emission_factor_bioelectricity', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
% ncwriteatt(filename, 'emission_factor_bioelectricity', 'standard_name', 'emission_factor_bioelectricity');
% ncwriteatt(filename, 'emission_factor_bioelectricity', 'long_name', 'emission_factor_bioelectricity');
% ncwriteatt(filename, 'emission_factor_bioelectricity', 'units', 'kgCO2eq GJ-1');
% ncwriteatt(filename, 'emission_factor_bioelectricity', 'missing_value', '-999');
%
% nccreate(filename, 'emission_factor_bioelectricity_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
% ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'standard_name', 'emission_factor_bioelectricity_ccs');
% ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'long_name', 'emission_factor_bioelectricity_ccs');
% ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'units', 'kgCO2eq GJ-1');
% ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'missing_value', '-999');
%
% nccreate(filename, 'emission_factor_fischer_tropsch_diesel', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'standard_name', 'emission_factor_fischer_tropsch_diesel');
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'long_name', 'emission_factor_fischer_tropsch_diesel');
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'units', 'kgCO2eq GJ-1');
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'missing_value', '-999');
%
% nccreate(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'standard_name', 'emission_factor_fischer_tropsch_diesel_ccs');
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'long_name', 'emission_factor_fischer_tropsch_diesel_ccs');
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'units', 'kgCO2eq GJ-1');
% ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'missing_value', '-999');

ncwrite(filename, 'lat', lat);
ncwrite(filename, 'lon', lon);
ncwrite(filename, 'abandoned_cropland_fractions', abandoned_cropland_fractions);
ncwrite(filename, 'primary_energy', pe_Norway);
ncwrite(filename, 'final_energy_bioelectricity', 10^-3*fe_el_Norway);
ncwrite(filename, 'final_energy_bioelectricity_ccs', 10^-3*fe_el_ccs_Norway);
ncwrite(filename, 'final_energy_fischer_tropsch_diesel', 10^-3*fe_ft_Norway);
ncwrite(filename, 'final_energy_fischer_tropsch_diesel_ccs', 10^-3*fe_ft_ccs_Norway);
ncwrite(filename, 'abandoned_cropland_hectare', abandoned_cropland_hectare);

%% TRØNDELAG
ncid = netcdf.open('municipal_ids.nc');
kommune_ids = netcdf.getVar(ncid,2);
fylke_ids = netcdf.getVar(ncid,3);

id_Trondelag = 50;

unique_kommune = unique(kommune_ids(kommune_ids > 0));
unique_fylke = unique(fylke_ids(fylke_ids > 0));

land_kommune = zeros(4320,2160);
land_fylke = zeros(4320,2160);

pe_kommune = zeros(4320,2160);
pe_fylke = zeros(4320,2160);

land_kommune_tot = zeros(1,length(unique_kommune));
pe_kommune_tot = zeros(1,length(unique_kommune));

land_fylke_tot = zeros(1,length(unique_fylke));
pe_fylke_tot = zeros(1,length(unique_fylke));

[~,~,lockup_table_municipalities] = xlsread('Data/lockup_table.xlsx', 'Kommune');
[~,~,lockup_table_fylke] = xlsread('Data/lockup_table.xlsx', 'Fylke');

mSize_mun = size(lockup_table_municipalities);
mSize_fylke= size(lockup_table_fylke);

out_municipal = lockup_table_municipalities;
out_fylke = lockup_table_fylke;
out_municipal{1,3} = 'Abandoned cropland (ha)';
out_municipal{1,4} = 'Bioenergy potential (TJ yr-1)';


for i = 1:length(unique_kommune)

    mask_this = kommune_ids == unique_kommune(i);
    land_kommune_tot(i) = sum(sum(abandoned_cropland_hectare(mask_this)));
    pe_kommune_tot(i) = sum(sum(pe_Norway(mask_this)));

    land_kommune(mask_this) = land_kommune_tot(i);
    pe_kommune(mask_this) = pe_kommune_tot(i);

    for j = 2:mSize_mun(1)
        if unique_kommune(i) == lockup_table_municipalities{j,2}
            out_municipal{j,3} = land_kommune_tot(i);
            out_municipal{j,4} = pe_kommune_tot(i)*10^-3;
        end

    end

end

for i = 1:length(unique_fylke)
    mask_this = fylke_ids == unique_fylke(i);
    land_fylke_tot(i) = sum(sum(abandoned_cropland_hectare(mask_this)));
    pe_fylke_tot(i) = sum(sum(pe_Norway(mask_this)));

    pe_fylke(mask_this) = pe_fylke_tot(i);
    land_fylke(mask_this) = land_fylke_tot(i);

    for j = 2:mSize_fylke(1)
        if unique_fylke(i) == lockup_table_fylke{j,2}
            out_fylke{j,3} = land_fylke_tot(i);
            out_fylke{j,4} = pe_fylke_tot(i)*10^-6;
        end

    end
end

pe_tot_Trondelag = sum(sum(pe_Norway(fylke_ids == id_Trondelag)))
land_Trondelag = sum(sum(abandoned_cropland_hectare(fylke_ids == id_Trondelag)))

nat_reg_Trondelag = sum(sum(obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year.*land_Trondelag))

land_kommune_Trondelag = land_kommune;
pe_kommune_Trondelag = pe_kommune;

land_kommune_Trondelag(fylke_ids ~= id_Trondelag) = -999;
pe_kommune_Trondelag(fylke_ids ~= id_Trondelag) = -999;


pe_Trondelag = 10^-3*pe_Norway;

pe_Trondelag(fylke_ids ~= id_Trondelag) = -999;

ac_Trondelag = abandoned_cropland_hectare;
ac_Trondelag(fylke_ids ~= id_Trondelag) = -999;


nccreate(filename, 'abandoned_cropland_fylke', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_fylke', 'standard_name', 'abandoned_cropland_fylke');
ncwriteatt(filename, 'abandoned_cropland_fylke', 'long_name', 'abandoned_cropland_fylke');
ncwriteatt(filename, 'abandoned_cropland_fylke', 'units', 'ha');
ncwriteatt(filename, 'abandoned_cropland_fylke', 'missing_value', '-999');

nccreate(filename, 'abandoned_cropland_municipalities', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_municipalities', 'standard_name', 'abandoned_cropland_municipalities');
ncwriteatt(filename, 'abandoned_cropland_municipalities', 'long_name', 'abandoned_cropland_municipalities');
ncwriteatt(filename, 'abandoned_cropland_municipalities', 'units', 'ha');
ncwriteatt(filename, 'abandoned_cropland_municipalities', 'missing_value', '-999');

nccreate(filename, 'primary_energy_fylke', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy_fylke', 'standard_name', 'primary_energy_fylke');
ncwriteatt(filename, 'primary_energy_fylke', 'long_name', 'primary_energy_fylke');
ncwriteatt(filename, 'primary_energy_fylke', 'units', 'PJ yr-1');
ncwriteatt(filename, 'primary_energy_fylke', 'missing_value', '-999');

nccreate(filename, 'primary_energy_municipalities', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy_municipalities', 'standard_name', 'primary_energy_municipalities');
ncwriteatt(filename, 'primary_energy_municipalities', 'long_name', 'primary_energy_municipalities');
ncwriteatt(filename, 'primary_energy_municipalities', 'units', 'TJ yr-1');
ncwriteatt(filename, 'primary_energy_municipalities', 'missing_value', '-999');

nccreate(filename, 'abandoned_cropland_Trondelag_municipalities', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_Trondelag_municipalities', 'standard_name', 'abandoned_cropland_Trondelag_municipalities');
ncwriteatt(filename, 'abandoned_cropland_Trondelag_municipalities', 'long_name', 'abandoned_cropland_Trondelag_municipalities');
ncwriteatt(filename, 'abandoned_cropland_Trondelag_municipalities', 'units', 'ha');
ncwriteatt(filename, 'abandoned_cropland_Trondelag_municipalities', 'missing_value', '-999');

nccreate(filename, 'primary_energy_Trondelag_municipalities', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy_Trondelag_municipalities', 'standard_name', 'primary_energy_Trondelag_municipalities');
ncwriteatt(filename, 'primary_energy_Trondelag_municipalities', 'long_name', 'primary_energy_Trondelag_municipalities');
ncwriteatt(filename, 'primary_energy_Trondelag_municipalities', 'units', 'TJ yr-1');
ncwriteatt(filename, 'primary_energy_Trondelag_municipalities', 'missing_value', '-999');

nccreate(filename, 'abandoned_cropland_Trondelag', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_Trondelag', 'standard_name', 'abandoned_cropland_Trondelag');
ncwriteatt(filename, 'abandoned_cropland_Trondelag', 'long_name', 'abandoned_cropland_Trondelag');
ncwriteatt(filename, 'abandoned_cropland_Trondelag', 'units', 'ha');
ncwriteatt(filename, 'abandoned_cropland_Trondelag', 'missing_value', '-999');

nccreate(filename, 'primary_energy_Trondelag', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy_Trondelag', 'standard_name', 'primary_energy_Trondelag');
ncwriteatt(filename, 'primary_energy_Trondelag', 'long_name', 'primary_energy_Trondelag');
ncwriteatt(filename, 'primary_energy_Trondelag', 'units', 'TJ yr-1');
ncwriteatt(filename, 'primary_energy_Trondelag', 'missing_value', '-999');

pe_kommune = pe_kommune*10^-3;
pe_fylke = pe_fylke*10^-6;

land_fylke(~Norway_mask) = -999;
land_kommune(~Norway_mask) = -999;
pe_fylke(~Norway_mask) = -999;
pe_kommune(~Norway_mask) = -999;

ncwrite(filename, 'abandoned_cropland_fylke', land_fylke);
ncwrite(filename, 'abandoned_cropland_municipalities', land_kommune);
ncwrite(filename, 'primary_energy_fylke', pe_fylke);
ncwrite(filename, 'primary_energy_municipalities', pe_kommune);
ncwrite(filename, 'abandoned_cropland_Trondelag_municipalities', land_kommune_Trondelag);
ncwrite(filename, 'primary_energy_Trondelag_municipalities', pe_kommune_Trondelag);
ncwrite(filename, 'abandoned_cropland_Trondelag', ac_Trondelag);
ncwrite(filename, 'primary_energy_Trondelag', pe_Trondelag);

writecell(out_municipal,'Norway.xlsx', 'Sheet', 'Municipalities');
writecell(out_fylke,'Norway.xlsx', 'Sheet', 'Regions');


end

