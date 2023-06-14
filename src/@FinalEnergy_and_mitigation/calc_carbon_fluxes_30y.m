function obj= calc_carbon_fluxes_30y(obj)

mSize = size(obj.land_after_abandonment_year);

% Preallocation
obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year = zeros(mSize(1),mSize(2),mSize(3));
obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year = zeros(mSize(1),mSize(2),mSize(3));
obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year = zeros(mSize(1),mSize(2),mSize(3));
obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year = zeros(mSize(1),mSize(2),mSize(3));

% Calcs

%Note that the variable obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year
%is on a basis per hectare. The variable name is misleading.

% Also note that there is no emission factors given as NaN with bioenergy
% productivity, but there are some cases where the combination of positive natural regrowth 
% and no bioenergy productivity are given as NaNs.


obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_electricity.*obj.fe_elec_after_aban_year)-(obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year.*obj.land_after_abandonment_year);
obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_FT.*obj.fe_ft_after_aban_year)-(obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year.*obj.land_after_abandonment_year);
obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_electricity_CCS.*obj.fe_elec_ccs_after_aban_year)-(obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year.*obj.land_after_abandonment_year);
obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_FT_CCS.*obj.fe_ft_ccs_after_aban_year)-(obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year.*obj.land_after_abandonment_year);


obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_electricity.*obj.fe_elec_after_aban_year);
obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_FT.*obj.fe_ft_after_aban_year);
obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_electricity_CCS.*obj.fe_elec_ccs_after_aban_year);
obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year = 10^-3*(obj.emission_factor_FT_CCS.*obj.fe_ft_ccs_after_aban_year);

% Totals
obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year(~isnan(obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year)))));
obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year(~isnan(obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year)))));
obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year(~isnan(obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year)))));
obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year(~isnan(obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year)))));

obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year(~isnan(obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year)))));
obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year(~isnan(obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year)))));
obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year(~isnan(obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year)))));
obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_tot = sum(sum(sum(obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year(~isnan(obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year)))));


%% 2D results at 1 degree

% Preallocation

n_lat_1deg = 180;
n_lon_1deg = 360;

obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);
obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);
obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);
obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);

obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);
obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);
obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_yr_1deg2d = zeros(n_lon_1deg,n_lat_1deg);
obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_1_deg_2d = zeros(n_lon_1deg,n_lat_1deg);

% cf_elec_2d = zeros(mSize(1),mSize(2));
% cf_elec_ccs_2d = zeros(mSize(1),mSize(2));
% cf_ft_2d = zeros(mSize(1),mSize(2));
% cf_ft_ccs_2d = zeros(mSize(1),mSize(2));
% 
% delta_cf_elec_2d = zeros(mSize(1),mSize(2));
% delta_cf_elec_ccs_2d = zeros(mSize(1),mSize(2));
% delta_ft_2d = zeros(mSize(1),mSize(2));
% delta_ft_ccs_2d =zeros(mSize(1),mSize(2));


binary_not_nan_cf_elec = ~isnan(obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year);
binary_not_nan_cf_ft = ~isnan(obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year);
binary_not_nan_cf_elec_ccs = ~isnan(obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year);
binary_not_nan_cf_ft_ccs = ~isnan(obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year);

cf_elec = obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year;
cf_elec(~binary_not_nan_cf_elec) = 0;
cf_ft = obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year;
cf_ft(~binary_not_nan_cf_ft) = 0;
cf_elec_ccs = obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year;
cf_elec_ccs(~binary_not_nan_cf_elec_ccs) = 0;
cf_ft_ccs = obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year;
cf_ft_ccs(~binary_not_nan_cf_ft_ccs) = 0;

delta_cf_elec = obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year;
delta_cf_elec(isnan(delta_cf_elec)) = 0;
delta_cf_ft = obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year;
delta_cf_ft(isnan(delta_cf_ft)) = 0;
delta_cf_elec_ccs = obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year;
delta_cf_elec_ccs(isnan(delta_cf_elec_ccs)) = 0;
delta_cf_ft_ccs = obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year;
delta_cf_ft_ccs(isnan(delta_cf_ft_ccs)) = 0;

cf_elec_2d = sum(cf_elec,3);
cf_ft_2d = sum(cf_ft,3);
cf_elec_ccs_2d = sum(cf_elec_ccs,3);
cf_ft_ccs_2d = sum(cf_ft_ccs,3);

delta_cf_elec_2d = sum(delta_cf_elec,3);
delta_cf_ft_2d = sum(delta_cf_ft,3);
delta_cf_elec_ccs_2d = sum(delta_cf_elec_ccs,3);
delta_cf_ft_ccs_2d = sum(delta_cf_ft_ccs,3);



step = 360/4320;
lon_5arcmin = [-180+step/2:step:180-step/2];
lat_5arcmin = [90-step/2:-step:-90+step/2];

[obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_1_deg_2d, obj.lat_1deg, obj.lon_1deg, ~, ~]  = aggregateMatrix2givenDimensions( cf_elec_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( cf_ft_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( cf_elec_ccs_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( cf_ft_ccs_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

[obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( delta_cf_elec_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( delta_cf_ft_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_yr_1deg2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( delta_cf_elec_ccs_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( delta_cf_ft_ccs_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

abc_co2eq_5arcmin = sum(obj.current_aboveground_carbon_after_abandonment_year,3)*3.67;
obj.current_aboveground_carbon_MgCO2eq_stock_2d = abc_co2eq_5arcmin;
obj.current_aboveground_carbon_MgCO2eq_stock_2d_tot = sum(sum(abc_co2eq_5arcmin));
[obj.current_aboveground_carbon_MgCO2eq_stock_1deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( abc_co2eq_5arcmin,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

land_be_prod = obj.land_availability;
land_be_prod(obj.pe_yield <= 0) = 0;
[obj.land_availability_bioenergy_productive_1deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( land_be_prod,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
obj.land_availability_bioenergy_productive_2d = land_be_prod;
obj.land_availability_bioenergy_productive_tot = sum(sum(land_be_prod));

nr_rate_MgCO2_per_year = obj.mean_natural_regrowth_seq_rate_Mg_CO2_per_year.*obj.land_after_abandonment_year.*(obj.fe_bioelectricity_with_css_GJ_per_year > 0);
nr_rate_MgCO2_per_year_2d = sum(nr_rate_MgCO2_per_year,3);
obj.carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_yr_tot_2d = nr_rate_MgCO2_per_year_2d;

[obj.natural_regrowth_carbon_flux_Mg_CO2eq_per_year_1_deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( nr_rate_MgCO2_per_year_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

obj.natural_regrowth_carbon_flux_Mg_CO2eq_per_ha_year_1_deg_2d = obj.natural_regrowth_carbon_flux_Mg_CO2eq_per_year_1_deg_2d./obj.land_availability_bioenergy_productive_1deg_2d;
obj.natural_regrowth_carbon_flux_Mg_CO2eq_per_ha_year_1_deg_2d(obj.land_availability_bioenergy_productive_1deg_2d == 0) = 0;

obj.carbon_flux_30y_avg_continued_natural_regrowth_Mg_CO2eq_pr_yr = -nr_rate_MgCO2_per_year;
obj.carbon_flux_30y_avg_continued_natural_regrowth_Mg_CO2eq_pr_yr(obj.fe_elec_after_aban_year <= 0) = 0;
obj.carbon_flux_30y_avg_continued_natural_regr_Mg_CO2eq_pr_yr_tot = -sum(nr_rate_MgCO2_per_year(obj.fe_elec_after_aban_year > 0));

obj.carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_ha_yr_tot_2d = obj.carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_yr_tot_2d./obj.land_availability_bioenergy_productive_tot;
obj.carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_ha_yr_tot_2d(obj.land_availability_bioenergy_productive_tot == 0) = 0;

% Final energy totals
obj.fe_elec_tot = sum(sum(sum(obj.fe_elec_after_aban_year(obj.fe_elec_after_aban_year >0))));
obj.fe_elec_ccs_tot = sum(sum(sum(obj.fe_elec_ccs_after_aban_year(obj.fe_elec_ccs_after_aban_year >0))));
obj.fe_ft_tot = sum(sum(sum(obj.fe_ft_after_aban_year(obj.fe_ft_after_aban_year >0))));
obj.fe_ft_ccs_tot = sum(sum(sum(obj.fe_ft_ccs_after_aban_year(obj.fe_ft_ccs_after_aban_year >0))));

% Save 2d matricies at 5arcmin for later Norway calcs
obj.carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_2d = cf_elec_2d;
obj.carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_2d = cf_ft_2d;
obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_2d = cf_elec_ccs_2d;
obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_2d = cf_ft_ccs_2d;

obj.delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_2d = delta_cf_elec_2d;
obj.delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_2d = delta_cf_ft_2d;
obj.delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_yr_2d = delta_cf_elec_ccs_2d;
obj.delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_2d = delta_cf_ft_ccs_2d;

% Save dimension vectors
obj.lat_5arcmin = lat_5arcmin;
obj.lon_5arcmin = lon_5arcmin;

% Primary energy
[obj.pe_1deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( obj.pe,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
obj.pe_yield_1deg_2d = obj.pe_1deg_2d./obj.land_availability_bioenergy_productive_1deg_2d;
obj.pe_yield_1deg_2d(isnan(obj.pe_yield_1deg_2d)) = 0;

% Final energy
fe_elec_ccs_3d = obj.fe_elec_ccs_after_aban_year;
fe_elec_ccs_3d(isnan(fe_elec_ccs_3d)) = 0;
fe_elec_ccs_2d = sum(fe_elec_ccs_3d,3);
obj.fe_bioelectricity_ccs_2d = fe_elec_ccs_2d;
[obj.fe_bioelectricity_ccs_1deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( fe_elec_ccs_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

fe_ft_ccs_3d = obj.fe_ft_after_aban_year; 
fe_ft_ccs_3d(isnan(fe_ft_ccs_3d)) = 0;
fe_ft_ccs_2d = sum(fe_ft_ccs_3d,3);
obj.fe_ft_ccs_2d = fe_ft_ccs_2d;
[obj.fe_ft_ccs_1deg_2d, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( fe_ft_ccs_2d,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

%Final energy yields
obj.fe_yield_bioelectricity_ccs_1deg_2d = obj.fe_bioelectricity_ccs_1deg_2d./obj.land_availability_bioenergy_productive_1deg_2d;
obj.fe_yield_bioelectricity_ccs_1deg_2d(obj.land_availability_bioenergy_productive_1deg_2d == 0) = 0;

obj.fe_yield_ft_ccs_1deg_2d= obj.fe_ft_ccs_1deg_2d./obj.land_availability_bioenergy_productive_1deg_2d;
obj.fe_yield_ft_ccs_1deg_2d(obj.land_availability_bioenergy_productive_1deg_2d == 0) = 0;

%% Emission factors, no LUC

kgco2eq_emissions_no_luc_beccs_bioelectricity_5arcmin = obj.fe_bioelectricity_ccs_2d.*obj.emission_factor_electricity_CCS_excluding_luc;
kgco2eq_emissions_no_luc_beccs_ft_5arcmin = obj.fe_ft_ccs_2d.*obj.emission_factor_FT_CCS_excluding_luc;


[kgco2eq_emissions_no_luc_beccs_bioelectricity_1deg, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( kgco2eq_emissions_no_luc_beccs_bioelectricity_5arcmin,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );
[kgco2eq_emissions_no_luc_beccs_ft_1deg, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( kgco2eq_emissions_no_luc_beccs_ft_5arcmin,lon_5arcmin,lat_5arcmin, n_lon_1deg, n_lat_1deg );

obj.emission_factor_electricity_CCS_excluding_luc_1_deg_2d = kgco2eq_emissions_no_luc_beccs_bioelectricity_1deg./obj.fe_bioelectricity_ccs_1deg_2d;
obj.emission_factor_electricity_CCS_excluding_luc_1_deg_2d(obj.fe_bioelectricity_ccs_1deg_2d == 0) = 0;

obj.emission_factor_FT_CCS_excluding_luc_1_deg_2d = kgco2eq_emissions_no_luc_beccs_ft_1deg./obj.fe_ft_ccs_1deg_2d;
obj.emission_factor_FT_CCS_excluding_luc_1_deg_2d(obj.fe_ft_ccs_1deg_2d == 0) = 0;



end

