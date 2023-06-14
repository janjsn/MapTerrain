diary log_out.txt
diary on
fprintf('----- NEW RUN -------- \n')

tic

addpath(genpath(pwd()));

run_all = 1;
run_gcam = 1;
export = 0;

carbon_to_co2_factor = 3.67;

if run_all == 1
    fprintf('Importing country masks.. \n');
    [ CountryMask ] = importCountryMasks(  );

    fprintf('Importing gridded world population.. \n')
    gpw_population_2020 = readgeoraster('Data/gpw-v4-population-count-rev11_2020_2pt5_min_tif/gpw_v4_population_count_rev11_2020_2pt5_min.tif');
    gpw_population_2020 = gpw_population_2020';
    gpw_population_2020(gpw_population_2020 <= 0) = 0;
    fprintf('Importing GCAM land use projections.. \n')
    GCAM_array = get_GCAM(  );
    GCAM_2020 = GCAM('Data/GCAM/GCAM_Demeter_LU_ssp2_rcp26_modelmean_2020.nc', 2020, 'SSP2-26, 2020', 2, 2.6);
    GCAM_2015 = GCAM('Data/GCAM/GCAM_Demeter_LU_ssp2_rcp26_modelmean_2015.nc', 2015, 'SSP2-26, 2015', 2, 2.6);

    fprintf('Importing abandoned cropland.. \n')
    ncid_ac = netcdf.open('Data/abandoned_cropland_1992_2020_5arcmin_timestamp_2022_07_18_1212.nc');
    ac_ha = netcdf.getVar(ncid_ac,3);
    ac_fractions = netcdf.getVar(ncid_ac,4);
    netcdf.close(ncid_ac);

    fprintf('Importing cropland cover.. \n')
    ncid_crop = netcdf.open('Data/Croplands_Global_2020_5arcmin_timestamp_2022_10_18_1002.nc');
    cropland = netcdf.getVar(ncid_crop,3);
    cropland_fractions = netcdf.getVar(ncid_crop,4);
    netcdf.close(ncid_crop);

    fprintf('Importing bioenergy potentials from GAEZ.. \n')
    ncid_bioenergy_cruts32 = netcdf.open('Data/GAEZ/bioenergy_ac_1992_2020_cruts32_1980_2010_hi_rf.nc');
    pe_cruts32 = netcdf.getVar(ncid_bioenergy_cruts32,13);
    pe_cruts32_opt = netcdf.getVar(ncid_bioenergy_cruts32,14);
    crops_cruts_32_opt = netcdf.getVar(ncid_bioenergy_cruts32,11);
    netcdf.close(ncid_bioenergy_cruts32);

    ncid_bioenergy_noresm_rcp45_2020s = netcdf.open('Data/GAEZ/bioenergy_ac_1992_2020_NorESM-M_RCP45_2011_2040_hi_rf.nc');
    pe_noresm_rcp45_2020s = netcdf.getVar(ncid_bioenergy_noresm_rcp45_2020s,13);
    pe_noresm_rcp45_2020s_opt = netcdf.getVar(ncid_bioenergy_noresm_rcp45_2020s,14);
    crops_noresm_rcp45_2020s_opt = netcdf.getVar(ncid_bioenergy_noresm_rcp45_2020s,11);
    dm_yield_all_noresm_rcp45_2020s = netcdf.getVar(ncid_bioenergy_noresm_rcp45_2020s,8);
    netcdf.close(ncid_bioenergy_noresm_rcp45_2020s);

    fprintf('Importing current aboveground carbon on abandoned cropland... \n')
    ncid_ac_carbon = netcdf.open('Data/standing_carbon_on_abandoned_cropland_2020.nc');
    time_ac_carbon = netcdf.getVar(ncid_ac_carbon,2);
    ac_after_abandonment_year = netcdf.getVar(ncid_ac_carbon,3);
    ac_abc_after_abandonment_year = netcdf.getVar(ncid_ac_carbon,5);
    netcdf.close(ncid_ac_carbon);

    ac_abc_intensity_tonC_per_ha_after_abandonment_year = ac_abc_after_abandonment_year./ac_after_abandonment_year;
    ac_abc_intensity_tonC_per_ha_after_abandonment_year(isnan(ac_abc_intensity_tonC_per_ha_after_abandonment_year)) = 0;

    ncid_ac_carbon_aggregated = netcdf.open('Data/standing_carbon_on_abandoned_cropland_2020.nc');
    ac_abc_grid_tot_2022 = netcdf.getVar(ncid_ac_carbon_aggregated,5);
    ac_abc_regrowth_rate_grid_mean_2022 = netcdf.getVar(ncid_ac_carbon_aggregated,6);
    netcdf.close(ncid_ac_carbon_aggregated);

end

ac_abc_regrowth_rate_grid_mean_2022(ac_abc_regrowth_rate_grid_mean_2022 == -999) = 0;

ac_per_cropland_gridded = ac_ha./cropland;
ac_per_cropland_gridded(cropland == 0) = 0;
ac_per_cropland_gridded(CountryMask.countryMask_5arcmin == CountryMask.waterValue) = NaN;
ac_per_cropland2020_country_averages = zeros(length(CountryMask.lon_5arcmin), length(CountryMask.lat_5arcmin));
ac_per_capita_country_averages = zeros(length(CountryMask.lon_5arcmin), length(CountryMask.lat_5arcmin));




SSPs = [1 2 4 5];
c_to_co2 = 3.67;

%Clean up

if run_gcam == 1

    for countries = 1:length(CountryMask.CountryArray)

        this_id = CountryMask.CountryArray(countries).GPW_country_ISO_numeric;
        this_mask_25 = CountryMask.countryMask == this_id;
        this_mask_5arcmin = CountryMask.countryMask_5arcmin == this_id;
        this_mask_05 = CountryMask.countryMask_05deg == this_id;

        CountryMask.CountryArray(countries).population_2020 = sum(sum(gpw_population_2020(this_mask_25)));
        CountryMask.CountryArray(countries).abandoned_cropland_1992_to_2020_tot = sum(sum(ac_ha(this_mask_5arcmin)));
        CountryMask.CountryArray(countries).cropland_2020_tot = sum(sum(cropland(this_mask_5arcmin)));

        CountryMask.CountryArray(countries).pe_ac_cruts32_tot = sum(sum(pe_cruts32_opt(this_mask_5arcmin)));
        CountryMask.CountryArray(countries).pe_ac_noresm_rcp45_2020s_tot = sum(sum(pe_noresm_rcp45_2020s_opt(this_mask_5arcmin)));

        CountryMask.CountryArray(countries).abandoned_cropland_as_share_of_2020_cropland = CountryMask.CountryArray(countries).abandoned_cropland_1992_to_2020_tot/CountryMask.CountryArray(countries).cropland_2020_tot;

        %CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2050_per_SSP

        % Per capita
        CountryMask.CountryArray(countries).ac_per_cap = CountryMask.CountryArray(countries).abandoned_cropland_1992_to_2020_tot/CountryMask.CountryArray(countries).population_2020;
        CountryMask.CountryArray(countries).cropland_per_cap = CountryMask.CountryArray(countries).cropland_2020_tot/CountryMask.CountryArray(countries).population_2020;

        CountryMask.CountryArray(countries).pe_ac_cruts32_per_cap = CountryMask.CountryArray(countries).pe_ac_cruts32_tot/CountryMask.CountryArray(countries).population_2020;
        CountryMask.CountryArray(countries).pe_ac_noresm_rcp45_2020s_per_cap = CountryMask.CountryArray(countries).pe_ac_noresm_rcp45_2020s_tot/CountryMask.CountryArray(countries).population_2020;

        for projections = 1:length(GCAM_array)
            idx = find(GCAM_array(projections).SSP == SSPs);

            if GCAM_array(projections).RCP ~= 2.6
                continue
            end

            if GCAM_array(projections).year == 2050
                CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2050_per_SSP(idx) = sum(sum(GCAM_array(projections).bioenergy_crops_rf_hectare(this_mask_05)+GCAM_array(projections).bioenergy_crops_ir_hectare(this_mask_05)));
            elseif GCAM_array(projections).year == 2100
                CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2100_per_SSP(idx) = sum(sum(GCAM_array(projections).bioenergy_crops_rf_hectare(this_mask_05)+GCAM_array(projections).bioenergy_crops_ir_hectare(this_mask_05)));

            end
        end

        CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2050_per_cap_per_SSP = CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2050_per_SSP/CountryMask.CountryArray(countries).population_2020;
        CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2100_per_cap_per_SSP = CountryMask.CountryArray(countries).bioenergy_land_use_GCAM_2100_per_SSP/CountryMask.CountryArray(countries).population_2020;

        ac_per_cropland2020_country_averages(this_mask_5arcmin) = CountryMask.CountryArray(countries).abandoned_cropland_as_share_of_2020_cropland;
        ac_per_capita_country_averages(this_mask_5arcmin) = CountryMask.CountryArray(countries).ac_per_cap;

        if strcmp(CountryMask.CountryArray(countries).country_name, 'Norway')
            Norway =CountryMask.CountryArray(countries);

        end

    end


    country_names = {CountryMask.CountryArray.country_name};
    country_ac_per_cap = [CountryMask.CountryArray.ac_per_cap];
    country_ac_per_cropland = [CountryMask.CountryArray.abandoned_cropland_as_share_of_2020_cropland];

    [country_ac_per_cap_sorted,country_ac_per_cap_sorted_index ] = sort(country_ac_per_cap, 'descend');
    country_names_sorted_ac_per_cap = {country_names{country_ac_per_cap_sorted_index}};
    [country_ac_per_cropland2020_sorted,country_ac_per_cropland2020_sorted_index ] = sort(country_ac_per_cropland, 'descend');
    country_names_sorted_ac_per_cropland = {country_names{country_ac_per_cropland2020_sorted_index}};
end

idx_Norway = find(strcmp({CountryMask.CountryArray.country_name},'Norway'));
Norway_iso = CountryMask.CountryArray(idx_Norway).GPW_country_ISO_numeric;
Norway_mask_5arcmin = CountryMask.countryMask_5arcmin == Norway_iso;

%% Calc final energy and mitigation
EnergyConversion_Hansen = EnergyConversion();
CropTypeArray = makeCropTypeArray();
willow_id = 5;
n_crops = 3;

FEM = FinalEnergy_and_mitigation;
FEM.time_period = 30;

FEM.land_availability = ac_ha;
FEM.crop_allocation = crops_noresm_rcp45_2020s_opt;
FEM.binary_is_grassy = (0 < crops_noresm_rcp45_2020s_opt) & (crops_noresm_rcp45_2020s_opt < 4);
FEM.binary_is_woody = crops_noresm_rcp45_2020s_opt == willow_id;
FEM.dm_yield = zeros(4320,2160);
FEM.carbon_yield = zeros(4320,2160);
FEM.pe = pe_noresm_rcp45_2020s_opt; %Set directly from input

FEM.pe_yield = FEM.pe./FEM.land_availability;


FEM.fe_bioelectricity = zeros(4320,2160);
FEM.fe_FT_diesel = zeros(4320,2160);
FEM.fe_FT_diesel_with_ccs_GJ_per_year = zeros(4320,2160);
FEM.fe_bioelectricity_with_css_GJ_per_year = zeros(4320,2160);

for crops = 1:n_crops
    this_id = CropTypeArray(crops).ID;
    binary_this_crop = FEM.crop_allocation == this_id;
    dm_yield_this = dm_yield_all_noresm_rcp45_2020s(:,:,crops);
    FEM.dm_yield(binary_this_crop) = dm_yield_this(binary_this_crop);
    FEM.carbon_yield(binary_this_crop) = dm_yield_this(binary_this_crop)*CropTypeArray(crops).carbon_content_of_dry_mass;
    %Conversion factors
    if this_id ~= willow_id
        lhv_Hansen = EnergyConversion_Hansen.lhv_grassy;
        CropTypeArray(crops).electricity_from_biomass_GJelec_per_ton_dm = EnergyConversion_Hansen.electricity_from_grassy_biomass_GJelec_per_ton_dm*CropTypeArray(crops).calorificValue_weighted_MJperKgYieldModelOutput_lhv/lhv_Hansen;
        CropTypeArray(crops).ft_diesel_from_biomass_GJfuel_per_ton_dm = EnergyConversion_Hansen.ft_diesel_from_grassy_biomass_GJfuel_per_ton_dm*CropTypeArray(crops).calorificValue_weighted_MJperKgYieldModelOutput_lhv/lhv_Hansen;
        CropTypeArray(crops).ethanol_from_biomass_GJfuel_per_ton_dm = EnergyConversion_Hansen.ethanol_from_woody_biomass_GJfuel_per_ton_dm*CropTypeArray(crops).calorificValue_weighted_MJperKgYieldModelOutput_lhv/lhv_Hansen;
    elseif this_id == willow_id
        lhv_Hansen = EnergyConversion_Hansen.lhv_woody;
        CropTypeArray(crops).electricity_from_biomass_GJelec_per_ton_dm = EnergyConversion_Hansen.electricity_from_woody_biomass_GJelec_per_ton_dm*CropTypeArray(crops).calorificValue_weighted_MJperKgYieldModelOutput_lhv/lhv_Hansen;
        CropTypeArray(crops).ft_diesel_from_biomass_GJfuel_per_ton_dm = EnergyConversion_Hansen.ft_diesel_from_woody_biomass_GJfuel_per_ton_dm*CropTypeArray(crops).calorificValue_weighted_MJperKgYieldModelOutput_lhv/lhv_Hansen;
        CropTypeArray(crops).ethanol_from_biomass_GJfuel_per_ton_dm = EnergyConversion_Hansen.ethanol_from_woody_biomass_GJfuel_per_ton_dm*CropTypeArray(crops).calorificValue_weighted_MJperKgYieldModelOutput_lhv/lhv_Hansen;
    end


    % Final energy, no ccs
    FEM.fe_bioelectricity(binary_this_crop) = CropTypeArray(crops).electricity_from_biomass_GJelec_per_ton_dm*FEM.dm_yield(binary_this_crop).*FEM.land_availability(binary_this_crop);
    FEM.fe_FT_diesel(binary_this_crop) = CropTypeArray(crops).ft_diesel_from_biomass_GJfuel_per_ton_dm*FEM.dm_yield(binary_this_crop).*FEM.land_availability(binary_this_crop);
    % Final energy, with ccs
    FEM.fe_bioelectricity_with_css_GJ_per_year(binary_this_crop) = FEM.fe_bioelectricity(binary_this_crop)-(EnergyConversion_Hansen.conversion_efficiency_penalty_ccs_electricity_GJelec_per_ton_dm.*FEM.dm_yield(binary_this_crop).*FEM.land_availability(binary_this_crop));
    FEM.fe_FT_diesel_with_ccs_GJ_per_year(binary_this_crop) = FEM.fe_FT_diesel(binary_this_crop)-(EnergyConversion_Hansen.conversion_efficiency_penalty_ccs_ft_diesel_GJfuel_per_ton_dm.*FEM.dm_yield(binary_this_crop).*FEM.land_availability(binary_this_crop));
end

FEM.dm = FEM.dm_yield.*FEM.land_availability;

%Fertilizer emissions
FEM.fertilizer_em_Mg_CO2eq = 10^-3*FEM.dm*EnergyConversion_Hansen.fertiliser_emissions_grassy_biomass_kgCO2eq_per_ton_dm;
FEM.emission_factor_fertilizer_electricity = 10^3*FEM.fertilizer_em_Mg_CO2eq./FEM.fe_bioelectricity;
FEM.emission_factor_fertilizer_FT = 10^3*FEM.fertilizer_em_Mg_CO2eq./FEM.fe_FT_diesel;
FEM.emission_factor_fertilizer_electricity_CCS = 10^3*FEM.fertilizer_em_Mg_CO2eq./FEM.fe_bioelectricity_with_css_GJ_per_year;
FEM.emission_factor_fertilizer_FT_CCS = 10^3*FEM.fertilizer_em_Mg_CO2eq./FEM.fe_FT_diesel_with_ccs_GJ_per_year;

% Supply chain emissions
FEM.sc_em_electricity_Mg_CO2eq = 10^-3*FEM.fe_bioelectricity*EnergyConversion_Hansen.supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ;
FEM.sc_em_FT_diesel_Mg_CO2eq = 10^-3*FEM.fe_FT_diesel*EnergyConversion_Hansen.supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ;

FEM.sc_em_ccs_electricity_Mg_CO2eq = 10^-3*FEM.fe_bioelectricity_with_css_GJ_per_year*EnergyConversion_Hansen.ccs_additional_sc_emissions_electricty_kgCO2eq_per_GJelec;
FEM.sc_em_ccs_FT_diesel_Mg_CO2eq = 10^-3*FEM.fe_FT_diesel_with_ccs_GJ_per_year*EnergyConversion_Hansen.ccs_additional_sc_emissions_fuel_kgCO2eq_per_GJfuel;

% CCS
FEM.co2_captured_electricty_Mg_CO2 = carbon_to_co2_factor*FEM.carbon_yield.*FEM.land_availability*EnergyConversion_Hansen.ccs_efficiency_electricity;
FEM.co2_captured_FT_diesel_Mg_CO2 = carbon_to_co2_factor*FEM.carbon_yield.*FEM.land_availability*EnergyConversion_Hansen.ccs_efficiency_ft_diesel;

Mg_co2_captured_per_GJelec_ccs = FEM.co2_captured_electricty_Mg_CO2./FEM.fe_bioelectricity_with_css_GJ_per_year;
Mg_co2_captured_per_GJft_ccs = FEM.co2_captured_FT_diesel_Mg_CO2./FEM.fe_FT_diesel_with_ccs_GJ_per_year;


% Emissions
FEM.gwp100_electricity_excluding_abc = FEM.sc_em_electricity_Mg_CO2eq + FEM.fertilizer_em_Mg_CO2eq;
FEM.gwp100_FT_diesel_excluding_abc = FEM.sc_em_FT_diesel_Mg_CO2eq + FEM.fertilizer_em_Mg_CO2eq;

FEM.gwp100_electricity_with_ccs_excluding_abc = FEM.sc_em_ccs_electricity_Mg_CO2eq + FEM.fertilizer_em_Mg_CO2eq + FEM.co2_captured_electricty_Mg_CO2;
FEM.gwp100_FT_diesel_with_ccs_excluding_abc = FEM.sc_em_ccs_FT_diesel_Mg_CO2eq + FEM.fertilizer_em_Mg_CO2eq + FEM.co2_captured_FT_diesel_Mg_CO2;

% Emission intensities
FEM.gwp100_intensity_electricity_excluding_abc = FEM.gwp100_electricity_excluding_abc./FEM.fe_bioelectricity;
FEM.gwp100_intensity_FT_diesel_excluding_abc = FEM.gwp100_FT_diesel_excluding_abc./FEM.fe_FT_diesel;


%% Calculate land use change emissions
FEM.land_after_abandonment_year = ac_after_abandonment_year;
FEM.current_aboveground_carbon_after_abandonment_year =  ac_abc_after_abandonment_year;
FEM.current_aboveground_carbon_per_hectare_after_abandonment_year = ac_abc_intensity_tonC_per_ha_after_abandonment_year;
FEM.time_ac = time_ac_carbon;
FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year = carbon_to_co2_factor*FEM.current_aboveground_carbon_after_abandonment_year/FEM.time_period;
FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year = ac_abc_regrowth_rate_grid_mean_2022*carbon_to_co2_factor; %% NB: THIS IS PER HECTARE. FIX VARIABLE NAME LATER.

% GJcarrier/ha
fe_yield_elec = FEM.fe_bioelectricity./FEM.land_availability;
fe_yield_elec_ccs = FEM.fe_bioelectricity_with_css_GJ_per_year./FEM.land_availability;
fe_yield_ft = FEM.fe_FT_diesel./FEM.land_availability;
fe_yield_ft_ccs = FEM.fe_FT_diesel_with_ccs_GJ_per_year./FEM.land_availability;

% Preallocation
mSize = size(FEM.land_after_abandonment_year);
FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec = zeros(mSize(1), mSize(2), mSize(3));
FEM.land_clearing_emission_intensity_kgCO2eq_per_GJft = zeros(mSize(1), mSize(2), mSize(3));
FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec_w_CCS = zeros(mSize(1), mSize(2), mSize(3));
FEM.land_clearing_emission_intensity_kgCO2eq_per_GJft_w_CCS = zeros(mSize(1), mSize(2), mSize(3));

fe_elec_after_aban_year = zeros(mSize(1), mSize(2), mSize(3));
fe_elec_ccs_after_aban_year = zeros(mSize(1), mSize(2), mSize(3));
fe_ft_after_aban_year = zeros(mSize(1), mSize(2), mSize(3));
fe_ft_ccs_after_aban_year = zeros(mSize(1), mSize(2), mSize(3));

share_of_standing_biomass_to_supply_chain = 0;

fprintf('Calculating land clearing emission intensities: ')
% Calc land clearing emission intensities (kgCO2eq per GJ final energy).
for t = 1:length(FEM.time_ac)
    fprintf([num2str(FEM.time_ac(t)) ', '])

    fe_elec_after_aban_year(:,:,t) = fe_yield_elec(:,:).*FEM.land_after_abandonment_year(:,:,t);
    fe_elec_ccs_after_aban_year(:,:,t) = fe_yield_elec_ccs(:,:).*FEM.land_after_abandonment_year(:,:,t);
    fe_ft_after_aban_year(:,:,t) = fe_yield_ft(:,:).*FEM.land_after_abandonment_year(:,:,t);
    fe_ft_ccs_after_aban_year(:,:,t) = fe_yield_ft_ccs(:,:).*FEM.land_after_abandonment_year(:,:,t);


    FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec(:,:,t) = 10^3*FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year(:,:,t)./fe_elec_after_aban_year(:,:,t);
    FEM.land_clearing_emission_intensity_kgCO2eq_per_GJft(:,:,t) = 10^3*FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year(:,:,t)./fe_ft_after_aban_year(:,:,t);

    FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec_w_CCS(:,:,t) = (10^3*FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year(:,:,t)-10^3*share_of_standing_biomass_to_supply_chain*EnergyConversion_Hansen.ccs_efficiency_electricity*FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year(:,:,t))./fe_elec_ccs_after_aban_year(:,:,t);
    FEM.land_clearing_emission_intensity_kgCO2eq_per_GJft_w_CCS(:,:,t)= (10^3*FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year(:,:,t)-10^3*share_of_standing_biomass_to_supply_chain*EnergyConversion_Hansen.ccs_efficiency_ft_diesel*FEM.initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year(:,:,t))./fe_ft_ccs_after_aban_year(:,:,t);
end
fprintf('\n')

FEM.fe_elec_after_aban_year = fe_elec_after_aban_year;
FEM.fe_elec_ccs_after_aban_year = fe_elec_ccs_after_aban_year;
FEM.fe_ft_after_aban_year = fe_ft_after_aban_year;
FEM.fe_ft_ccs_after_aban_year = fe_ft_ccs_after_aban_year;


% Preallocation
FEM.emission_factor_electricity = zeros(mSize(1), mSize(2), mSize(3));
FEM.emission_factor_FT = zeros(mSize(1), mSize(2), mSize(3));
FEM.emission_factor_electricity_CCS = zeros(mSize(1), mSize(2), mSize(3));
FEM.emission_factor_FT_CCS = zeros(mSize(1), mSize(2), mSize(3));


fprintf('Calculating net emission factors: ')
for t = 1:length(FEM.time_ac)

    % em_foregone_nr_kgCO2_per_GJelec= zeros(mSize(1), mSize(2), mSize(3));
    % em_foregone_nr_kgCO2_per_GJft= zeros(mSize(1), mSize(2), mSize(3));
    % em_foregone_nr_kgCO2_per_GJelec_CCS= zeros(mSize(1), mSize(2), mSize(3));
    % em_foregone_nr_kgCO2_per_GJft_CCS= zeros(mSize(1), mSize(2), mSize(3));
    %
    %
    % em_foregone_nr_kgCO2_per_GJelec(:,:,t) = FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./FEM.fe_bioelectricity;
    % em_foregone_nr_kgCO2_per_GJft(:,:,t) = FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./FEM.fe_FT_diesel;
    %
    % em_foregone_nr_kgCO2_per_GJelec_CCS(:,:,t) = FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./FEM.fe_bioelectricity_with_css_GJ_per_year;
    % em_foregone_nr_kgCO2_per_GJft_CCS(:,:,t) =  FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./FEM.fe_FT_diesel_with_ccs_GJ_per_year;


    fprintf([num2str(FEM.time_ac(t)) ', '])
    FEM.emission_factor_electricity(:,:,t) = (10^3*FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./fe_yield_elec(:,:))+FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec(:,:,t) + FEM.emission_factor_fertilizer_electricity(:,:) + (EnergyConversion_Hansen.supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ);
    FEM.emission_factor_FT(:,:,t) = (10^3*FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./fe_yield_ft(:,:))+FEM.land_clearing_emission_intensity_kgCO2eq_per_GJft(:,:,t) + FEM.emission_factor_fertilizer_FT(:,:) + (EnergyConversion_Hansen.supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ);

    FEM.emission_factor_electricity_CCS(:,:,t) = (10^3*FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./fe_yield_elec_ccs(:,:))+ ...
        FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec_w_CCS(:,:,t) + FEM.emission_factor_fertilizer_electricity_CCS(:,:) + ...
        (EnergyConversion_Hansen.supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ)-(10^3*Mg_co2_captured_per_GJelec_ccs(:,:))...
        +EnergyConversion_Hansen.ccs_additional_sc_emissions_electricty_kgCO2eq_per_GJelec;

    FEM.emission_factor_FT_CCS(:,:,t) = (10^3*FEM.mean_natural_regrowth_seq_rate_Mg_CO2_per_year(:,:,t)./fe_yield_ft_ccs(:,:))+FEM.land_clearing_emission_intensity_kgCO2eq_per_GJft_w_CCS(:,:,t) ...
        + FEM.emission_factor_fertilizer_FT_CCS(:,:) + (EnergyConversion_Hansen.supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ)-(10^3*Mg_co2_captured_per_GJft_ccs(:,:)) ...
        + EnergyConversion_Hansen.ccs_additional_sc_emissions_fuel_kgCO2eq_per_GJfuel;
end

FEM.emission_factor_electricity_excluding_luc = FEM.land_clearing_emission_intensity_kgCO2eq_per_GJelec(:,:,t) + FEM.emission_factor_fertilizer_electricity(:,:) + (EnergyConversion_Hansen.supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ);
FEM.emission_factor_FT_excluding_luc = FEM.emission_factor_fertilizer_FT(:,:) + (EnergyConversion_Hansen.supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ);

FEM.emission_factor_electricity_CCS_excluding_luc = FEM.emission_factor_fertilizer_electricity_CCS(:,:) + ...
    (EnergyConversion_Hansen.supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ)-(10^3*Mg_co2_captured_per_GJelec_ccs(:,:))...
    +EnergyConversion_Hansen.ccs_additional_sc_emissions_electricty_kgCO2eq_per_GJelec;

FEM.emission_factor_FT_CCS_excluding_luc = FEM.emission_factor_fertilizer_FT_CCS(:,:) + (EnergyConversion_Hansen.supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ)-(10^3*Mg_co2_captured_per_GJft_ccs(:,:)) ...
    + EnergyConversion_Hansen.ccs_additional_sc_emissions_fuel_kgCO2eq_per_GJfuel;


fprintf('\n')

%% CALCING carbon fluxes and upscaling to 1 deg

FEM = FEM.calc_carbon_fluxes_30y();

%% Design future scenarios
FEM = FEM.add_FutureAC;
FEM.FutureArray(1)

%% Export Future array to xls
export_future_to_xls = 1;
if export_future_to_xls == 1
    FEM.export_future_to_xls()
end


%% Get NOR/TRØND masks
land_NOR = FEM.land_availability_bioenergy_productive_2d;
land_NOR(~Norway_mask_5arcmin) = 0;


ncid = netcdf.open('municipal_ids.nc');
fylke_ids = netcdf.getVar(ncid,3);
id_Trondelag = 50;
netcdf.close(ncid);

mask_Trondelag = fylke_ids == id_Trondelag;
land_Trondelag = FEM.land_availability_bioenergy_productive_2d;
land_Trondelag(~mask_Trondelag) = 0;

% Fractions of bioenergy productive land availability
n_lon_1deg = 360;
n_lat_1deg = 180;
[land_1deg, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( FEM.land_availability_bioenergy_productive_2d,FEM.lon_5arcmin,FEM.lat_5arcmin, n_lon_1deg, n_lat_1deg );
[land_NOR_1deg, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( land_NOR,FEM.lon_5arcmin,FEM.lat_5arcmin, n_lon_1deg, n_lat_1deg );
[land_Trondelag_1deg, ~,  ~, ~, ~]  = aggregateMatrix2givenDimensions( land_Trondelag,FEM.lon_5arcmin,FEM.lat_5arcmin, n_lon_1deg, n_lat_1deg );

fractions_NOR = land_NOR_1deg./land_1deg;
fractions_NOR(land_1deg == 0) = 0;
fractions_Trondelag = land_Trondelag_1deg./land_1deg;
fractions_Trondelag(land_1deg == 0) = 0;

%% Plot future results using masks
FEM.plot_using_mask(fractions_NOR, 'Norway')
FEM.plot_using_mask(fractions_Trondelag, 'Trøndelag')




%% EXPORT
if export == 1

    filename = 'mapTerrain_output_5arcmin.nc';

    if exist(filename, 'file')
        delete(filename);
    end

    lat = CountryMask.lat_5arcmin;
    lon = CountryMask.lon_5arcmin;

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

    nccreate(filename, 'cropland_fractions', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'cropland_fractions', 'standard_name', 'cropland_fractions');
    ncwriteatt(filename, 'cropland_fractions', 'long_name', 'cropland_fractions');
    ncwriteatt(filename, 'cropland_fractions', 'units', 'fractions');
    ncwriteatt(filename, 'cropland_fractions', 'missing_value', '-999');

    nccreate(filename, 'abandoned_cropland_fractions', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'standard_name', 'abandoned_cropland_fractions');
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'long_name', 'abandoned_cropland_fractions');
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'units', 'fractions');
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'missing_value', '-999');

    nccreate(filename, 'ac_per_cropland', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'ac_per_cropland', 'standard_name', 'ac_per_cropland');
    ncwriteatt(filename, 'ac_per_cropland', 'long_name', 'ac_per_cropland');
    ncwriteatt(filename, 'ac_per_cropland', 'units', 'fractions');
    ncwriteatt(filename, 'ac_per_cropland', 'missing_value', '-999');

    nccreate(filename, 'ac_per_cropland_country_averages', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'ac_per_cropland_country_averages', 'standard_name', 'ac_per_cropland_country_averages');
    ncwriteatt(filename, 'ac_per_cropland_country_averages', 'long_name', 'ac_per_cropland_country_averages');
    ncwriteatt(filename, 'ac_per_cropland_country_averages', 'units', 'fractions');
    ncwriteatt(filename, 'ac_per_cropland_country_averages', 'missing_value', '-999');

    nccreate(filename, 'ac_per_capita_country_averages', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'ac_per_capita_country_averages', 'standard_name', 'ac_per_capita_country_averages');
    ncwriteatt(filename, 'ac_per_capita_country_averages', 'long_name', 'ac_per_capita_country_averages');
    ncwriteatt(filename, 'ac_per_capita_country_averages', 'units', 'fractions');
    ncwriteatt(filename, 'ac_per_capita_country_averages', 'missing_value', '-999');


    ncwrite(filename, 'lat', lat);
    ncwrite(filename, 'lon', lon);
    ncwrite(filename, 'cropland_fractions', cropland_fractions);
    ncwrite(filename, 'abandoned_cropland_fractions', ac_fractions);
    ncwrite(filename, 'ac_per_cropland', ac_per_cropland_gridded);
    ncwrite(filename, 'ac_per_cropland_country_averages', ac_per_cropland2020_country_averages)
    ncwrite(filename, 'ac_per_capita_country_averages', ac_per_capita_country_averages);




    %% 0.05deg
    %     filename = 'mapTerrain_output_005_deg.nc';
    %     if exist(filename, 'file')
    %         delete(filename);
    %     end
    %
    %     lat = CountryMask.lat_005deg;
    %     lon = CountryMask.lon_005deg;
    %
    %     nccreate(filename,'lat','Dimensions',{'lat' length(lat)});
    %     ncwriteatt(filename, 'lat', 'standard_name', 'latitude');
    %     ncwriteatt(filename, 'lat', 'long_name', 'latitude');
    %     ncwriteatt(filename, 'lat', 'units', 'degrees_north');
    %     ncwriteatt(filename, 'lat', '_CoordinateAxisType', 'Lat');
    %     ncwriteatt(filename, 'lat', 'axis', 'Y');
    %
    %     nccreate(filename,'lon','Dimensions',{'lon' length(lon)});
    %     ncwriteatt(filename, 'lon', 'standard_name', 'longitude');
    %     ncwriteatt(filename, 'lon', 'long_name', 'longitude');
    %     ncwriteatt(filename, 'lon', 'units', 'degrees_east');
    %     ncwriteatt(filename, 'lon', '_CoordinateAxisType', 'Lon');
    %     ncwriteatt(filename, 'lon', 'axis', 'X');
    %
    %     nccreate(filename, 'bioenergy_lu_gcam_2050', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050', 'standard_name', 'bioenergy_lu_gcam_2050');
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050', 'long_name', 'bioenergy_lu_gcam_2050');
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050', 'units', 'fractions');
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050', 'missing_value', '-999');
    %
    %     nccreate(filename, 'bioenergy_lu_gcam_2050_per_2020_cropland', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050_per_2020_cropland', 'standard_name', 'bioenergy_lu_gcam_2050_per_2020_cropland');
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050_per_2020_cropland', 'long_name', 'bioenergy_lu_gcam_2050_per_2020_cropland');
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050_per_2020_cropland', 'units', 'fractions');
    %     ncwriteatt(filename, 'bioenergy_lu_gcam_2050_per_2020_cropland', 'missing_value', '-999');
    %
    %

end

FEM.export_Norway


fprintf('Finished script in: ');
toc
fprintf('\n');

diary off
