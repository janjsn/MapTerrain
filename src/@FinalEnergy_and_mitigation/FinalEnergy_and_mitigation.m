classdef FinalEnergy_and_mitigation
    %FINALENERGY_AND_MITIGATION Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        time_period 
        lat_5arcmin
        lon_5arcmin
        land_availability
        crop_allocation
        binary_is_grassy
        binary_is_woody
        dm_yield
        carbon_yield
        pe_yield
        dm
        pe
        unit_fe = 'GJ year-1';
        fe_bioelectricity
        fe_FT_diesel
        
        aboveground_carbon_on_land % Mg C
        time_ac
        land_after_abandonment_year % ha
        current_aboveground_carbon_after_abandonment_year % Mg C
        current_aboveground_carbon_per_hectare_after_abandonment_year % Mg C ha-1
        mean_natural_regrowth_seq_rate_Mg_CO2_per_year %% NB: THIS IS PER HECTARE. FIX VARIABLE NAME LATER.
        initial_em_pulse_Mg_co2_abc_clearance_avg_30y_after_aban_year
        
        fertilizer_em_Mg_CO2eq

        sc_em_electricity_Mg_CO2eq
        sc_em_FT_diesel_Mg_CO2eq

        gwp100_electricity_excluding_abc %kg CO2eq
        gwp100_FT_diesel_excluding_abc %kg CO2eq

        gwp100_intensity_electricity_excluding_abc %kg CO2eq
        gwp100_intensity_FT_diesel_excluding_abc %kg CO2eq
        
        fe_bioelectricity_with_css_GJ_per_year
        fe_FT_diesel_with_ccs_GJ_per_year

        sc_em_ccs_electricity_Mg_CO2eq %kg CO2eq
        sc_em_ccs_FT_diesel_Mg_CO2eq %kg CO2eq

        co2_captured_electricty_Mg_CO2 % Mg CO2
        co2_captured_FT_diesel_Mg_CO2 % Mg CO2

        gwp100_electricity_with_ccs_excluding_abc %kg CO2eq
        gwp100_FT_diesel_with_ccs_excluding_abc %kg CO2eq

        land_clearing_emission_intensity_kgCO2eq_per_GJelec
        land_clearing_emission_intensity_kgCO2eq_per_GJft
        
        land_clearing_emission_intensity_kgCO2eq_per_GJelec_w_CCS
        land_clearing_emission_intensity_kgCO2eq_per_GJft_w_CCS

        lost_natural_regrowth_em_intensity_kgCO2eq_per_GJelec
        lost_natural_regrowth_em_intensity_kgCO2eq_per_GJft
        lost_natural_regrowth_em_intensity_kgCO2eq_per_GJelec_w_CCS
        lost_natural_regrowth_em_intensity_kgCO2eq_per_GJft_w_CCS
        
        emission_factor_fertilizer_electricity
        emission_factor_fertilizer_FT
        emission_factor_fertilizer_electricity_CCS
        emission_factor_fertilizer_FT_CCS

        emission_factor_electricity_excluding_luc
        emission_factor_FT_excluding_luc
        emission_factor_electricity_CCS_excluding_luc
        emission_factor_FT_CCS_excluding_luc

        emission_factor_electricity
        emission_factor_FT
        emission_factor_electricity_CCS
        emission_factor_FT_CCS
        

        fe_elec_after_aban_year
        fe_elec_ccs_after_aban_year
        fe_ft_after_aban_year
        fe_ft_ccs_after_aban_year

        fe_elec_tot
        fe_elec_ccs_tot
        fe_ft_tot
        fe_ft_ccs_tot
        
        land_availability_bioenergy_productive_2d
        current_aboveground_carbon_MgCO2eq_stock_2d
        fe_ft_ccs_2d
        fe_bioelectricity_ccs_2d

        carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year
        carbon_flux_30y_avg_FT_Mg_CO2eq_per_year
        carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year
        carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year
        carbon_flux_30y_avg_continued_natural_regrowth_Mg_CO2eq_pr_yr
        note_cf_coninuted_natural_regrowth = 'Only areas with bioenergy productivity';

        delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year
        delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year
        delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year
        delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year
        
        land_availability_bioenergy_productive_tot 
        current_aboveground_carbon_MgCO2eq_stock_2d_tot
        
        carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_tot
        carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_tot
        carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_tot
        carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_tot
        carbon_flux_30y_avg_continued_natural_regr_Mg_CO2eq_pr_yr_tot

        delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_tot
        delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_tot
        delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_year_tot
        delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_tot

        carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_2d
        carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_2d
        carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_2d
        carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_2d
        carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_yr_tot_2d
        carbon_flux_30y_avg_continued_nat_regr_Mg_CO2eq_pr_ha_yr_tot_2d

        delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_2d
        delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_2d
        delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_yr_2d
        delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_2d

        
        lat_1deg
        lon_1deg
        land_availability_bioenergy_productive_1deg_2d
        pe_1deg_2d
        pe_yield_1deg_2d 
        current_aboveground_carbon_MgCO2eq_stock_1deg_2d
        fe_bioelectricity_ccs_1deg_2d
        fe_ft_ccs_1deg_2d
        fe_yield_bioelectricity_ccs_1deg_2d
        fe_yield_ft_ccs_1deg_2d

        carbon_flux_30y_avg_electricity_Mg_CO2eq_per_year_1_deg_2d
        carbon_flux_30y_avg_FT_Mg_CO2eq_per_year_1_deg_2d
        carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_1_deg_2d
        carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_1_deg_2d

        delta_carbon_flux_be_nr_electricity_Mg_CO2eq_per_year_1_deg_2d
        delta_carbon_flux_be_nr_FT_Mg_CO2eq_per_year_1_deg_2d
        delta_carbon_flux_be_nr_electricity_CCS_Mg_CO2eq_per_yr_1deg2d
        delta_carbon_flux_be_nr_FT_CCS_Mg_CO2eq_per_year_1_deg_2d

        natural_regrowth_carbon_flux_Mg_CO2eq_per_year_1_deg_2d
        natural_regrowth_carbon_flux_Mg_CO2eq_per_ha_year_1_deg_2d

        emission_factor_electricity_excluding_luc_1_deg_2d
        emission_factor_FT_excluding_luc_1_deg_2d
        emission_factor_electricity_CCS_excluding_luc_1_deg_2d
        emission_factor_FT_CCS_excluding_luc_1_deg_2d

        FutureArray

    end
    
    methods
        upscale_and_export_maps_1deg(obj)
        export_Norway(obj, Norway_mask)
        plot_carbon_fluxes(obj)
        plot_final_energy(obj)
        obj = calc_carbon_fluxes_30y(obj)
        plot_using_mask(obj, fraction_of_cell_is_region_1deg, region_name)
    end
end

