classdef FutureAC
    %FUTUREAC Summary of this class goes here
    %   Detailed explanation goes here

    properties

        scenario_description
        lon
        lat
        abandonment_half_life
        share_of_recultivation_to_food_production
        share_of_recultivation_to_beccs_electricity
        share_of_recultivation_to_beccs_ft

        avoided_emissions_factor_kgCO2eq_per_GJ_electricity
        avoided_emissions_factor_kgCO2eq_per_GJ_ft

        land_availability_2d_t0
        Mgco2_stock_on_ac_2d_t0
        natural_regrowth_rate_Mgco2_per_ha_year_2d
        fe_yield_bioelectricity_beccs_2d
        fe_yield_ft_beccs_2d

        time
        time_since_t0

        share_of_initial_land_remaining
        share_of_initial_land_converted_to_food_production
        share_of_initial_land_converted_to_beccs
        share_of_initial_land_converted_to_beccs_electricity
        share_of_initial_land_converted_to_beccs_ft

        land_regrowing
        land_recultivated_food_or_bioenergy
        land_converted_to_food_production
        land_converted_to_bioenergy
        land_converted_to_beccs_electricity
        land_converted_to_beccs_ft

        totals_vec_land_regrowing
        totals_vec_land_recultivated_food_or_bioenergy
        totals_vec_land_converted_to_food_production
        totals_vec_land_converted_to_bioenergy
        totals_vec_land_converted_to_beccs_electricity
        totals_vec_land_converted_to_beccs_ft
        
        tco2eq_aboveground_carbon_stocks
        
        % FLUXES in a given year
        tco2eq_flux_natural_regrowth_added
        tco2eq_flux_natural_regrowth_removed
        tco2eq_flux_natural_regrowth_net

        tco2eq_flux_beccs_electricity
        tco2eq_flux_beccs_ft
        tco2eq_flux_beccs_tot

        totals_vec_tco2eq_flux_natural_regrowth_added
        totals_vec_tco2eq_flux_natural_regrowth_removed
        totals_vec_tco2eq_flux_natural_regrowth_net

        totals_vec_tco2eq_flux_beccs_electricity
        totals_vec_tco2eq_flux_beccs_ft
        totals_vec_tco2eq_flux_beccs_tot

        totals_vec_tco2eq_flux_net_nr_beccs

        %Final energy produced per year
        fe_bioelectricity_ccs
        fe_ft_ccs
        fe_ccs_tot

        totals_vec_fe_bioelectricity_ccs
        totals_vec_fe_ft_ccs
        totals_vec_fe_ccs_tot

        % Cumulative final energy
        cumulative_fe_bioelectricity_ccs
        cumulative_fe_ft_ccs
        cumulative_fe_bioelectricity_ccs_tot
        cumulative_fe_ft_ccs_tot
        
        % CUMULATIVE OVER TIME PERIOD
        cumulative_tco2_no_substitution_natural_regrowth
        cumulative_tco2_no_substitution_beccs

        cumulative_tco2_no_substitution
        cumulative_tco2_with_substitution

        cumulative_tco2_no_substitution_natural_regrowth_tot
        cumulative_tco2_no_substitution_beccs_tot
        
        cumulative_tco2_no_substitution_tot
        cumulative_tco2_with_substitution_diesel_tot
        
        cumulative_tco2_with_substitution_elec_fos_ccs_tot
        cumulative_tco2_with_substitution_elec_natural_gas_tot
        cumulative_tco2_with_substitution_elec_coal_tot

        %delta_cumulative_tco2_no_substitution
        %delta_cumulative_tco2_with_substitution
        
        cumulative_tco2_with_substitution_diesel
        
        cumulative_tco2_with_substitution_elec_fos_ccs
        cumulative_tco2_with_substitution_elec_natural_gas
        cumulative_tco2_with_substitution_elec_coal
    end

    methods

    end
end

