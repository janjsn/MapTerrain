function obj= add_FutureAC(obj)

n_scens = 1;

Future(1:n_scens) = FutureAC;
% Baseline
this_scen = 1;
Future(this_scen).scenario_description = 'R23_F100_EL0_FT0';
Future(this_scen).abandonment_half_life = 23;
Future(this_scen).share_of_recultivation_to_food_production = 1;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 0;
Future(this_scen).share_of_recultivation_to_beccs_ft = 0;

% No recultivation
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 10^100;
Future(this_scen).scenario_description = 'Rinf_F100_EL0_FT0';
Future(this_scen).share_of_recultivation_to_food_production = 0.5;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 0.5;
Future(this_scen).share_of_recultivation_to_beccs_ft = 0;

% halflife 23y, to beccs electricity
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 23;
Future(this_scen).scenario_description = 'R23_F0_EL100_FT0';
Future(this_scen).share_of_recultivation_to_food_production = 0;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 1;
Future(this_scen).share_of_recultivation_to_beccs_ft = 0;

% halflife 23y, to beccs ft
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 23;
Future(this_scen).scenario_description = 'R23_F0_EL0_FT100';
Future(this_scen).share_of_recultivation_to_food_production = 0;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 0;
Future(this_scen).share_of_recultivation_to_beccs_ft = 1;

% halflife 23y, to beccs 50-50 electricity and ft
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 23;
Future(this_scen).scenario_description = 'R23_F0_EL50_FT50';
Future(this_scen).share_of_recultivation_to_food_production = 0;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 0.5;
Future(this_scen).share_of_recultivation_to_beccs_ft = 0.5;

% halflife 11y, to beccs electricity
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 11;
Future(this_scen).scenario_description = 'R11_F0_EL100_FT0';
Future(this_scen).share_of_recultivation_to_food_production = 0;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 1;
Future(this_scen).share_of_recultivation_to_beccs_ft = 0;

% halflife 11y, to beccs ft
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 11;
Future(this_scen).scenario_description = 'R11_F0_EL0_FT100';
Future(this_scen).share_of_recultivation_to_food_production = 0;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 0;
Future(this_scen).share_of_recultivation_to_beccs_ft = 1;

% halflife 11y, to beccs 50-50 electricity and ft
this_scen = this_scen+1;
Future(this_scen).abandonment_half_life = 11;
Future(this_scen).scenario_description = 'R11_F0_EL50_FT50';
Future(this_scen).share_of_recultivation_to_food_production = 0;
Future(this_scen).share_of_recultivation_to_beccs_electricity = 0.5;
Future(this_scen).share_of_recultivation_to_beccs_ft = 0.5;

for i = 1:length(Future)
    Future(i).lat = obj.lat_1deg;
    Future(i).lon = obj.lon_1deg;
    Future(i).time = [2022:1:2053];
    Future(i).time_since_t0 = Future(i).time-Future(i).time(1);

    Future(i).land_availability_2d_t0 = obj.land_availability_bioenergy_productive_1deg_2d;
    Future(i).Mgco2_stock_on_ac_2d_t0 = obj.current_aboveground_carbon_MgCO2eq_stock_1deg_2d;
    Future(i).natural_regrowth_rate_Mgco2_per_ha_year_2d = obj.natural_regrowth_carbon_flux_Mg_CO2eq_per_ha_year_1_deg_2d;
    for t = 1:length(Future(i).time)
        Future(i).share_of_initial_land_remaining(t) = 1*0.5^(Future(i).time_since_t0(t)/Future(i).abandonment_half_life);
    end
    Future(i).share_of_initial_land_converted_to_food_production = (1-Future(i).share_of_initial_land_remaining)*Future(i).share_of_recultivation_to_food_production;
    Future(i).share_of_initial_land_converted_to_beccs = (1-Future(i).share_of_initial_land_remaining)*(Future(i).share_of_recultivation_to_beccs_electricity+Future(i).share_of_recultivation_to_beccs_ft);
    Future(i).share_of_initial_land_converted_to_beccs_electricity = (1-Future(i).share_of_initial_land_remaining)*Future(i).share_of_recultivation_to_beccs_electricity;
    Future(i).share_of_initial_land_converted_to_beccs_ft = (1-Future(i).share_of_initial_land_remaining)*Future(i).share_of_recultivation_to_beccs_ft;

    % Preallocation
    Future(i).land_regrowing = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).land_recultivated_food_or_bioenergy = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).land_converted_to_food_production = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).land_converted_to_bioenergy = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).land_converted_to_beccs_electricity = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).land_converted_to_beccs_ft = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).tco2eq_aboveground_carbon_stocks = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).tco2eq_flux_natural_regrowth_added = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).tco2eq_flux_beccs_electricity = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).tco2eq_flux_beccs_ft = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).tco2eq_flux_beccs_tot = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).cumulative_tco2_no_substitution = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).cumulative_tco2_with_substitution = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).fe_bioelectricity_ccs = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).fe_ft_ccs = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).tco2eq_flux_natural_regrowth_removed= zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).tco2eq_flux_natural_regrowth_net= zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).cumulative_tco2_no_substitution_natural_regrowth = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).cumulative_tco2_no_substitution_beccs = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));

    Future(i).cumulative_fe_bioelectricity_ccs = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));
    Future(i).cumulative_fe_ft_ccs = zeros(length(Future(i).lon), length(Future(i).lat), length(Future(i).time));


    % Calc final energy yields
    Future(i).fe_yield_bioelectricity_beccs_2d = obj.fe_yield_bioelectricity_ccs_1deg_2d;
    Future(i).fe_yield_ft_beccs_2d = obj.fe_yield_ft_ccs_1deg_2d;

    for t = 1:length(Future(i).time)
        Future(i).land_regrowing(:,:,t) = Future(i).land_availability_2d_t0.*Future(i).share_of_initial_land_remaining(t);
        Future(i).land_recultivated_food_or_bioenergy(:,:,t) = Future(i).land_availability_2d_t0-Future(i).land_regrowing(:,:,t);
        Future(i).land_converted_to_food_production(:,:,t) = Future(i).land_recultivated_food_or_bioenergy(:,:,t)*Future(i).share_of_recultivation_to_food_production;
        Future(i).land_converted_to_beccs_electricity(:,:,t) = Future(i).land_recultivated_food_or_bioenergy(:,:,t)*Future(i).share_of_recultivation_to_beccs_electricity;
        Future(i).land_converted_to_beccs_ft(:,:,t) = Future(i).land_recultivated_food_or_bioenergy(:,:,t)*Future(i).share_of_recultivation_to_beccs_ft ;
        Future(i).land_converted_to_bioenergy(:,:,t) = Future(i).land_converted_to_beccs_electricity(:,:,t)+Future(i).land_converted_to_beccs_ft(:,:,t);
        if t == 1
            Future(i).tco2eq_aboveground_carbon_stocks(:,:,t) = Future(i).Mgco2_stock_on_ac_2d_t0;
        else
            share_of_land_remaining_from_previous_time_step = Future(i).share_of_initial_land_remaining(t)/Future(i).share_of_initial_land_remaining(t-1);

            Future(i).tco2eq_aboveground_carbon_stocks(:,:,t) = Future(i).tco2eq_aboveground_carbon_stocks(:,:,t-1)*share_of_land_remaining_from_previous_time_step + ...
                (Future(i).land_regrowing(:,:,t).*Future(i).natural_regrowth_rate_Mgco2_per_ha_year_2d );

            Future(i).tco2eq_flux_natural_regrowth_removed(:,:,t) = Future(i).tco2eq_aboveground_carbon_stocks(:,:,t-1)*(1-share_of_land_remaining_from_previous_time_step);
        end

        Future(i).tco2eq_flux_natural_regrowth_added(:,:,t) = -Future(i).land_regrowing(:,:,t).*Future(i).natural_regrowth_rate_Mgco2_per_ha_year_2d ;
        Future(i).tco2eq_flux_natural_regrowth_net(:,:,t) = Future(i).tco2eq_flux_natural_regrowth_added(:,:,t)+Future(i).tco2eq_flux_natural_regrowth_removed(:,:,t);

        Future(i).fe_bioelectricity_ccs(:,:,t) = Future(i).land_converted_to_beccs_electricity(:,:,t).*Future(i).fe_yield_bioelectricity_beccs_2d ;
        Future(i).fe_ft_ccs(:,:,t) = Future(i).land_converted_to_beccs_ft(:,:,t).*Future(i).fe_yield_ft_beccs_2d ;
        Future(i).fe_ccs_tot(:,:,t) = Future(i).fe_bioelectricity_ccs(:,:,t) + Future(i).fe_ft_ccs(:,:,t);

        %Future(i).tco2eq_flux_beccs_electricity(:,:,t) = obj.carbon_flux_30y_avg_electricity_CCS_Mg_CO2eq_per_year_1_deg_2d*Future(i).share_of_initial_land_converted_to_beccs_electricity(t);
        %Future(i).tco2eq_flux_beccs_ft(:,:,t) = obj.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_1_deg_2d*Future(i).share_of_initial_land_converted_to_beccs_ft(t);
        %Future(i).tco2eq_flux_beccs_tot(:,:,t) = Future(i).tco2eq_flux_beccs_electricity(:,:,t) + Future(i).tco2eq_flux_beccs_ft(:,:,t);


        % Calc co2 flux from beccs
        Future(i).tco2eq_flux_beccs_electricity(:,:,t) = 10^-3*Future(i).land_converted_to_beccs_electricity(:,:,t).*obj.fe_yield_bioelectricity_ccs_1deg_2d(:,:) ...
            .*obj.emission_factor_electricity_CCS_excluding_luc_1_deg_2d(:,:);

        Future(i).tco2eq_flux_beccs_ft(:,:,t) = 10^-3*Future(i).land_converted_to_beccs_ft(:,:,t).*obj.fe_yield_ft_ccs_1deg_2d(:,:) ...
            .*obj.emission_factor_FT_CCS_excluding_luc_1_deg_2d(:,:);

        Future(i).tco2eq_flux_beccs_tot(:,:,t) = Future(i).tco2eq_flux_beccs_electricity(:,:,t) + Future(i).tco2eq_flux_beccs_ft(:,:,t);

        Future(i).totals_vec_tco2eq_flux_natural_regrowth_added(t) = 10^-9*sum(sum(Future(i).tco2eq_flux_natural_regrowth_added(:,:,t)));
        Future(i).totals_vec_tco2eq_flux_natural_regrowth_removed(t) = 10^-9*sum(sum(Future(i).tco2eq_flux_natural_regrowth_removed(:,:,t)));
        Future(i).totals_vec_tco2eq_flux_natural_regrowth_net(t) = 10^-9*sum(sum(Future(i).tco2eq_flux_natural_regrowth_net(:,:,t)));

        Future(i).totals_vec_tco2eq_flux_beccs_electricity(t) = 10^-9*sum(sum(Future(i).tco2eq_flux_beccs_electricity(:,:,t)));
        Future(i).totals_vec_tco2eq_flux_beccs_ft(t) = 10^-9*sum(sum(Future(i).tco2eq_flux_beccs_ft(:,:,t)));
        Future(i).totals_vec_tco2eq_flux_beccs_tot(t) = 10^-9*sum(sum(Future(i).tco2eq_flux_beccs_tot(:,:,t)));

        Future(i).totals_vec_tco2eq_flux_net_nr_beccs(t) = Future(i).totals_vec_tco2eq_flux_beccs_tot(t)+Future(i).totals_vec_tco2eq_flux_natural_regrowth_net(t);

        %Cumulative

        if t == 1
            Future(i).cumulative_tco2_no_substitution(:,:,t) = Future(i).tco2eq_flux_natural_regrowth_net(:,:,t)+Future(i).tco2eq_flux_beccs_tot(:,:,t);
            Future(i).cumulative_tco2_no_substitution_natural_regrowth(:,:,t) = Future(i).tco2eq_flux_natural_regrowth_net(:,:,t);
            Future(i).cumulative_tco2_no_substitution_beccs(:,:,t) = Future(i).tco2eq_flux_beccs_tot(:,:,t);
        elseif t > 1
            Future(i).cumulative_tco2_no_substitution(:,:,t) = Future(i).cumulative_tco2_no_substitution(:,:,t-1) + Future(i).tco2eq_flux_natural_regrowth_net(:,:,t)+Future(i).tco2eq_flux_beccs_tot(:,:,t);

            Future(i).cumulative_tco2_no_substitution_beccs(:,:,t) = Future(i).cumulative_tco2_no_substitution_beccs(:,:,t-1) + Future(i).tco2eq_flux_beccs_tot(:,:,t);
            Future(i).cumulative_tco2_no_substitution_natural_regrowth(:,:,t) = Future(i).cumulative_tco2_no_substitution_natural_regrowth(:,:,t-1) + Future(i).tco2eq_flux_natural_regrowth_net(:,:,t);
        end

        Future(i).cumulative_tco2_no_substitution_natural_regrowth_tot(t) = 10^-9*sum(sum(sum(Future(i).cumulative_tco2_no_substitution_natural_regrowth(:,:,t))));
        Future(i).cumulative_tco2_no_substitution_beccs_tot(t) = 10^-9*sum(sum(sum(Future(i).cumulative_tco2_no_substitution_beccs(:,:,t))));
        Future(i).cumulative_tco2_no_substitution_tot(t) = 10^-9*sum(sum(sum(Future(i).cumulative_tco2_no_substitution(:,:,t))));
        



        %% Calculate totals
        Future(i).totals_vec_fe_bioelectricity_ccs(t) = 10^-9*sum(sum(sum(Future(i).fe_bioelectricity_ccs(:,:,t))));
        Future(i).totals_vec_fe_ft_ccs(t) = 10^-9*sum(sum(sum(Future(i).fe_ft_ccs(:,:,t))));
        Future(i).totals_vec_fe_ccs_tot(t) = 10^-9*sum(sum(sum(Future(i).fe_ccs_tot(:,:,t))));

        Future(i).totals_vec_land_regrowing(t) = 10^-6*sum(sum(Future(i).land_regrowing(:,:,t)));
        Future(i).totals_vec_land_recultivated_food_or_bioenergy(t) = 10^-6*sum(sum(Future(i).land_recultivated_food_or_bioenergy(:,:,t)));
        Future(i).totals_vec_land_converted_to_food_production(t) = 10^-6*sum(sum(Future(i).land_converted_to_food_production(:,:,t)));
        Future(i).totals_vec_land_converted_to_bioenergy(t) = 10^-6*sum(sum(Future(i).land_converted_to_bioenergy(:,:,t)));
        Future(i).totals_vec_land_converted_to_beccs_electricity(t) = 10^-6*sum(sum(Future(i).land_converted_to_beccs_electricity(:,:,t)));
        Future(i).totals_vec_land_converted_to_beccs_ft(t) = 10^-6*sum(sum(Future(i).land_converted_to_beccs_ft(:,:,t)));

        %Future(i).totals_vec_land_converted_to_beccs_ft(t) 

        %% Cumulative final energy
        if t == 1
            Future(i).cumulative_fe_bioelectricity_ccs(:,:,t) = Future(i).fe_bioelectricity_ccs(:,:,t);
            Future(i).cumulative_fe_ft_ccs(:,:,t) = Future(i).fe_ft_ccs(:,:,t);
        elseif t > 1
            Future(i).cumulative_fe_bioelectricity_ccs(:,:,t) = Future(i).cumulative_fe_bioelectricity_ccs(:,:,t-1) + Future(i).fe_bioelectricity_ccs(:,:,t);
            Future(i).cumulative_fe_ft_ccs(:,:,t) = Future(i).cumulative_fe_ft_ccs(:,:,t-1) + Future(i).fe_ft_ccs(:,:,t);
        end

        Future(i).cumulative_fe_bioelectricity_ccs_tot(t) = sum(sum(sum(Future(i).cumulative_fe_bioelectricity_ccs(:,:,t) )));
        Future(i).cumulative_fe_ft_ccs_tot(t) = sum(sum(sum(Future(i).cumulative_fe_ft_ccs(:,:,t) )));

        %% CUMULATIVES WITH SUBSTITUTION
        % REF Scarlat et al. (2022), Applied energy
        emf_Norway_elec_prod = 7.8;
        emf_EU27_elec_prod = 86.1;
        % Emission factors, alternative techs
        
        emf_fossil_with_ccs_elec = [44 73];
        emf_natural_gas_elec = [136 146];
        emf_coal_elec = [220 259];
        
        emf_diesel_fuel = 93.9;

        Future(i).cumulative_tco2_with_substitution_diesel_tot(t) = Future(i).cumulative_tco2_no_substitution_tot(t)-10^-9*Future(i).cumulative_fe_ft_ccs_tot(t)*(10^-3*emf_diesel_fuel);
        Future(i).cumulative_tco2_with_substitution_elec_fos_ccs_tot(t) = Future(i).cumulative_tco2_no_substitution_tot(t)-10^-9*Future(i).cumulative_fe_bioelectricity_ccs_tot(t)*mean(10^-3*emf_fossil_with_ccs_elec);
        Future(i).cumulative_tco2_with_substitution_elec_natural_gas_tot(t) = Future(i).cumulative_tco2_no_substitution_tot(t)-10^-9*Future(i).cumulative_fe_bioelectricity_ccs_tot(t)*mean(10^-3*emf_natural_gas_elec);
        Future(i).cumulative_tco2_with_substitution_elec_coal_tot(t) = Future(i).cumulative_tco2_no_substitution_tot(t)-10^-9*Future(i).cumulative_fe_bioelectricity_ccs_tot(t)*mean(10^-3*emf_coal_elec);


    end
    Future(i).scenario_description
    Future(i).cumulative_fe_ft_ccs_tot
    Future(i).cumulative_fe_bioelectricity_ccs_tot
end

obj.FutureArray = Future;


end

