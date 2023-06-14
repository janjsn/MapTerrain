classdef EnergyConversion
    %ENERGYCONVERSION Summary of this class goes here
    %   Detailed explanation goes here

    properties

        lhv_woody
        lhv_grassy

        electricity_from_woody_biomass_GJelec_per_ton_dm
        electricity_from_grassy_biomass_GJelec_per_ton_dm
        ft_diesel_from_woody_biomass_GJfuel_per_ton_dm
        ft_diesel_from_grassy_biomass_GJfuel_per_ton_dm
        ethanol_from_woody_biomass_GJfuel_per_ton_dm
        ethanol_from_grassy_biomass_GJfuel_per_ton_dm

        conversion_efficiency_penalty_ccs_electricity_GJelec_per_ton_dm
        conversion_efficiency_penalty_ccs_ft_diesel_GJfuel_per_ton_dm
        conversion_efficiency_penalty_ccs_ethanol_GJfuel_per_ton_dm

        fertiliser_emissions_woody_biomass_kgCO2eq_per_ton_dm
        fertiliser_emissions_grassy_biomass_kgCO2eq_per_ton_dm

        supply_chain_emissions_electricity_woody_biomass_kgCO2eq_per_GJ
        supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ
        supply_chain_emissions_ft_diesel_woody_biomass_kgCO2eq_per_GJ
        supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ
        supply_chain_emissions_ft_ethanol_woody_biomass_kgCO2eq_per_GJ
        supply_chain_emissions_ft_ethanol_grass_biomass_kgCO2eq_per_GJ

        ccs_additional_sc_emissions_electricty_kgCO2eq_per_GJelec
        ccs_additional_sc_emissions_fuel_kgCO2eq_per_GJfuel

        ccs_efficiency_electricity
        ccs_efficiency_ft_diesel
        ccs_efficiency_ethanol
    end

    methods
        function obj = EnergyConversion()

            % THESE NUMBERS ARE BASED ON THE LITERATURE REVIEW FROM Hansen
            % et al. (2020) in Nature Climate Change.

            obj.lhv_grassy = 18.4;
            obj.lhv_woody = 18.6;

            obj.electricity_from_woody_biomass_GJelec_per_ton_dm = 5.8;
            obj.electricity_from_grassy_biomass_GJelec_per_ton_dm = 5.7;
            obj.ft_diesel_from_woody_biomass_GJfuel_per_ton_dm = 8.1;
            obj.ft_diesel_from_grassy_biomass_GJfuel_per_ton_dm = 8.0;

            obj.ethanol_from_woody_biomass_GJfuel_per_ton_dm = 7.1;
            obj.ethanol_from_grassy_biomass_GJfuel_per_ton_dm = 7.0;

            obj.conversion_efficiency_penalty_ccs_electricity_GJelec_per_ton_dm = 1.8;
            obj.conversion_efficiency_penalty_ccs_ft_diesel_GJfuel_per_ton_dm = 0;
            obj.conversion_efficiency_penalty_ccs_ethanol_GJfuel_per_ton_dm = 1.9;

            obj.fertiliser_emissions_woody_biomass_kgCO2eq_per_ton_dm = 55;
            obj.fertiliser_emissions_grassy_biomass_kgCO2eq_per_ton_dm = 54;

            obj.supply_chain_emissions_electricity_woody_biomass_kgCO2eq_per_GJ = 13;
            obj.supply_chain_emissions_electricity_grass_biomass_kgCO2eq_per_GJ = 16;
            obj.supply_chain_emissions_ft_diesel_woody_biomass_kgCO2eq_per_GJ = 19;
            obj.supply_chain_emissions_ft_diesel_grass_biomass_kgCO2eq_per_GJ = 18;
            obj.supply_chain_emissions_ft_ethanol_woody_biomass_kgCO2eq_per_GJ = 14;
            obj.supply_chain_emissions_ft_ethanol_grass_biomass_kgCO2eq_per_GJ = 20;

            obj.ccs_additional_sc_emissions_electricty_kgCO2eq_per_GJelec = 11;
            obj.ccs_additional_sc_emissions_fuel_kgCO2eq_per_GJfuel = 3;

            obj.ccs_efficiency_electricity = 0.90;
            obj.ccs_efficiency_ft_diesel = 0.52;
            obj.ccs_efficiency_ethanol = 0.12;


        end

        function outputArg = method1(obj,inputArg)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            outputArg = obj.Property1 + inputArg;
        end
    end
end

