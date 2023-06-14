classdef Biorefinery
    %BIOREFINERY Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        description
        year
        energy_efficiency
        ethanol_efficiency
        fischer_tropsch_efficiency
        share_ft_to_diesel
        share_ft_to_gasoline
        electricity_efficiency
        fraction_biomass_C_in_fuel
        fraction_biomass_C_emitted_at_refinery
        fraction_biomass_C_sequestered_through_CCS_or_char
        blue_water_footprint_m3_per_kg_biomass
        missing_value = -999;
        references
        
    end
    
    methods
        function obj = Biorefinery(description, year, energy_efficiency, ethanol_efficiency, fischer_tropsch_efficiency, electricity_efficiency, fraction_biomass_C_in_fuel, fraction_biomass_C_emitted_at_refinery, fraction_biomass_C_sequestered_through_CCS_or_char, blue_water_footprint_m3_per_kg_biomass, share_ft_to_diesel, share_ft_to_gasoline, references)
            if nargin > 0
                obj.description = description;
                obj.year = year;
                obj.energy_efficiency = energy_efficiency;
                obj.ethanol_efficiency = ethanol_efficiency;
                obj.fischer_tropsch_efficiency = fischer_tropsch_efficiency;
                obj.electricity_efficiency = electricity_efficiency;
                obj.fraction_biomass_C_in_fuel =  fraction_biomass_C_in_fuel;
                obj.fraction_biomass_C_emitted_at_refinery = fraction_biomass_C_emitted_at_refinery;
                obj.fraction_biomass_C_sequestered_through_CCS_or_char = fraction_biomass_C_sequestered_through_CCS_or_char;
                obj.blue_water_footprint_m3_per_kg_biomass = blue_water_footprint_m3_per_kg_biomass;
                obj.share_ft_to_diesel = share_ft_to_diesel;
                obj.share_ft_to_gasoline = share_ft_to_gasoline;
                obj.references = references;
            end
        end
    end
    
end

