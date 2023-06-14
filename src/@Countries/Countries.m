classdef Countries
    %COUNTRIES Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        country_name
        GPW_country_ISO_numeric = -999;
        GPW_country_ISO_string = '-';
        GPW_continent_string = '-';
        
        population_2020

        bioenergy_land_use_GCAM_2050_per_SSP
        bioenergy_land_use_GCAM_2100_per_SSP

        bioenergy_land_use_GCAM_2050_per_cap_per_SSP
        bioenergy_land_use_GCAM_2100_per_cap_per_SSP

        bioenergy_land_use_GCAM_ssp1_rcp26_vec
        bioenergy_land_use_GCAM_ssp2_rcp26_vec
        bioenergy_land_use_GCAM_ssp4_rcp26_vec
        bioenergy_land_use_GCAM_ssp5_rcp26_vec
        time_vec = [2015:5:2100];

        abandoned_cropland_1992_to_2020_tot
        ac_per_cap

        cropland_2020_tot
        cropland_per_cap

        abandoned_cropland_as_share_of_2020_cropland

        pe_ac_cruts32_tot 
        pe_ac_noresm_rcp45_2020s_tot

        pe_ac_cruts32_per_cap
        pe_ac_noresm_rcp45_2020s_per_cap
        

        


    end
    
    methods
    end
    
end

