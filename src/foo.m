% test = FEM.calc_carbon_fluxes_30y
%
%
% test.carbon_flux_30y_avg_FT_CCS_Mg_CO2eq_per_year_tot+sum(sum(sum(seq_co2_nr(test.fe_ft_ccs_after_aban_year > 0))))



time = [2015:5:2100];

global_ssp1_proj = zeros(1,length(time));
global_ssp2_proj = zeros(1,length(time));
global_ssp4_proj = zeros(1,length(time));
global_ssp5_proj = zeros(1,length(time));

unique_continents = unique({CountryMask.CountryArray.GPW_continent_string});
n_continents = length(unique_continents);

continental_ssp1_proj = zeros(n_continents,length(time));
continental_ssp2_proj = zeros(n_continents,length(time));
continental_ssp4_proj = zeros(n_continents,length(time));
continental_ssp5_proj = zeros(n_continents,length(time));

continental_pe = zeros(1,n_continents);

do_gcam_calcs = 0;
if do_gcam_calcs == 1
    for i = 1:length(CountryMask.CountryArray)
        global_ssp1_proj = global_ssp1_proj + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
        global_ssp2_proj = global_ssp2_proj + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp2_rcp26_vec;
        global_ssp4_proj = global_ssp4_proj + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp4_rcp26_vec;
        global_ssp5_proj = global_ssp5_proj + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp5_rcp26_vec;

        for j = 1:length(unique_continents)
            if strcmp(unique_continents{j}, CountryMask.CountryArray(i).GPW_continent_string)
                continental_ssp1_proj(j,:) = continental_ssp1_proj(j,:) + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
                continental_ssp2_proj(j,:) = continental_ssp1_proj(j,:) + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
                continental_ssp4_proj(j,:) = continental_ssp1_proj(j,:) + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
                continental_ssp5_proj(j,:) = continental_ssp1_proj(j,:) + CountryMask.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;

                continental_pe(j) = continental_pe(j) + FEM.pe(CountryMask.countryMask_5arcmin == CountryMask.CountryArray(i).GPW_country_ISO_numeric);

            end


        end
    end
end

continental_pe = zeros(1,n_continents);
for i = 1:length(CountryMask.CountryArray)
    for j = 1:length(unique_continents)
        if strcmp(unique_continents{j}, CountryMask.CountryArray(i).GPW_continent_string)
            continental_pe(j) = continental_pe(j) + sum(sum(FEM.pe(CountryMask.countryMask_5arcmin == CountryMask.CountryArray(i).GPW_country_ISO_numeric)));
        end
    end
end


productive_land_3d = FEM.land_after_abandonment_year;

for i = 1:29
    this = zeros(4320,2160);
    this(:,:) = productive_land_3d(:,:,i);
    this(FEM.pe <=0) = 0;
    this(~Norway_mask_5arcmin) = 0;
    productive_land_3d(:,:,i) = this;

end

sum(sum(sum(productive_land_3d(FEM.emission_factor_electricity_CCS <0 ))))

sum(sum(sum(productive_land_3d(FEM.emission_factor_FT_CCS <0) )))


