function plot_ssp_be_lu_per_cap(CountryMask)

for i = 1:length(CountryMask.CountryArray)
    be_lu_2050_ssp1(i) = CountryMask.CountryArray(i).bioenergy_land_use_GCAM_2050_per_cap_per_SSP(1);
    be_lu_2050_ssp2(i) = CountryMask.CountryArray(i).bioenergy_land_use_GCAM_2050_per_cap_per_SSP(2);
    be_lu_2050_ssp4(i) = CountryMask.CountryArray(i).bioenergy_land_use_GCAM_2050_per_cap_per_SSP(3);
    be_lu_2050_ssp5(i) = CountryMask.CountryArray(i).bioenergy_land_use_GCAM_2050_per_cap_per_SSP(4);
end


continents = {CountryMask.CountryArray.GPW_continent_string};
country_name = {CountryMask.CountryArray.country_name};

unique_continents = unique(continents);

figure
for i = 1:4
    


    subplot(2,2,i)
end

end

