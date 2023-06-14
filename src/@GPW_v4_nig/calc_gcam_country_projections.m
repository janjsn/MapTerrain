function obj = calc_gcam_country_projections(obj)
%CALC_GCAM_COUNTRY_PROJECTIONS Summary of this function goes here

n_countries = length(obj.CountryArray);

years = [2015:5:2100];
n_years = length(years);
ssps = {'ssp1', 'ssp2', 'ssp4', 'ssp5'};
ssps_int = [1 2 4 5];
rcp = 2.6;


for cntr = 1:n_countries
    obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp1_rcp26_vec = zeros(1,n_years);
    obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp2_rcp26_vec = zeros(1,n_years);
    obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp4_rcp26_vec = zeros(1,n_years);
    obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp5_rcp26_vec = zeros(1,n_years);
end

filename_base = 'Data/GCAM/GCAM_Demeter_LU_';

for yr = 1:length(years)
    
    for scens = 1:length(ssps)
       fprintf([num2str(years(yr)) ': ' ssps{scens} '\n']);
        GCAM_this = GCAM([filename_base ssps{scens} '_rcp26_modelmean_' num2str(years(yr)) '.nc'], years(yr), [ssps{scens} '_' num2str(rcp)], ssps_int, rcp);



        for cntr = 1:n_countries
            mask_this = obj.countryMask_05deg == obj.CountryArray(cntr).GPW_country_ISO_numeric;

            be_lu_this = sum(sum(GCAM_this.bioenergy_crops_rf_hectare(mask_this))) + sum(sum(GCAM_this.bioenergy_crops_ir_hectare(mask_this)));

            if ssps_int(scens) == 1
                obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp1_rcp26_vec(yr) = be_lu_this;
            elseif ssps_int(scens) == 2
                obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp2_rcp26_vec(yr) = be_lu_this;
            elseif ssps_int(scens) == 4
                obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp4_rcp26_vec(yr) = be_lu_this;
            elseif ssps_int(scens) == 5
                obj.CountryArray(cntr).bioenergy_land_use_GCAM_ssp5_rcp26_vec(yr) = be_lu_this;
            end

        end


    end
end
end
