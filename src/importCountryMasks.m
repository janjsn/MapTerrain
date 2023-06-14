function [ CountryMask ] = importCountryMasks(  )
%IMPORTCOUNTRYMASKS Summary of this function goes here
%   Detailed explanation goes here
filename_country_masks = 'Data/gpw-v4-national-identifier-grid-rev11_2pt5_min_tif/gpw_v4_national_identifier_grid_rev11_2pt5_min.tif';


[countryMasks_2p5_arcmin, GeoInfo_countryMask] = getCountryMasks_gpw(filename_country_masks);
countryMasks_2p5_arcmin = countryMasks_2p5_arcmin';

CountryArray = getCountryMaskIDs( 'Data/gpw-v4-population-count-rev11_2020_2pt5_min_tif/gpw-v4-country-level-summary-rev11.xlsx', 'GPWv4 Rev11 Summary' );
    
CountryMask = GPW_v4_nig(CountryArray, countryMasks_2p5_arcmin, GeoInfo_countryMask);
CountryMask = CountryMask.generate_5arcmin_mask;
end






