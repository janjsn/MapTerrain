function upscale_and_export_maps_1deg(obj)


filename = 'Output/output_1deg.nc';

if exist(filename, 'file')
delete(filename)
end

step = 360/4320;
lon_5arcmin = [-180+step/2:step:180-step/2];
lat_5arcmin = [90-step/2:-step:-90+step/2];

[land_1deg, lat_1deg, lon_1deg, lat_bnds_1deg,  ~]  = aggregateMatrix2givenDimensions( obj.land_availability,lon_5arcmin,lat_5arcmin, 360, 180 );


[ area_after_lat_1deg ] = get_cell_area_per_latitude( lat_bnds_1deg, abs(lon_1deg(1)-lon_1deg(2)));

land_fraction_1deg = zeros(length(lon_1deg),length(lat_1deg));


for lats = 1:length(lat_1deg)
    land_fraction_1deg(:,lats) = land_1deg(:,lats)/area_after_lat_1deg(lats);
end

fprintf('Primary energy.. \n')
% Primary energy
[pe_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( obj.pe,lon_5arcmin,lat_5arcmin, 360, 180 );
fprintf('Final energy.. \n')

% Final energy
[fe_bioelectricity_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( obj.fe_bioelectricity,lon_5arcmin,lat_5arcmin, 360, 180 );
[fe_bioelectricity_ccs_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( obj.fe_bioelectricity_with_css_GJ_per_year,lon_5arcmin,lat_5arcmin, 360, 180 );
[fe_ft_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( obj.fe_FT_diesel,lon_5arcmin,lat_5arcmin, 360, 180 );
[fe_ft_ccs_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( obj.fe_FT_diesel_with_ccs_GJ_per_year,lon_5arcmin,lat_5arcmin, 360, 180 );

% Emission factors

%Emissions (kg CO2eq)
fprintf('Emission factors .. \n')
emissions_electricity_5arcmin = obj.emission_factor_electricity.*obj.fe_elec_after_aban_year;
emissions_electricity_ccs_5arcmin = obj.emission_factor_electricity_CCS.*obj.fe_elec_ccs_after_aban_year;
emissions_ft_5arcmin = obj.emission_factor_FT.*obj.fe_ft_after_aban_year;
emissions_ft_ccs_5arcmin = obj.emission_factor_FT_CCS.*obj.fe_ft_ccs_after_aban_year;

emissions_electricity_5arcmin(isnan(emissions_electricity_5arcmin)) = 0;
emissions_electricity_ccs_5arcmin(isnan(emissions_electricity_ccs_5arcmin)) = 0;
emissions_ft_5arcmin(isnan(emissions_ft_5arcmin)) = 0;
emissions_ft_ccs_5arcmin(isnan(emissions_ft_ccs_5arcmin)) =0 ;

% Make 2D
emissions_electricity_5arcmin = sum(emissions_electricity_5arcmin,3);
emissions_electricity_ccs_5arcmin = sum(emissions_electricity_ccs_5arcmin,3);
emissions_ft_5arcmin = sum(emissions_ft_5arcmin,3);
emissions_ft_ccs_5arcmin = sum(emissions_ft_ccs_5arcmin,3);

[emissions_electricity_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( emissions_electricity_5arcmin,lon_5arcmin,lat_5arcmin, 360, 180 );
[emissions_electricity_ccs_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( emissions_electricity_ccs_5arcmin,lon_5arcmin,lat_5arcmin, 360, 180 );

[emissions_ft_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( emissions_ft_5arcmin,lon_5arcmin,lat_5arcmin, 360, 180 );
[emissions_ft_ccs_1deg, ~, ~, ~,  ~]  = aggregateMatrix2givenDimensions( emissions_ft_ccs_5arcmin,lon_5arcmin,lat_5arcmin, 360, 180 );

emission_factor_electricity_1deg = emissions_electricity_1deg./fe_bioelectricity_1deg;
emission_factor_electricity_ccs_1deg = emissions_electricity_ccs_1deg./fe_bioelectricity_ccs_1deg;
emission_factor_ft_1deg = emissions_ft_1deg./fe_ft_1deg;
emission_factor_ccs_1deg = emissions_ft_ccs_1deg./fe_ft_ccs_1deg;



%% EXPORT

lat = lat_1deg;
lon = lon_1deg;

nccreate(filename,'lat','Dimensions',{'lat' length(lat)});
ncwriteatt(filename, 'lat', 'standard_name', 'latitude');
ncwriteatt(filename, 'lat', 'long_name', 'latitude');
ncwriteatt(filename, 'lat', 'units', 'degrees_north');
ncwriteatt(filename, 'lat', '_CoordinateAxisType', 'Lat');
ncwriteatt(filename, 'lat', 'axis', 'Y');

nccreate(filename,'lon','Dimensions',{'lon' length(lon)});
ncwriteatt(filename, 'lon', 'standard_name', 'longitude');
ncwriteatt(filename, 'lon', 'long_name', 'longitude');
ncwriteatt(filename, 'lon', 'units', 'degrees_east');
ncwriteatt(filename, 'lon', '_CoordinateAxisType', 'Lon');
ncwriteatt(filename, 'lon', 'axis', 'X');


nccreate(filename, 'abandoned_cropland_fractions', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'abandoned_cropland_fractions', 'standard_name', 'abandoned_cropland_fractions');
ncwriteatt(filename, 'abandoned_cropland_fractions', 'long_name', 'abandoned_cropland_fractions');
ncwriteatt(filename, 'abandoned_cropland_fractions', 'units', 'fractions');
ncwriteatt(filename, 'abandoned_cropland_fractions', 'missing_value', '-999');

nccreate(filename, 'primary_energy', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy', 'standard_name', 'primary_energy');
ncwriteatt(filename, 'primary_energy', 'long_name', 'primary_energy');
ncwriteatt(filename, 'primary_energy', 'units', 'TJ year-1');
ncwriteatt(filename, 'primary_energy', 'missing_value', '-999');

nccreate(filename, 'primary_energy_yield', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'primary_energy_yield', 'standard_name', 'primary_energy_yield');
ncwriteatt(filename, 'primary_energy_yield', 'long_name', 'primary_energy_yield');
ncwriteatt(filename, 'primary_energy_yield', 'units', 'TJ year-1');
ncwriteatt(filename, 'primary_energy_yield', 'missing_value', '-999');

nccreate(filename, 'final_energy_bioelectricity', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_bioelectricity', 'standard_name', 'final_energy_bioelectricity');
ncwriteatt(filename, 'final_energy_bioelectricity', 'long_name', 'final_energy_bioelectricity');
ncwriteatt(filename, 'final_energy_bioelectricity', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_bioelectricity', 'missing_value', '-999');

nccreate(filename, 'final_energy_bioelectricity_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'standard_name', 'final_energy_bioelectricity_ccs');
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'long_name', 'final_energy_bioelectricity_ccs');
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_bioelectricity_ccs', 'missing_value', '-999');

nccreate(filename, 'final_energy_fischer_tropsch_diesel', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'standard_name', 'final_energy_fischer_tropsch_diesel');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'long_name', 'final_energy_fischer_tropsch_diesel');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel', 'missing_value', '-999');

nccreate(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'standard_name', 'final_energy_fischer_tropsch_diesel_ccs');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'long_name', 'final_energy_fischer_tropsch_diesel_ccs');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'units', 'TJ year-1');
ncwriteatt(filename, 'final_energy_fischer_tropsch_diesel_ccs', 'missing_value', '-999');

nccreate(filename, 'emission_factor_bioelectricity', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'emission_factor_bioelectricity', 'standard_name', 'emission_factor_bioelectricity');
ncwriteatt(filename, 'emission_factor_bioelectricity', 'long_name', 'emission_factor_bioelectricity');
ncwriteatt(filename, 'emission_factor_bioelectricity', 'units', 'kgCO2eq GJ-1');
ncwriteatt(filename, 'emission_factor_bioelectricity', 'missing_value', '-999');

nccreate(filename, 'emission_factor_bioelectricity_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'standard_name', 'emission_factor_bioelectricity_ccs');
ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'long_name', 'emission_factor_bioelectricity_ccs');
ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'units', 'kgCO2eq GJ-1');
ncwriteatt(filename, 'emission_factor_bioelectricity_ccs', 'missing_value', '-999');

nccreate(filename, 'emission_factor_fischer_tropsch_diesel', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'standard_name', 'emission_factor_fischer_tropsch_diesel');
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'long_name', 'emission_factor_fischer_tropsch_diesel');
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'units', 'kgCO2eq GJ-1');
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel', 'missing_value', '-999');

nccreate(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'standard_name', 'emission_factor_fischer_tropsch_diesel_ccs');
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'long_name', 'emission_factor_fischer_tropsch_diesel_ccs');
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'units', 'kgCO2eq GJ-1');
ncwriteatt(filename, 'emission_factor_fischer_tropsch_diesel_ccs', 'missing_value', '-999');

ncwrite(filename, 'lat', lat);
ncwrite(filename, 'lon', lon);
ncwrite(filename, 'abandoned_cropland_fractions', land_fraction_1deg);
ncwrite(filename, 'primary_energy', 10^-3*pe_1deg);
ncwrite(filename, 'primary_energy_yield', obj.pe_yield_1deg_2d);
ncwrite(filename, 'final_energy_bioelectricity', 10^-3*fe_bioelectricity_1deg);
ncwrite(filename, 'final_energy_bioelectricity_ccs', 10^-3*fe_bioelectricity_ccs_1deg);
ncwrite(filename, 'final_energy_fischer_tropsch_diesel', 10^-3*fe_ft_1deg);
ncwrite(filename, 'final_energy_fischer_tropsch_diesel_ccs', 10^-3*fe_ft_ccs_1deg);
ncwrite(filename, 'emission_factor_bioelectricity', emission_factor_electricity_1deg);
ncwrite(filename, 'emission_factor_bioelectricity_ccs', emission_factor_electricity_ccs_1deg);
ncwrite(filename, 'emission_factor_fischer_tropsch_diesel',emission_factor_ft_1deg);
ncwrite(filename, 'emission_factor_fischer_tropsch_diesel_ccs', emission_factor_ccs_1deg);



end

