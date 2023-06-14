

[ac_1deg, lat_1deg, lon_1deg, ~, ~]  = aggregateMatrix2givenDimensions( old_matrix,old_longitudeVector_centered,old_latitudeVector_centered, new_nRows_longitude, new_nCols_latitude );


    lat = CountryMask.lat_5arcmin;
    lon = CountryMask.lon_5arcmin;

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

    nccreate(filename, 'cropland_fractions', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'cropland_fractions', 'standard_name', 'cropland_fractions');
    ncwriteatt(filename, 'cropland_fractions', 'long_name', 'cropland_fractions');
    ncwriteatt(filename, 'cropland_fractions', 'units', 'fractions');
    ncwriteatt(filename, 'cropland_fractions', 'missing_value', '-999');

    nccreate(filename, 'abandoned_cropland_fractions', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'standard_name', 'abandoned_cropland_fractions');
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'long_name', 'abandoned_cropland_fractions');
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'units', 'fractions');
    ncwriteatt(filename, 'abandoned_cropland_fractions', 'missing_value', '-999');
