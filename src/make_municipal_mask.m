[~,~,in_data] = xlsread('Data/municipality_mask.xlsx');

id_col = 11;
name_col = 8;


n_lat = 2160;
n_lon = 4320;
%Dimension steps
lat_step  = 180/n_lat;
lon_step = 360/n_lon;

%Dimension vectors
lat = [90-(lat_step/2):-lat_step:-90+(lat_step/2)];
lon = [-180+(lon_step/2):lon_step:180-(lon_step/2)];

municipal_mask = zeros(n_lon,n_lat);
fylke_mask= zeros(n_lon,n_lat);

lat_points = [in_data{2:end,lat_col}];
lon_points = [in_data{2:end,lon_col}];
kommune_id = [in_data{2:end,id_col}];
kommune_name = {in_data{2:end,name_col}};

idx_lat_points = zeros(1,length(lat_points));
idx_lon_points = zeros(1,length(lon_points));

for i = 1:length(idx_lat_points)
[~,idx_lat_points(i)] = min(abs(lat_points(i)-lat)); 
[~,idx_lon_points(i)] = min(abs(lon_points(i)-lon));

municipal_mask(idx_lon_points(i),idx_lat_points(i)) = kommune_id(i);

fylke_mask(idx_lon_points(i),idx_lat_points(i)) = floor((kommune_id(i)-mod(kommune_id(i),100))/100);


end

filename = 'municipal_ids.nc';

if exist(filename, 'file')
    delete(filename)
end

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

nccreate(filename, 'municipal_id', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'municipal_id', 'standard_name', 'municipal_id');
ncwriteatt(filename, 'municipal_id', 'long_name', 'municipal_id');
ncwriteatt(filename, 'municipal_id', 'units', 'id');
ncwriteatt(filename, 'municipal_id', 'missing_value', '-999');

nccreate(filename, 'fylke_id', 'Dimensions', {'lon' length(lon) 'lat' length(lat)}, 'DeflateLevel', 4);
ncwriteatt(filename, 'fylke_id', 'standard_name', 'fylke_id');
ncwriteatt(filename, 'fylke_id', 'long_name', 'fylke_id');
ncwriteatt(filename, 'fylke_id', 'units', 'id');
ncwriteatt(filename, 'fylke_id', 'missing_value', '-999');

ncwrite(filename, 'lat', lat);
ncwrite(filename, 'lon', lon);
ncwrite(filename, 'municipal_id', municipal_mask);
ncwrite(filename, 'fylke_id', fylke_mask);