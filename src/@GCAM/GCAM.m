classdef GCAM
    %GCAM Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        lat
        lon
        filename
        description
        year
        SSP
        RCP
        bioenergy_crops_rf_hectare
        bioenergy_crops_ir_hectare
        bioenergy_crops_rf_fractions
        bioenergy_crops_ir_fractions
        
        bioenergy_crops_rf_hectare_tot
        bioenergy_crops_ir_hectare_tot
        
    end
    
    methods
        function obj = GCAM(filename, year, description, SSP_int, RCP_int)
            if nargin > 0
                obj.year = year;
                obj.description = description;
                obj.SSP = SSP_int;
                obj.RCP = RCP_int;
                
                
                ncid = netcdf.open(filename);
                
                [~,nvars,~,~] = netcdf.inq(ncid);
                
                for i = 1:nvars
                    varname = netcdf.inqVar(ncid,i-1);
                    if strcmp('latitude', varname)
                        obj.lat = double(netcdf.getVar(ncid,i-1));
                    elseif strcmp('longitude', varname)
                        obj.lon = double(netcdf.getVar(ncid,i-1));
                    elseif strcmp('PFT29', varname)
                        obj.bioenergy_crops_rf_fractions = double(netcdf.getVar(ncid,i-1)')/100;
                        obj.bioenergy_crops_rf_fractions(isnan(obj.bioenergy_crops_rf_fractions)) = 0;
                    elseif strcmp('PFT30', varname)
                        obj.bioenergy_crops_ir_fractions = double(netcdf.getVar(ncid,i-1)')/100;
                        obj.bioenergy_crops_ir_fractions(isnan(obj.bioenergy_crops_ir_fractions)) = 0;
                    end
                end
                
                lat_step = abs(obj.lat(1)-obj.lat(2));
                
                lat_bnds = zeros(2,length(obj.lat));
                lat_bnds(1,:) = [obj.lat+lat_step/2];
                lat_bnds(2,:) = [obj.lat-lat_step/2];
                
                cell_extent_in_longitude = abs(obj.lon(1)-obj.lon(2));
                
                [ areaPerLatitude_hectare ] = get_cell_area_per_latitude( lat_bnds, cell_extent_in_longitude);
                
                obj.bioenergy_crops_rf_hectare = zeros(length(obj.lon), length(obj.lat));
                obj.bioenergy_crops_ir_hectare = zeros(length(obj.lon), length(obj.lat));
                
                for i = 1:length(obj.lat)
                    obj.bioenergy_crops_rf_hectare(:,i) = obj.bioenergy_crops_rf_fractions(:,i)*areaPerLatitude_hectare(i);
                    obj.bioenergy_crops_ir_hectare(:,i) = obj.bioenergy_crops_ir_fractions(:,i)*areaPerLatitude_hectare(i);
                end
               
                obj.bioenergy_crops_rf_hectare_tot = 10^-6*sum(sum(obj.bioenergy_crops_rf_hectare));
                obj.bioenergy_crops_ir_hectare_tot = 10^-6*sum(sum(obj.bioenergy_crops_ir_hectare));
                
            end
            
        end
    end
    
end

