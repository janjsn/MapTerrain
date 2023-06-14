classdef GPW_v4_nig
    %GPW_V4_NIG Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        CountryArray
        countryMask
        GeoRef_countryMask
        latitudeVector_mask_centered
        longitudeVector_mask_centered
        latitude_bounds
        longitude_bounds
        cellAreaPerLatitude_hectare
        waterValue = 32767;
        countryMask_5arcmin
        lat_5arcmin
        lon_5arcmin

        countryMask_05deg
        lat_005deg
        lon_005deg

    end
    
    methods
        
        function obj = GPW_v4_nig(CountryArray, countryMask, GeoRef)
           obj.CountryArray = CountryArray;
           obj.countryMask = countryMask;
           obj.GeoRef_countryMask = GeoRef;
           %dimensions, center
           obj.latitudeVector_mask_centered = [90-GeoRef.CellExtentInLatitude/2:-GeoRef.CellExtentInLatitude:-90+GeoRef.CellExtentInLatitude/2];
           obj.longitudeVector_mask_centered = [-180+GeoRef.CellExtentInLatitude/2:GeoRef.CellExtentInLatitude:180-GeoRef.CellExtentInLatitude/2];
           %bounds
           obj.latitude_bounds = zeros(2,length(obj.latitudeVector_mask_centered));
           obj.latitude_bounds(1,:) = [90:-GeoRef.CellExtentInLatitude:-90+GeoRef.CellExtentInLatitude];
           obj.latitude_bounds(2,:) = [90-GeoRef.CellExtentInLatitude:-GeoRef.CellExtentInLatitude:-90];
           obj.longitude_bounds = zeros(2,length(obj.longitudeVector_mask_centered));
           obj.longitude_bounds(1,:) = [-180:GeoRef.CellExtentInLongitude:180-GeoRef.CellExtentInLongitude];
           obj.longitude_bounds(2,:) = [-180+GeoRef.CellExtentInLongitude:GeoRef.CellExtentInLongitude:180];
           %cell area per latitude
           Earth = referenceSphere('earth','m');
           surfaceArea_hectare = Earth.SurfaceArea/10000;
                      
           obj.cellAreaPerLatitude_hectare = zeros(1,length(obj.latitudeVector_mask_centered));
               for i = 1:length(obj.cellAreaPerLatitude_hectare)
                   obj.cellAreaPerLatitude_hectare(i) =  areaquad(obj.latitude_bounds(1,i),0,obj.latitude_bounds(2,i),GeoRef.CellExtentInLongitude)*surfaceArea_hectare;           
               end
           
        end

        function obj = generate_5arcmin_mask(obj)
            step = 360/4320;
            lat_5arcmin = [90-step/2:-step:-90+step/2];
            lon_5arcmin = [-180+step/2:step:180-step/2];
            mask_5arcmin = zeros(length(lon_5arcmin), length(lat_5arcmin));

            idx_lat = zeros(1,length(lat_5arcmin));
            idx_lon = zeros(1,length(lon_5arcmin));

            for i = 1:length(lat_5arcmin)
                [~,idx_lat(i)] = min(abs(lat_5arcmin(i)-obj.latitudeVector_mask_centered));
            end
            
            for j = 1:length(lon_5arcmin)
                [~,idx_lon(j)] = min(abs(lon_5arcmin(j)-obj.longitudeVector_mask_centered));
            end
            
            for i = 1:length(lat_5arcmin)
                for j = 1:length(lon_5arcmin)
                    mask_5arcmin(j,i) = obj.countryMask(idx_lon(j),idx_lat(i));
                end
            end

            obj.countryMask_5arcmin = mask_5arcmin;
            obj.lat_5arcmin = lat_5arcmin;
            obj.lon_5arcmin = lon_5arcmin;

            %% 0.05 DEGREES
            step = 0.05
            lat_005 = [90-step/2:-step:-90+step/2];
            lon_005 = [-180+step/2:step:180-step/2];
            mask_005 = zeros(length(lon_005),length(lat_005));

            idx_lat = zeros(1,length(lat_005));
            idx_lon = zeros(1,length(lon_005));

            for i = 1:length(lat_005)
                [~,idx_lat(i)] = min(abs(lat_005(i)-obj.latitudeVector_mask_centered));
            end
            
            for j = 1:length(lon_005)
                [~,idx_lon(j)] = min(abs(lon_005(j)-obj.longitudeVector_mask_centered));
            end
            
            for i = 1:length(lat_005)
                for j = 1:length(lon_005)
                    mask_005(j,i) = obj.countryMask(idx_lon(j),idx_lat(i));
                end
            end
               
            obj.countryMask_05deg = mask_005;
            obj.lat_005deg = lat_005;
            obj.lon_005deg = lon_005;


        end

        obj = calc_gcam_country_projections(obj)
        plot_gcam_be_proj_globe_nor(obj)



    end
    
end

