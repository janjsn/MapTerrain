function plot_gcam_be_proj_globe_nor(obj)

time = [2015:5:2100];

global_ssp1_proj = zeros(1,length(time));
global_ssp2_proj = zeros(1,length(time));
global_ssp4_proj = zeros(1,length(time));
global_ssp5_proj = zeros(1,length(time));

unique_continents = unique({obj.CountryArray.GPW_continent_string});
n_continents = length(unique_continents);

continental_ssp1_proj = zeros(n_continents,length(time));
continental_ssp2_proj = zeros(n_continents,length(time));
continental_ssp4_proj = zeros(n_continents,length(time));
continental_ssp5_proj = zeros(n_continents,length(time));


for i = 1:length(obj.CountryArray)
    global_ssp1_proj = global_ssp1_proj + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
    global_ssp2_proj = global_ssp2_proj + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp2_rcp26_vec;
    global_ssp4_proj = global_ssp4_proj + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp4_rcp26_vec;
    global_ssp5_proj = global_ssp5_proj + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp5_rcp26_vec;

    for j = 1:length(unique_continents)
        if strcmp(unique_continents{j}, obj.CountryArray(i).GPW_continent_string)
            continental_ssp1_proj(j,:) = continental_ssp1_proj(j,:) + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
            continental_ssp2_proj(j,:) = continental_ssp1_proj(j,:) + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
            continental_ssp4_proj(j,:) = continental_ssp1_proj(j,:) + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
            continental_ssp5_proj(j,:) = continental_ssp1_proj(j,:) + obj.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
        end


    end

    if strcmp(obj.CountryArray(i).country_name, 'Norway')
        nor_ssp1_proj = obj.CountryArray(i).bioenergy_land_use_GCAM_ssp1_rcp26_vec;
        nor_ssp2_proj = obj.CountryArray(i).bioenergy_land_use_GCAM_ssp2_rcp26_vec;
        nor_ssp4_proj = obj.CountryArray(i).bioenergy_land_use_GCAM_ssp4_rcp26_vec;
        nor_ssp5_proj = obj.CountryArray(i).bioenergy_land_use_GCAM_ssp5_rcp26_vec;

        nor_fao_agricultural_land_2020(1:length(time)) = 0.986;
        nor_fao_land_under_temperorary_crops_2020(1:length(time)) = 0.328;

        figure
        plot(time, nor_ssp1_proj*10^-6,'LineWidth',2);
        hold on
        plot(time, nor_ssp2_proj*10^-6,'LineWidth',2),;
        plot(time, nor_ssp4_proj*10^-6,'LineWidth',2);
        plot(time, nor_ssp5_proj*10^-6,'LineWidth',2);
        plot(time, nor_fao_agricultural_land_2020, '--','LineWidth',2)
        plot(time, nor_fao_land_under_temperorary_crops_2020, '--','LineWidth',2)
        ylabel('Mha')
        xlabel('Year');
        xlim([time(1) time(end)]);
        legend({'SSP1-2.6','SSP2-2.6','SSP4-2.6','SSP5-2.6', 'Agricultural land use 2020', 'Annual crops 2020'})

        filename = 'Global_GCAM_be_lu.pdf';
        print('-vector','-dpdf', '-r1000', filename)

    end

end

%Plot global

global_fao_agricultural_land_2020(1:length(time)) = 4745;
global_fao_land_under_temperorary_crops_2020(1:length(time)) = 1562;

figure
plot(time, global_ssp1_proj*10^-6,'LineWidth',2);
hold on
plot(time, global_ssp2_proj*10^-6,'LineWidth',2);
plot(time, global_ssp4_proj*10^-6,'LineWidth',2);
plot(time, global_ssp5_proj*10^-6,'LineWidth',2);
ylabel('Mha');
xlabel('Year');
plot(time, global_fao_agricultural_land_2020, '--','LineWidth',2)
plot(time, global_fao_land_under_temperorary_crops_2020, '--','LineWidth',2)
xlim([time(1) time(end)]);


ylabel('Mha')
xlabel('Year');
legend({'SSP1-2.6','SSP2-2.6','SSP4-2.6','SSP5-2.6', 'Agricultural LU 2020', 'Annual crops LU 2020'})

filename = 'Global_GCAM_be_lu.pdf';
print('-vector','-dpdf', '-r1000', filename)
end

