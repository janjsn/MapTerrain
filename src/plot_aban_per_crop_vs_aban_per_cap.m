function plot_aban_per_crop_vs_aban_per_cap(CountryMask)


ac_per_cap = [CountryMask.CountryArray.ac_per_cap];

ac_per_crop = [CountryMask.CountryArray.abandoned_cropland_as_share_of_2020_cropland];

continents = {CountryMask.CountryArray.GPW_continent_string};

unique_continents = unique(continents);

figure
hold on
grid on

n_years = 28-3;

continent_ac = zeros(1,length(unique_continents));
continent_cap = zeros(1,length(unique_continents));
continent_cropland = zeros(1,length(unique_continents));


for i = 1:length(unique_continents)
    binary = strcmp(unique_continents{i},continents);
    scatter(100*ac_per_crop(binary),ac_per_cap(binary)/n_years )

    continent_ac(i) = sum([CountryMask.CountryArray(binary).abandoned_cropland_1992_to_2020_tot]);
    continent_cap(i) = sum([CountryMask.CountryArray(binary).population_2020]);
    continent_cropland(i) = sum([CountryMask.CountryArray(binary).cropland_2020_tot]);

end



legend(unique_continents);

%ylim([0 0.004])
xlim([0 1])
ylabel('Cropland abandonment (1992-2020) as share of cropland (2020)')
xlabel('Cropland abandonment per capita per year')

plot_continents = 0;
if plot_continents == 1

continent_ac_per_cap = continent_ac./continent_cap;
continent_ac_per_cap_year = continent_ac./(continent_cap*n_years);
continent_ac_per_cropland2020 = continent_ac./continent_cropland;
continent_ac_per_cropland2020_year = continent_ac./(continent_cropland*n_years);

figure
hold on
colors = {'b','r','g','y','k', 'c'}; 
legend(unique_continents)

for i = 1:length(continent_ac_per_cropland2020)
viscircles([continent_ac_per_cap(i) 100*continent_ac_per_cropland2020(i)], continent_ac(i)*10^-10,'Color',colors{i});
end
end


end

