function plot_ac_per_cap_countries(CountryMask)


ac_per_cap = [CountryMask.CountryArray.ac_per_cap];
continents = {CountryMask.CountryArray.GPW_continent_string};
%country_name = {CountryMask.CountryArray.country_name};

unique_continents = unique(continents);

intervals = [0 0.02 0.04 0.06 0.08 0.10
    0.02 0.04 0.06 0.08 0.10 inf];

%interval_vec_1 = [0:0.00]

mSize = size(intervals);

n_countries_per_interval = zeros(1,mSize(2));
x_vec = [1:1:mSize(2)];

for i = 1:mSize(2)
    n_countries_per_interval(i) = sum(ac_per_cap >= intervals(1,i) & ac_per_cap < intervals(2,i));
end


figure
bar(x_vec,n_countries_per_interval)


plot_matrix = zeros(mSize(2),length(unique_continents));

for i = 1:length(unique_continents)
    binary = strcmp(unique_continents{i},continents);
    ac_per_cap_this  = ac_per_cap(binary);

    for j = 1:mSize(2)
        plot_matrix(j,i) = sum(ac_per_cap_this >= intervals(1,j) & ac_per_cap_this < intervals(2,j));
    end

end

figure
bar(plot_matrix, 'stacked');
legend(unique_continents);

end

