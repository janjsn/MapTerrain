function plot_pe_potential_per_region()

pe_Trondelag = 3.10;
pe_More = 1.14;
pe_Vestland = 1.02;
pe_Nordland = 0.87;
pe_Troms = 0.49;
pe_other = 0.39;

figure
bar([pe_Trondelag pe_More pe_Vestland pe_Nordland pe_Troms pe_other]);
xticklabels({'Trøndelag', 'Møre og Romsdal', 'Vestland', 'Nordland', 'Troms & Finnmark', 'Other'});
ylabel('PJ year^{-1}')

end

