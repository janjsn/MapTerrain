SSPs = {'ssp1', 'ssp2', 'ssp4', 'ssp5'};
years = [2015:5:2100];
rcp = {'rcp26'};
filename_base = 'Data/GCAM/GCAM_Demeter_LU_';

for scens = 1:length(SSPs)

    for yr = 1:length(years)
        filename_this = [filename_base SSPs{scens} '_modelmean' ];
    end
end