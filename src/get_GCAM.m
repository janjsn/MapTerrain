function GCAM_array = get_GCAM(  )

n_files = 8;

GCAM_array(1:n_files) = GCAM;

GCAM_array(1) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp1_rcp26_modelmean_2050.nc', 2050, 'SSP1-26, 2050', 1, 2.6);
GCAM_array(2) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp1_rcp26_modelmean_2100.nc', 2100, 'SSP1-26, 2100', 1, 2.6);

GCAM_array(3) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp2_rcp26_modelmean_2050.nc', 2050, 'SSP2-26, 2050', 2, 2.6);
GCAM_array(4) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp2_rcp26_modelmean_2100.nc', 2100, 'SSP2-26, 2100', 2, 2.6);

GCAM_array(5) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp4_rcp26_modelmean_2050.nc', 2050, 'SSP4-26, 2050', 4, 2.6);
GCAM_array(6) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp4_rcp26_modelmean_2100.nc', 2100, 'SSP4-26, 2100', 4, 2.6);

GCAM_array(7) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp5_rcp26_modelmean_2050.nc', 2050, 'SSP5-26, 2050', 5, 2.6);
GCAM_array(8) = GCAM('Data/GCAM/GCAM_Demeter_LU_H_ssp5_rcp26_modelmean_2100.nc', 2100, 'SSP5-26, 2100', 5, 2.6);


end

