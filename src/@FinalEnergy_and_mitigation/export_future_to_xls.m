function export_future_to_xls(obj)

FutureArray = obj.FutureArray;

%% SET FILE NAME
outfile = 'Output/Future_mitigation.xlsx';

%Delete if existing
if exist(outfile, 'file')
    delete(outfile);
end


n_variables = 5;

first_num_row = 3;

header_l1 = {'year', 'time_vec', 'land_regrowing', 'land_food_prod', 'land_bioenergy', 'land_beccs_electricity', 'land_beccs_FT', ...
    'final_energy_beccs_tot','final_energy_beccs_bioelectricity', 'final_energy_beccs_FT_diesel'...
    'CO2_flux_sequestration_natural_regrowth', 'CO2_flux_recultivation_loss_natural_regrowth', 'CO2_flux_net_aboveground_carbon', ...
    'CO2_flux_BECCS_bioelectricity', 'CO2_flux_BECCS_FT', 'CO2_flux_BECCS_net',...
    'CO2_flux_NR_BECCS_net', ...
    'CO2_net_NR_accumulated', 'CO2_net_BECCS_accumulated','CO2_net_accumulated'};

header_l2 = {'yr', 'yr', 'Mha', 'Mha', 'Mha', 'Mha', 'Mha',...
    'EJ yr-1', 'EJ yr-1', 'EJ yr-1'...
    'GtCO2eq yr-1', 'GtCO2eq yr-1', 'GtCO2eq yr-1', ...
    'GtCO2eq yr-1',  'GtCO2eq yr-1',  'GtCO2eq yr-1', ...
    'GtCO2eq yr-1' , ...
    'GtCO2eq' ,  'GtCO2eq' ,  'GtCO2eq' };





for i = 1:length(FutureArray)
    outMatrix = cell(length(FutureArray(i).time)+2,n_variables);

    outMatrix(first_num_row:end,1) = num2cell(FutureArray(i).time);
    outMatrix(first_num_row:end,2) = num2cell(FutureArray(i).time_since_t0);

    outMatrix(first_num_row:end,3) = num2cell(FutureArray(i).totals_vec_land_regrowing);
    outMatrix(first_num_row:end,4) = num2cell(FutureArray(i).totals_vec_land_converted_to_food_production);
    outMatrix(first_num_row:end,5) = num2cell(FutureArray(i).totals_vec_land_converted_to_bioenergy);
    outMatrix(first_num_row:end,6) = num2cell(FutureArray(i).totals_vec_land_converted_to_beccs_electricity);
    outMatrix(first_num_row:end,7) = num2cell(FutureArray(i).totals_vec_land_converted_to_beccs_ft);

    outMatrix(first_num_row:end,8) = num2cell(FutureArray(i).totals_vec_fe_ccs_tot);
    outMatrix(first_num_row:end,9) = num2cell(FutureArray(i).totals_vec_fe_bioelectricity_ccs);
    outMatrix(first_num_row:end,10) = num2cell(FutureArray(i).totals_vec_fe_ft_ccs);

    outMatrix(first_num_row:end,11) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_natural_regrowth_added);
    outMatrix(first_num_row:end,12) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_natural_regrowth_removed);
    outMatrix(first_num_row:end,13) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_natural_regrowth_net);

    outMatrix(first_num_row:end,14) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_beccs_electricity);
    outMatrix(first_num_row:end,15) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_beccs_ft);
    outMatrix(first_num_row:end,16) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_beccs_tot);

    outMatrix(first_num_row:end,17) = num2cell(FutureArray(i).totals_vec_tco2eq_flux_net_nr_beccs);

    outMatrix(first_num_row:end,18) = num2cell(FutureArray(i).cumulative_tco2_no_substitution_natural_regrowth_tot);
    outMatrix(first_num_row:end,19) = num2cell(FutureArray(i).cumulative_tco2_no_substitution_beccs_tot);
    outMatrix(first_num_row:end,20) = num2cell(FutureArray(i).cumulative_tco2_no_substitution_tot);


    outMatrix(1,:) = header_l1;
    outMatrix(2,:) = header_l2;

    writecell(outMatrix,outfile, 'Sheet', FutureArray(i).scenario_description)


end


for i = 1:length(FutureArray)
    figure

    plot(FutureArray(i).time(2:end),FutureArray(i).totals_vec_tco2eq_flux_natural_regrowth_net(2:end),'LineWidth',2);
    hold on
    plot(FutureArray(i).time(2:end), FutureArray(i).totals_vec_tco2eq_flux_beccs_tot(2:end),'LineWidth',2);
    plot(FutureArray(i).time(2:end), FutureArray(i).totals_vec_tco2eq_flux_net_nr_beccs(2:end),'LineWidth',2);
    xlabel('Year');
    ylabel('GtCO_{2}eq yr^{-1}');
    legend({'Natural vegetation', 'BECCS', 'Net'})

    filename = ['Output/Future/' FutureArray(i).scenario_description '.pdf'];
    if exist(filename, 'file')
        delete(filename);
    end

    print('-vector','-dpdf', '-r1000', filename)
    %save('Output/src_data_ft_diesel_ccs_Norway.mat', 'fe_ft_ccs_potential_Nor', 'x_vec_emission_factor');

end

idx_of_interest = 0;
c=1;

color_codes = zeros(6,3);

for i = 1:length(FutureArray)
    if strcmp(FutureArray(i).scenario_description, 'R23_F100_EL0_FT0')
        idx_of_interest(c) = i;

        idx_baseline = i;
        color_codes(c,:) = [0.4940 0.1840 0.5560]; % Purple
        legend_labels{c} = 'Baseline';

        c=c+1;
    elseif strcmp(FutureArray(i).scenario_description, 'Rinf_F100_EL0_FT0')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Continued regrowth';

        color_codes(c,:) = [34 139 34]/255	;	% GREEN
        idx_nat_reg = i;
        c=c+1;
    elseif strcmp(FutureArray(i).scenario_description, 'R23_F0_EL100_FT0')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Recultivation, BECCS bioelectricity';
        color_codes(c,:) = [0.8500 0.3250 0.0980];
        c=c+1;
    elseif strcmp(FutureArray(i).scenario_description, 'R23_F0_EL0_FT100')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Recultivation, BECCS FT-diesel';
        color_codes(c,:) = [0 0.4470 0.7410];

        c=c+1;
    elseif strcmp(FutureArray(i).scenario_description, 'R11_F0_EL100_FT0')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Accelerated recultivation, BECCS bioelectricity';
        color_codes(c,:) = [0.9290 0.6940 0.1250];
        c=c+1;
    elseif strcmp(FutureArray(i).scenario_description, 'R11_F0_EL0_FT100')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Accelerated recultivation, BECCS FT-diesel';
        color_codes(c,:) = [0.3010 0.7450 0.9330];
        c=c+1;
    end



end

%% Plot final energy
src_data_fe_global = zeros(length(idx_of_interest),length(FutureArray(1).time));

time_this = FutureArray(1).time;

figure
line_width = 3;
c_fe = 1;
for idx = 1:length(idx_of_interest)
    if FutureArray(idx_of_interest(idx)).totals_vec_fe_ccs_tot(end) > 0
        plot(FutureArray(idx_of_interest(idx)).time(1:end), FutureArray(idx_of_interest(idx)).totals_vec_fe_ccs_tot(1:end),'LineWidth',line_width, 'Color', color_codes(idx,:));
        hold on
        labels_this{c_fe} = legend_labels{idx};
        c_fe = c_fe+1;
        %legend_labels{idx} = FutureArray(idx_of_interest(idx)).scenario_description;
        src_data_fe_global(idx,:) = FutureArray(idx_of_interest(idx)).totals_vec_fe_ccs_tot(1:end);
    end
end


xlim([2022 2050]);
xlabel('Year');
ylabel('EJ yr^{-1}');
legend(labels_this, 'location', 'northwest');

filename = ['Output/Future/global_final_energy.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)

save('Output/Future/src_data/src_data_fe_global.mat', 'src_data_fe_global', 'time_this');

%% Plot fluxes
src_data_fluxes_global = zeros(length(idx_of_interest),length(FutureArray(1).time));
figure
for idx = 1:length(idx_of_interest)
    plot(FutureArray(idx_of_interest(idx)).time(2:end), FutureArray(idx_of_interest(idx)).totals_vec_tco2eq_flux_net_nr_beccs(2:end),'LineWidth',line_width,'Color', color_codes(idx,:));
    hold on

    src_data_fluxes_global(idx,:) = FutureArray(idx_of_interest(idx)).totals_vec_tco2eq_flux_net_nr_beccs(1:end);

    %legend_labels{idx} = FutureArray(idx_of_interest(idx)).scenario_description;
end

xlim([2023 2050]);
xlabel('Year');
ylabel('GtCO_{2}eq yr^{-1}');
legend(legend_labels, 'location', 'southwest');

filename = ['Output/Future/global_net_carbon_fluxes.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)


save('Output/Future/src_data/src_data_fluxes_global.mat', 'src_data_fluxes_global', 'time_this');


%% Plot cumulatives
src_data_cumulatives_global = zeros(length(idx_of_interest),length(FutureArray(1).time));
figure
for idx = 1:length(idx_of_interest)
    tmp = FutureArray(idx_of_interest(idx)).cumulative_tco2_no_substitution_tot(1:end);
    tmp = tmp-tmp(1);
    plot(FutureArray(idx_of_interest(idx)).time(1:end), tmp(1:end),'LineWidth',line_width, 'Color', color_codes(idx,:));
    hold on
    %legend_labels{idx} = FutureArray(idx_of_interest(idx)).scenario_description;
    src_data_cumulatives_global(idx,:) = tmp;
end

xlim([2022 2050]);
xlabel('Year');
ylabel('GtCO_{2}eq');
legend(legend_labels, 'location', 'southwest');


filename = ['Output/Future/global_cumulative_co2.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)

save('Output/Future/src_data/src_data_cumulatives_global.mat', 'src_data_cumulatives_global', 'time_this');

%% Plot mitigation
figure

years_to_consider = [2030 2040 2050];
idx_years2consider = zeros(1,length(years_to_consider));

src_data = cell(1,5);
src_color_codes = cell(1,5);
src_data_signs = cell(1,5);
src_scenario_description = cell(1,5);
c=1;

for i = 1:length(years_to_consider)
    for j = 1:length(FutureArray(1).time)
        if years_to_consider(i) == FutureArray(1).time(j)
            idx_years2consider(i) = j;
        end
    end
end

time_since_t0 = FutureArray(1).time_since_t0(idx_years2consider);

x_tick_labels = cell(1,length(years_to_consider)+2);
x_tick_labels{1} = ' ';
x_tick_labels{end} = ' ';
%c = 1;
for t = 1:length(years_to_consider)
    time_since_t0_this = time_since_t0(t);
    idx_time_this = idx_years2consider(t);

    x_tick_labels{t+1} = num2str(years_to_consider(t));

    for idx = 1:length(idx_of_interest)
        color_this = color_codes(idx,:);
        
        % Harmonize to make t0 equal zero
        FA_this = FutureArray(idx_of_interest(idx));
        this_cum_no_sub = FA_this.cumulative_tco2_no_substitution_tot;
        
        this_cum_w_sub_ft_ccs = FA_this.cumulative_tco2_with_substitution_diesel_tot;

        this_cum_w_sub_elec_ccs_sub_fos_ccs = FA_this.cumulative_tco2_with_substitution_elec_fos_ccs_tot;
        this_cum_w_sub_elec_ccs_sub_nat_gas = FA_this.cumulative_tco2_with_substitution_elec_natural_gas_tot;
        this_cum_w_sub_elec_ccs_sub_coal = FA_this.cumulative_tco2_with_substitution_elec_coal_tot;

        this_cum_no_sub = this_cum_no_sub-this_cum_no_sub(1);
        this_cum_w_sub_ft_ccs = this_cum_w_sub_ft_ccs-this_cum_w_sub_ft_ccs(1);
        this_cum_w_sub_elec_ccs_sub_fos_ccs = this_cum_w_sub_elec_ccs_sub_fos_ccs-this_cum_w_sub_elec_ccs_sub_fos_ccs(1);
        this_cum_w_sub_elec_ccs_sub_nat_gas = this_cum_w_sub_elec_ccs_sub_nat_gas-this_cum_w_sub_elec_ccs_sub_nat_gas(1);
        this_cum_w_sub_elec_ccs_sub_coal = this_cum_w_sub_elec_ccs_sub_coal-this_cum_w_sub_elec_ccs_sub_coal(1);
        
        mit_this_avg_no_sub = this_cum_no_sub(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_ft_ccs = this_cum_w_sub_ft_ccs(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_elec_ccs_sub_fos_ccs = this_cum_w_sub_elec_ccs_sub_fos_ccs(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_elec_ccs_sub_nat_gas = this_cum_w_sub_elec_ccs_sub_nat_gas(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_elec_ccs_sub_coal = this_cum_w_sub_elec_ccs_sub_coal(idx_time_this)/time_since_t0_this;

        if strcmp(FA_this.scenario_description, 'Rinf_F100_EL0_FT0')
            scatter(t, mit_this_avg_no_sub, 'o', 'Color', color_this);
            src_data{c} = [t mit_this_avg_no_sub];
            src_data_signs{c} = 'o';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
            hold on

        elseif strcmp(FA_this.scenario_description, 'R23_F100_EL0_FT0')
            scatter(t, mit_this_avg_no_sub, 'o', 'Color', color_this);
            src_data{c} = [t mit_this_avg_no_sub];
            src_data_signs{c} = 'o';
            src_color_codes{c} = color_this;

            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
            hold on

        elseif FA_this.share_of_recultivation_to_beccs_ft == 1
            scatter(t, mit_this_avg_no_sub, 'o', 'Color', color_this);
            src_data{c} = [t mit_this_avg_no_sub];
            src_data_signs{c} = 'o';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
            hold on
            scatter(t, mit_this_avg_w_sub_ft_ccs, 'd', 'Color', color_this); 
            src_data{c} = [t mit_this_avg_w_sub_ft_ccs];
            src_data_signs{c} = 'd';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;

        elseif FA_this.share_of_recultivation_to_beccs_electricity == 1
            scatter(t, mit_this_avg_no_sub, 'o', 'Color', color_this);
            src_data{c} = [t mit_this_avg_no_sub];
            src_data_signs{c} = 'o';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
            hold on

            scatter(t, mit_this_avg_w_sub_elec_ccs_sub_fos_ccs, '+', 'Color', color_this); 
            src_data{c} = [t mit_this_avg_w_sub_elec_ccs_sub_fos_ccs];
            src_data_signs{c} = '+';
            src_color_codes{c} = color_this;
            c=c+1;
            scatter(t, mit_this_avg_w_sub_elec_ccs_sub_nat_gas, '*', 'Color', color_this); 
            src_data{c} = [t mit_this_avg_w_sub_elec_ccs_sub_nat_gas];
            src_data_signs{c} = '*';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
            scatter(t, mit_this_avg_w_sub_elec_ccs_sub_coal, 'x', 'Color', color_this); 
            src_data{c} = [t mit_this_avg_w_sub_elec_ccs_sub_coal];
            src_data_signs{c} = 'x';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
        end
    end
end

xlabel('Year');
ylabel('Average GtCO_{2}eq year^{-1}');
%xlim([0 length(years_to_consider)+1]);
%xticklabels(x_tick_labels);
box on
grid on

filename = ['Output/Future/global_mitigation.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)

fprintf('TEST \n')

save('src_data_test.mat', 'src_data', 'src_color_codes', 'src_data_signs', 'src_scenario_description');

end % FUNCTION









