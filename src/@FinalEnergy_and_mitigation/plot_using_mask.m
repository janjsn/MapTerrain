function plot_using_mask(obj, fraction_of_cell_is_region_1deg, region_name)


idx_of_interest = 0;
c=1;
% for i = 1:length(obj.FutureArray)
%     if strcmp(obj.FutureArray(i).scenario_description, 'R23_F100_EL0_FT0')
%         idx_of_interest(c) = i;
%
%         idx_baseline = i;
%
%         legend_labels{c} = 'Baseline';
%         c=c+1;
%     elseif strcmp(obj.FutureArray(i).scenario_description, 'Rinf_F100_EL0_FT0')
%         idx_of_interest(c) = i;
%         legend_labels{c} = 'Continued regrowth';
%         c=c+1;
%         idx_nat_reg = i;
%
%     elseif strcmp(obj.FutureArray(i).scenario_description, 'R23_F0_EL100_FT0')
%         idx_of_interest(c) = i;
%         legend_labels{c} = 'Recultivation, BECCS bioelectricity';
%         c=c+1;
%     elseif strcmp(obj.FutureArray(i).scenario_description, 'R23_F0_EL0_FT100')
%         idx_of_interest(c) = i;
%         legend_labels{c} = 'Recultivation, BECCS FT-diesel';
%         c=c+1;
%     elseif strcmp(obj.FutureArray(i).scenario_description, 'R11_F0_EL100_FT0')
%         idx_of_interest(c) = i;
%         legend_labels{c} = 'Accelerated recultivation, BECCS bioelectricity';
%         c=c+1;
%     elseif strcmp(obj.FutureArray(i).scenario_description, 'R11_F0_EL0_FT100')
%         idx_of_interest(c) = i;
%         legend_labels{c} = 'Accelerated recultivation, BECCS FT-diesel';
%         c=c+1;
%     end
%
%
%
% end

color_codes = zeros(6,3);

for i = 1:length(obj.FutureArray)
    if strcmp(obj.FutureArray(i).scenario_description, 'R23_F100_EL0_FT0')
        idx_of_interest(c) = i;

        idx_baseline = i;

        legend_labels{c} = 'Baseline';
        color_codes(c,:) = [0.4940 0.1840 0.5560]; % Purple

        c=c+1;
    elseif strcmp(obj.FutureArray(i).scenario_description, 'Rinf_F100_EL0_FT0')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Continued regrowth';

        color_codes(c,:) = [34 139 34]/255	; % GREEN
        idx_nat_reg = i;
        c=c+1;

    elseif strcmp(obj.FutureArray(i).scenario_description, 'R23_F0_EL100_FT0')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Recultivation, BECCS bioelectricity';
        color_codes(c,:) = [0.8500 0.3250 0.0980];
        c=c+1;
    elseif strcmp(obj.FutureArray(i).scenario_description, 'R23_F0_EL0_FT100')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Recultivation, BECCS FT-diesel';
        color_codes(c,:) = [0 0.4470 0.7410]; % Dark blue
        c=c+1;
    elseif strcmp(obj.FutureArray(i).scenario_description, 'R11_F0_EL100_FT0')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Accelerated recultivation, BECCS bioelectricity';
        color_codes(c,:) = [0.9290 0.6940 0.1250];
        c=c+1;
    elseif strcmp(obj.FutureArray(i).scenario_description, 'R11_F0_EL0_FT100')
        idx_of_interest(c) = i;
        legend_labels{c} = 'Accelerated recultivation, BECCS FT-diesel';
        color_codes(c,:) = [0.3010 0.7450 0.9330]; % Light blue
        c=c+1;
    end

end


FA_interest = obj.FutureArray(idx_of_interest);

n_scens = length(FA_interest);

time = FA_interest(1).time;
time_since_t0 = FA_interest(1).time_since_t0;


land_to_bioenergy_tot = zeros(n_scens,length(time));
land_regrowing_tot = zeros(n_scens,length(time));

fe_ccs_tot = zeros(n_scens,length(time));
co2_flux = zeros(n_scens,length(time));

cumulative_co2 = zeros(n_scens,length(time));

for i = 1:length(FA_interest)
    for t = 1:length(time)
        land_be_this = FA_interest(i).land_converted_to_bioenergy(:,:,t);
        land_to_bioenergy_tot(i,t) = sum(sum(land_be_this.*fraction_of_cell_is_region_1deg));

        fe_ccs_tot_this = FA_interest(i).fe_ccs_tot(:,:,t);
        fe_ccs_tot(i,t) = sum(sum(fe_ccs_tot_this.*fraction_of_cell_is_region_1deg));

        co2_flux_this = FA_interest(i).tco2eq_flux_beccs_tot(:,:,t) + FA_interest(i).tco2eq_flux_natural_regrowth_net(:,:,t);
        co2_flux(i,t) = sum(sum(co2_flux_this.*fraction_of_cell_is_region_1deg));

        cumulative_co2_no_substitution_this = FA_interest(i).cumulative_tco2_no_substitution(:,:,t);
        cumulative_co2(i,t) = sum(sum(cumulative_co2_no_substitution_this.*fraction_of_cell_is_region_1deg));
        %cumulative_co2(i,t) = cumulative_co2(i,t)-cumulative_co2(i,1);


    end
end


%% Plot final energy
src_data_fe = zeros(length(FA_interest),length(time));

figure
c_fe = 1;
line_width = 3;
for i = 1:length(FA_interest)
    if fe_ccs_tot(i,end) > 0
        plot(time(1:end), fe_ccs_tot(i,:),'LineWidth',line_width, 'Color', color_codes(i,:));
        hold on
        labels_this{c_fe} = legend_labels{i};
        c_fe = c_fe+1;

        src_data_fe(i,:) = fe_ccs_tot(i,:);
        %legend_labels{idx} = FutureArray(idx_of_interest(idx)).scenario_description;
    end
end


xlim([2022 2050]);
xlabel('Year');
ylabel('EJ yr^{-1}');
legend(labels_this, 'location', 'northwest');

filename = ['Output/Future/' region_name '_final_energy.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)


filename_src = ['Output/Future/' region_name '_fe_src_data.mat'];
save(filename_src, 'src_data_fe', 'time')

%% Plot fluxes
src_data_fluxes = zeros(length(FA_interest),length(time));

figure

for i =1:length(FA_interest)
    plot(time(2:end), 10^-6*co2_flux(i,2:end),'LineWidth',line_width, 'Color', color_codes(i,:));
    hold on
    
    src_data_fluxes(i,:) = 10^-6*co2_flux(i,:);
end

xlim([2023 2050]);
xlabel('Year');
ylabel('MtCO_{2}eq yr^{-1}');
legend(legend_labels, 'location', 'southwest');

filename = ['Output/Future/' region_name '_co2_fluxes.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)

filename_src = ['Output/Future/' region_name '_fluxes_co2eq_src_data.mat'];
save(filename_src, 'src_data_fluxes', 'time')

%% Plot cumulatives
src_data_cumulative = zeros(length(FA_interest),length(time));

figure

for i =1:length(FA_interest)
    tmp = cumulative_co2(i,1:end);

    tmp = tmp-tmp(1);

    plot(time(1:end), 10^-6*tmp(1:end),'LineWidth',line_width, 'Color', color_codes(i,:));
    hold on
    src_data_cumulative(i,:) = tmp;
end

xlim([2022 2050]);
xlabel('Year');
ylabel('MtCO_{2}eq');
legend(legend_labels, 'location', 'southwest');

filename = ['Output/Future/' region_name '_cumulative_co2.pdf'];
if exist(filename, 'file')
    delete(filename);
end

print('-vector','-dpdf', '-r1000', filename)

filename_src = ['Output/Future/' region_name '_cumulative_co2_src_data.mat'];
save(filename_src, 'src_data_cumulative', 'time')


%% Produce source data for mitigation scatter plot


% REF Scarlat et al. (2022), Applied energy
emf_Norway_elec_prod = 7.8;
emf_EU27_elec_prod = 86.1;
% Emission factors, alternative techs
emf_natural_gas_elec = [136 146];
emf_diesel_fuel = 93.9;

years_to_consider = [2030 2040 2050];
idx_years2consider = zeros(1,length(years_to_consider));

src_data = cell(1,5);
src_color_codes = cell(1,5);
src_data_signs = cell(1,5);
src_scenario_description = cell(1,5);
c=1;

for i = 1:length(years_to_consider)
    for j = 1:length(obj.FutureArray(1).time)
        if years_to_consider(i) == obj.FutureArray(1).time(j)
            idx_years2consider(i) = j;
        end
    end
end

time_since_t0 = obj.FutureArray(1).time_since_t0(idx_years2consider);

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

        FA_this = obj.FutureArray(idx_of_interest(idx));

        % GET CUMULATIVES
        for th = 1:length(FA_this.time)
            % UNIT is tCO2eq yr-1
            this_cum_gwp_no_sub_global = FA_this.cumulative_tco2_no_substitution(:,:,th);
            this_cum_gwp_no_sub_mask_tot(th) = sum(sum(sum(this_cum_gwp_no_sub_global.*fraction_of_cell_is_region_1deg)));

            %diesel substitution
            this_cum_w_sub_ft_ccs_global = this_cum_gwp_no_sub_global-(FA_this.cumulative_fe_ft_ccs(:,:,th)*(10^-3*emf_diesel_fuel));
            this_cum_w_sub_ft_ccs_mask_tot(th) = sum(sum(sum(this_cum_w_sub_ft_ccs_global.*fraction_of_cell_is_region_1deg)));

            % Norwegian grid mix electricity substitution
            this_cum_w_sub_NOR_elec_global = this_cum_gwp_no_sub_global-(FA_this.cumulative_fe_bioelectricity_ccs(:,:,th)*(10^-3*emf_Norway_elec_prod));
            this_cum_w_sub_NOR_elec_mask_tot(th) = sum(sum(sum(this_cum_w_sub_NOR_elec_global.*fraction_of_cell_is_region_1deg)));

            % EU grid mix electricity substitution
            this_cum_w_sub_EU_elec_global = this_cum_gwp_no_sub_global-(FA_this.cumulative_fe_bioelectricity_ccs(:,:,th)*(10^-3*emf_EU27_elec_prod));
            this_cum_w_sub_EU_elec_mask_tot(th) = sum(sum(sum(this_cum_w_sub_EU_elec_global.*fraction_of_cell_is_region_1deg)));

            % Natural gas electricity substitution
            this_cum_w_sub_NG_elec_global = this_cum_gwp_no_sub_global-(FA_this.cumulative_fe_bioelectricity_ccs(:,:,th)*(10^-3*mean(emf_natural_gas_elec)));
            this_cum_w_sub_NG_elec_mask_tot(th) = sum(sum(sum(this_cum_w_sub_NG_elec_global.*fraction_of_cell_is_region_1deg)));

        end

        % Align to set initial time step to zero for visualization (no nat reg.).
        this_cum_gwp_no_sub_mask_tot = this_cum_gwp_no_sub_mask_tot-this_cum_gwp_no_sub_mask_tot(1);
        this_cum_w_sub_ft_ccs_mask_tot = this_cum_w_sub_ft_ccs_mask_tot-this_cum_w_sub_ft_ccs_mask_tot(1);
        this_cum_w_sub_NOR_elec_mask_tot = this_cum_w_sub_NOR_elec_mask_tot-this_cum_w_sub_NOR_elec_mask_tot(1);
        this_cum_w_sub_EU_elec_mask_tot = this_cum_w_sub_EU_elec_mask_tot-this_cum_w_sub_EU_elec_mask_tot(1);
        this_cum_w_sub_NG_elec_mask_tot = this_cum_w_sub_NG_elec_mask_tot-this_cum_w_sub_NG_elec_mask_tot(1);

        % Get average mitigation
        mit_this_avg_no_sub = this_cum_gwp_no_sub_mask_tot(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_ft_ccs = this_cum_w_sub_ft_ccs_mask_tot(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_NOR_elec = this_cum_w_sub_NOR_elec_mask_tot(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_EU_elec = this_cum_w_sub_EU_elec_mask_tot(idx_time_this)/time_since_t0_this;
        mit_this_avg_w_sub_elec_ccs_sub_nat_gas = this_cum_w_sub_NG_elec_mask_tot(idx_time_this)/time_since_t0_this;

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
            
            % Norwegian electricity
            scatter(t, mit_this_avg_w_sub_NOR_elec, '^', 'Color', color_this);
            src_data{c} = [t mit_this_avg_w_sub_NOR_elec];
            src_data_signs{c} = '^';
            src_color_codes{c} = color_this;
            c=c+1;
            % Natural gas
            scatter(t, mit_this_avg_w_sub_elec_ccs_sub_nat_gas, '*', 'Color', color_this);
            src_data{c} = [t mit_this_avg_w_sub_elec_ccs_sub_nat_gas];
            src_data_signs{c} = '*';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
            % EU27 electricity
            scatter(t, mit_this_avg_w_sub_EU_elec, 'v', 'Color', color_this);
            src_data{c} = [t mit_this_avg_w_sub_EU_elec];
            src_data_signs{c} = 'v';
            src_color_codes{c} = color_this;
            src_scenario_description{c} = FA_this.scenario_description;
            c=c+1;
        end
    end
end

%% NOTE: Plot script written to make use of the source data. See individual file.

save(['Output/Future/src_data/src_data_mitigation_' region_name '.mat'], 'src_data', 'src_data_signs', 'src_color_codes', 'src_scenario_description');

end % FUNCTION

