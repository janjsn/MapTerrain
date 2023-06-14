filename = 'src_data_test.mat';
filename_nor = 'Output/Future/src_data/src_data_mitigation_Norway.mat';
filename_Trondelag = 'Output/Future/src_data/src_data_mitigation_Trøndelag.mat';

Src_data = open(filename);

figure

sz=500;
lw = 2;

for i = 1:length(Src_data.src_data)
data_this = Src_data.src_data{i};
x(i) = data_this(1);
y(i) = data_this(2);

s = scatter(x(i),y(i),sz, Src_data.src_data_signs{i}, 'LineWidth',lw);
hold on
s.MarkerEdgeColor = Src_data.src_color_codes{i};

end

xlim([0.5 3.5]);
ylim([-2.5 0.5]);
xticks([1 2 3]);
ylabel('GtCO2eq yr-1')
xlabel('Year');
xticklabels({'2030', '2040', '2050'});
box on


%% Norway
Src_data = open(filename_nor);

figure

sz=500;
lw = 2;

for i = 1:length(Src_data.src_data)
data_this = Src_data.src_data{i};
x(i) = data_this(1);
y(i) = 10^-6*data_this(2);

s = scatter(x(i),y(i),sz, Src_data.src_data_signs{i}, 'LineWidth',lw);
hold on
s.MarkerEdgeColor = Src_data.src_color_codes{i};

end

xlim([0.5 3.5]);
ylim([-0.4 0.1]);
xticks([1 2 3]);
yticks([-0.4:0.1:0.1])
ylabel('MtCO2eq yr-1')
xlabel('Year');
xticklabels({'2030', '2040', '2050'});
title('Norway')
box on


%% Trøndelag
Src_data = open(filename_Trondelag);

figure

sz=500;
lw = 2;

for i = 1:length(Src_data.src_data)
data_this = Src_data.src_data{i};
x(i) = data_this(1);
y(i) = 10^-6*data_this(2);

s = scatter(x(i),y(i),sz, Src_data.src_data_signs{i}, 'LineWidth',lw);
hold on
s.MarkerEdgeColor = Src_data.src_color_codes{i};

end

xlim([0.5 3.5]);
ylim([-0.16 0.04]);
xticks([1 2 3]);
yticks([-0.16:0.04:0.04])
ylabel('MtCO2eq yr-1')
xlabel('Year');
xticklabels({'2030', '2040', '2050'});
title('Trøndelag')
box on

