%%%%%%%%%
% ARRAY %
%%%%%%%%%

all_points = g.points();

PIG02 = dsearchn(all_points, [217.7, -115.8]);

fileID = fopen('illapel_amundsen_ice_q15_binary.bin', 'r');
Size = n_water;

A = fread(fileID, [2*n_water iter], 'single');

PIG02_data = A(PIG02, :);

t_axis = k_max:k_max:time_in_sec;
t_axis_hours = t_axis / 3600;

figure()

plot(t_axis_hours, PIG02_data*100, '-red', 'LineWidth', 2)
% ylim([-5, 5])

set(gca, 'FontSize', 20)
xlabel('$t$ (hours)', 'Interpreter', 'latex', 'visible', 'on')
ylabel('Vertical displacement, $w$ (cm)', 'Interpreter', 'latex', 'visible', 'on')