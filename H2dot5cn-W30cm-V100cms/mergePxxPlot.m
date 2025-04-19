clc; clear; close all;
mergeData = load('merged_pxx.mat');

H = 0.0251;

% Measuring points along vertical centre line
xv = mergeData.f_merge; 
yv = repmat(double(mergeData.Y_merge), 1, size(xv, 2)); 
vv = xv .* mergeData.pxx_merge; % pre-PSD
smooth_vv = zeros(size(vv));

% Plot and save the merged PSD
outputFolder = 'merge-prePSD-figure';
if ~exist(outputFolder, 'dir')
    mkdir(outputFolder);
end

for ii = 1:length(mergeData.Y_merge)
    figure('Position', [100, 100 560 420]);
    % EMD smooth data
    smooth_window = {"gaussian", 50};
    % nmode = 4;
    % [imf,residual] = emd(vv(ii, :));
    % reconstrct = sum(imf(:, end-nmode:end), 2) + residual;
    % smooth_vv(ii, :) = smoothdata(reconstrct, smooth_window{:});
    % smooth_vv(ii, :) = reconstrct;
    
    smooth_vv(ii, :) = smoothdata(vv(ii, :), smooth_window{:});

    plot(xv(ii, :), vv(ii, :), 'LineWidth', 1, 'Color', 'b'); hold on;
    plot(xv(ii, :), smooth_vv(ii, :), 'LineWidth', 1.5, 'Color', 'r');
    set(gca, 'XScale', 'log'); % set(gca, 'YScale', 'log');
    set(gca, 'FontSize', 16);
    set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex');
    set(ylabel("$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
    set(title(sprintf("PSD at $z = %.3f$ m", mergeData.Y_merge(ii))), 'Interpreter', 'latex');
    saveas(gcf, fullfile(outputFolder, sprintf("merged_psd_z%.3f.png", mergeData.Y_merge(ii))));
end

%% 
% Create a single figure
figure('Position', [100, 100, 800, 600]);

% Define the indices for the five data series
indices = 5:5:length(mergeData.Y_merge)-2; % Select index at given step, ignore near wall data (max y index)

% Extract the corresponding z values for these indices
z_values = mergeData.Y_merge(indices); % Vertical positions for the selected indices

% Normalize the selected z values to [0, 1] for color mapping
z_min = min(z_values);
z_max = max(z_values);
z_normalized = (z_values - z_min) / (z_max - z_min); % Normalize to [0, 1]

% Define a colormap (e.g., 'jet') for the five selected series
cmap = sky(length(indices)); % Colormap with 5 colors (one for each selected series)

% Plot the selected data series with a color gradient
hold on;
for idx = 1:length(indices)
    ii = indices(idx); % Get the actual index
    % Get the color for this vertical position
    color_idx = round(z_normalized(idx) * (size(cmap, 1) - 1)) + 1; % Map to colormap index
    plot(xv(ii, :), smooth_vv(ii, :), 'LineWidth', 1.5, 'Color', cmap(color_idx, :));
end
hold off;

% Set the axes properties
set(gca, 'XScale', 'log');
% set(gca, 'YScale', 'log');
set(gca, 'FontSize', 16);

% Set labels and title with LaTeX interpreter
set(xlabel('$f$ (Hz)'), 'Interpreter', 'latex', 'FontSize', 16);
set(ylabel('$fS_{uu}(f) (\rm m^2/s^2)$'), 'Interpreter', 'latex', 'FontSize', 16);
set(title('Pre-multiplied PSD along z direction'), 'Interpreter', 'latex', 'FontSize', 16);

% Add a colorbar to show the mapping of colors to z values
colormap(cmap);
cbar = colorbar;
set(cbar, 'FontSize', 16);
set(get(cbar, 'Label'), 'String', '$z$ (m)', 'Interpreter', 'latex', 'FontSize', 16);
clim([z_min z_max]); % Set the colorbar limits to the range of selected z values

% Save the figure
base_filename = fullfile(outputFolder, 'merged_psd_selected');
% Save as MATLAB figure (.fig)
savefig(gcf, [base_filename '.fig']);

% Save as EMF (.emf)
print(gcf, [base_filename '.emf'], '-dmeta');

% Save as high-DPI JPEG using print (e.g., 500 DPI)
print(gcf, [base_filename '.jpg'], '-djpeg', '-r500');

%%
% Interpolate data
[grid_row, grid_col] = deal(max(400, size(yv, 1)*10), max(1200, round(size(xv, 2)/100)));
xq = logspace(...
    log10(min(xv(xv > 0))), log10(max(xv(:))), ...
    grid_col);
yq = linspace(...
    min(yv(:)), max(yv(:)), ...
    grid_row);
[xq, yq] = meshgrid(xq, yq);

vq = griddata(xv, yv, smooth_vv, xq, yq);

%% Plot contour
% limite y axis to water depth and ignore near wall region
figure('Position', [100, 100, 800, 600]);
region_index = yq >= z_values(end-3) & yq <= H;
xq_valid = xq .* region_index;
yq_valid = yq .* region_index;
vq_valid = vq .* region_index;
contourf(xq_valid, yq_valid / H, vq_valid, 10, 'LineStyle', '--');
set(gca, 'XScale', 'log'); set(gca, 'FontSize', 16); % set(gca, 'YScale', 'log');
set(xlabel("$f$ (Hz)"), 'Interpreter', 'latex');
set(ylabel("$z/H$"), 'Interpreter', 'latex');
colormap("sky");
col = colorbar();
set(ylabel(col,"$fS_{uu}(f) (\rm m^2/s^2)$"), 'Interpreter', 'latex');
axis([(min(xv(xv > 0))) (max(xv(:))) z_values(end-3)/H 1]); % limite y axis to water depth and ignore near wall region

% Save the figure
contour_filename = fullfile(outputFolder, 'PSD_contour');
% Save as MATLAB figure (.fig)
savefig(gcf, [contour_filename '.fig']);
% Save as EMF (.emf)
print(gcf, [contour_filename '.emf'], '-dmeta');
% Save as high-DPI JPEG using print (e.g., 500 DPI)
print(gcf, [contour_filename '.jpg'], '-djpeg', '-r500');

