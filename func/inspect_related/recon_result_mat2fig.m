% April 3, 2025
% Yiwei Jia
% Standardize the visualization of mat file
clear;clc;
nImg=1;

%%
% Manually select two results
fig_name = 't1Libre_woBin';
fig_path = '/home/debi/yiwei/recon_results/250423/t1_libre_240_woBin/';
if isfolder(fig_path)==0
    mkdir(fig_path)
end

xpath{1} = ['/home/debi/yiwei/recon_results/250423/Sub001/T1_LIBRE_woBinning/output/mask_t1w_libre/' ...
    'x_nIter15_delta_1.000.mat']; 
x{1} = load(xpath{1}, 'x');
tt{1} = 'libre t1';

xpath{2} = ['/home/debi/yiwei/recon_results/250417/Sub003/T1_LIBRE_woBinning/output/mask_libre_t2w/' ...
    'x_nIter15_delta_1.000.mat']; 
x{2} = load(xpath{2}, 'x');
tt{2} = 'libre t2';


xpath{3} = ['/home/debi/yiwei/recon_results/250417/Sub004/T1_LIBRE_woBinning/output/mask_gre_original/' ...
    'x_nIter15_delta_1.000.mat']; 
x{3} = load(xpath{3}, 'x');
tt{3} = 'gre orig traj';

% xpath{3} = ['/home/debi/yiwei/recon_results/250402/' ...
%     'Sub003/T1_LIBRE_woBinning/output/mask_libre_rfsp_qua50_adc_ph/x0.mat'];
% x{3} = load(xpath{3}, 'x0');
% tt{3} = 'libre-rf50-adc-pi/2';
% xpath{4} = ['/home/debi/yiwei/recon_results/250402/' ...
%     'Sub004/T1_LIBRE_woBinning/output/mask_rect_rfsp_qua50_adc_ph/x0.mat'];
% x{4} = load(xpath{4}, 'x0');
% tt{4} = 'rect-rf50';
%%
% fig_name = 'coronal_sub9-sub14_0404_bern';
% tt = {'libre sp+50adc+ph', 'libre sp-50adc-ph', 'libre sp+50adc-ph', 'libre sp+50adc0', ...
%     'libre wosp', 'rf wosp'};
%% set the array list

for idx = 100:140
    fig=figure;set(gcf, 'Color', 'w');  
    
    % Determine subplot layout (square-like arrangement)
    % rows = ceil(sqrt(nImg));
    % cols = ceil(nImg / rows);
    rows = nImg;
    cols = 1;
    tiledlayout(rows, cols, 'Padding', 'none', 'TileSpacing', 'none');
    for i = 1:nImg
        if isa(x{1}, 'single')
            x_i = x{i};
            x_i = permute(x_i, [2, 1, 3]);  % size = [30, 20, 10]
            sl = abs(x_i(:,:,idx));
        else
            fields = fieldnames(x{i});
            x_i = x{i}.(fields{1});
            
            if isa(x_i, 'cell')
                x_i = x_i{1};
            end
            % x_i = permute(x_i, [2, 1, 3]);
            x_i = flip(x_i, 1);
            sl = mat2gray(abs(x_i(145:230,60:180,idx)));
            % sl = mat2gray(abs(x_i(:,:,idx)));
            if i==1
                sl = equal_func(sl, 0.01, 0.9);
            elseif i==2 
                sl = equal_func(sl, 0.01, 0.9);
            else
                sl = equal_func(sl, 0.01, 0.9);
            end
        end
        % subplot(rows, cols, i);
        ax = nexttile;
        imshow(sl, []);  % Display image
        axis off;  % Hide axes
        axis tight;
        text(ax, -0.05, 0.5, [tt{i},': ', num2str(idx)], ...
            'Units', 'normalized', ...
            'HorizontalAlignment', 'right', ...
            'VerticalAlignment', 'middle', ...
            'FontWeight', 'bold', ...
            'FontSize', 12);
        % title([tt{i},': ', num2str(idx)], 'FontSize', 12);  % Add title
    end
    disp(strcat(fig_path,fig_name,'_', num2str(idx),'.png'));
    exportgraphics(fig, strcat(fig_path,fig_name,'_', num2str(idx),'.png'), 'Resolution', 400);
    
end
