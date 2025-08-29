%%
addpath(genpath('/Users/cag/Documents/forclone/mfuderer-colorResources-05c224c'));
%%
clear x;
nImg = 4;
tt{1} = 'IDEA LIBRE';
tt{2} = 'IDEA Rect';
tt{3} = 'LIBRE wo grad sp';
tt{4} = 'Rect wo grad sp';
x{1} =flip_3axis(x_sub6);
x{2} = flip_3axis(x_sub7);
x{3} = x_sub4;
x{4} = x_sub5;

fig_path = '/Users/cag/Documents/Dataset/recon_results/250502/sub4567/';
fig_name = 'sub4567';
%%
nImg = 2;
clear x;
tt{1} = 'grad sp';
tt{2} = 'wo grad sp';

x{1} = x_sub3;
x{2} = x_sub4;
fig_path = '/Users/cag/Documents/Dataset/recon_results/250502/sub34/';
if isfolder(fig_path)==0
    mkdir(fig_path)
end
fig_name = 'sub34';
%%
save_fig = 1;
for idx = 110:140
    fig=figure('Position', [200, 300, 800, 800], 'Resize', 'on');set(gcf, 'Color', 'w'); 
    % Determine subplot layout (square-like arrangement)
    rows = ceil(sqrt(nImg));
    cols = ceil(nImg / rows);
    
    t= tiledlayout(rows, cols, 'Padding', 'none', 'TileSpacing', 'none'); %compact or none
    
    for i = 1:nImg
       
        if isa(x{1}, 'single')
            x_i = x{i};
            disp('not a cell')
        else
            fields = fieldnames(x{i});
            x_i = x{i}.(fields{1});
            disp('loading from a cell')
        end

            x_i = permute(x_i, [2,3,1]);
            sl = mat2gray(abs(x_i(:,:,idx)));
            if i==1
                sl = equal_func(sl, 0, 0.75);
            elseif i==4
                sl = equal_func(sl, 0, 0.9);
            end
            
 
        ax = nexttile;
        loLev = 0.05;upLev = 2;[imClip, rgb_vec] = relaxationColorMap('T1', sl , loLev, upLev);
        imshow(imClip, 'DisplayRange', [loLev, upLev], 'InitialMagnification', 'fit'); 
        colormap(ax, rgb_vec);
        if mod(i,cols)==0
             colorbar;
        end
        axis off;  % Hide axes
        axis tight;
        title([tt{i}, ' sl', num2str(idx)], 'FontSize', 12);  % Add title
    end

    if save_fig
        exportgraphics(fig, strcat(fig_path,fig_name,'_', num2str(idx),'.png'), 'Resolution', 300);
    end

end


function x_flipped=flip_3axis(x)
x_flipped = flip(flip(flip(x, 1), 2), 3);
end