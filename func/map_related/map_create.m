interactive_box_mapper(abs(x2_rot));

function interactive_box_mapper(volume)
    % INPUT: volume - a 3D matrix (e.g. size 120x120x120)
    % OUTPUT: saved 'mappedVolume' and 'boxMask' when user clicks Save
    
    sz = size(volume);
    cx = round(sz(2)/2); cy = round(sz(1)/2); cz = round(sz(3)/2);

    % Initial box is half size of each dimension centered
    box = [floor(sz(2)/4), floor(sz(1)/4), floor(sz(3)/4), ...
           floor(sz(2)/2), floor(sz(1)/2), floor(sz(3)/2)];

    f = figure('Name','Interactive Volume Box Mapper','NumberTitle','off','Position',[100 100 1600 800]);

    % Axes for original views
    ax1 = subplot(2,3,1); im_xy = imshow(volume(:,:,cz),[]); title('XY View');
    ax2 = subplot(2,3,2); im_xz = imshow(squeeze(volume(cx,:,:))',[]); title('XZ View');
    ax3 = subplot(2,3,3); im_yz = imshow(squeeze(volume(:,cy,:))',[]); title('YZ View');

 

    % Rectangles on first row
    r1 = drawrectangle(ax1, 'Position',[box(1), box(2), box(4), box(5)]);
    r2 = drawrectangle(ax2, 'Position',[box(1), box(3), box(4), box(6)]);
    r3 = drawrectangle(ax3, 'Position',[box(2), box(3), box(5), box(6)]);

    % Axes for mapped volume views
    ax4 = subplot(2,3,4); im_m_xy = imshow(volume(:,:,cz),[]); title('Mapped XY');
    ax5 = subplot(2,3,5); im_m_yz = imshow(squeeze(volume(cx,:,:))',[]); title('Mapped YZ');
    ax6 = subplot(2,3,6); im_m_xz = imshow(squeeze(volume(:,cy,:))',[]); title('Mapped XZ');

    % Button
    btn = uicontrol('Style','pushbutton','String','Save','Position',[20 20 100 40],...
                    'Callback',@saveMask);

    % Update everything initially
    update();

    % Add listeners
    addlistener(r1, 'ROIMoved', @(~,~) roiMoved(1));
    addlistener(r2, 'ROIMoved', @(~,~) roiMoved(2));
    addlistener(r3, 'ROIMoved', @(~,~) roiMoved(3));

    function roiMoved(viewId)
        % Read positions
        pos1 = round(r1.Position);
        pos2 = round(r2.Position);
        pos3 = round(r3.Position);

        % Update box from ROIs
        box(1) = pos1(1);
        box(2) = pos1(2);
        box(3) = pos2(2);
        box(4) = pos1(3);
        box(5) = pos1(4);
        box(6) = pos2(4);

        % Sync ROIs
        r1.Position = [box(1), box(2), box(4), box(5)];
        r2.Position = [box(1), box(3), box(4), box(6)];
        r3.Position = [box(2), box(3), box(5), box(6)];

        update();
    end

    function update()
        % Clamp
        x_rng = max(1,box(2)) : min(sz(1), box(2)+box(5)-1);
        y_rng = max(1,box(1)) : min(sz(2), box(1)+box(4)-1);
        z_rng = max(1,box(3)) : min(sz(3), box(3)+box(6)-1);

        boxMask = false(sz);
        
        boxMask(x_rng, y_rng, z_rng) = true;

        mappedVolume = volume;
        mappedVolume(~boxMask) = 0;

        % --- Center indices within box (for mapped views)
        cx_m = round(mean(x_rng));
        cy_m = round(mean(y_rng));
        cz_m = round(mean(z_rng));

        % Show mapped views
        im_m_xy.CData = mappedVolume(:,:,cz_m);
        im_m_yz.CData = squeeze(mappedVolume(cx_m,:,:))';
        im_m_xz.CData = squeeze(mappedVolume(:,cy_m,:))';
    end

    function saveMask(~,~)
        % Final ranges
        x_rng = max(1,box(2)) : min(sz(1), box(2)+box(5)-1);
        y_rng = max(1,box(1)) : min(sz(2), box(1)+box(4)-1);
        z_rng = max(1,box(3)) : min(sz(3), box(3)+box(6)-1);
        range ={x_rng, y_rng, z_rng};
        boxMask = false(sz);
        boxMask(x_rng, y_rng, z_rng) = true;
        mappedVolume = volume;
        mappedVolume(~boxMask) = 0;
        croppedVolume = mappedVolume(x_rng, y_rng,z_rng);
        assignin('base','boxMask',boxMask);
        assignin('base','mappedVolume',mappedVolume);
        assignin('base','croppedVolume',croppedVolume);
        assignin('base','range',range);
        disp('Saved to workspace: boxMask, mappedVolume, croppedVolume, range');
        close(f);
    end
end