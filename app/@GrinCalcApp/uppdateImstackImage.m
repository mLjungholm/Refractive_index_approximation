function uppdateImstackImage(app)
% Check if to show cropped
if app.ShowCropedIm.Value    

    imshow(app.lens.imStack(:,:,app.LayerSpinner.Value),'Parent',app.UIAxes)
    if app.ShowLines.Value
        % Check if there are any support lines
        if ~isempty(app.lens.supLines)
            app.lens.AppPlotLines(app.UIAxes,'support')
            app.lens.AppPlotLines(app.UIAxes,'support',app.supLineInd)
        end
        % Check if there are any sampling lines
%         if app.SamplingLinesCheckBox.Value
        if ~isnan(app.lens.lineNums)
            app.lens.AppPlotLines(app.UIAxes,'sampling')
            app.lens.AppPlotLines(app.UIAxes,'sampling',app.SamplingLineSpinner.Value)
        end
    end
    if app.RefractiveIndexOverlayCheckBox.Value && ~isempty(app.lens.n_overlay) 
        app.UIAxes.NextPlot = 'add';
        h = imshow(app.lens.n_overlay,'Parent',app.UIAxes);
        set(h, 'AlphaData', app.lens.n_overlay_map)
    end
    
else
    imshow(app.lens.imStackRaw(:,:,app.LayerSpinner.Value),'Parent',app.UIAxes)
    if app.ShowCropRec.Value && ~isempty(app.lens.cropRect)
        drawrectangle(app.UIAxes,'Position',app.lens.cropRect,'Color',[0.4660 0.6740 0.1880],'InteractionsAllowed',"none");
    end
end
