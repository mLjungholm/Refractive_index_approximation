function overwriteValue(app, cVar)

switch cVar
    
    case 'minLayer'
        if  app.signalFromApp
            app.signalFromApp = false;
            return
        end
        if rem(app.EdgeMinLayerEditField.Value,1) ~= 0
            app.signalFromApp = true;
            app.EdgeMinLayerEditField.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMin;
            return
        end
            answer = uiconfirm(app.UIFigure,'Overwrite min layer value?','Overwrite min layer','Icon','warning');
        if isequal(answer,'Cancel')
            app.signalFromApp = true;
            app.EdgeMinLayerEditField.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMin;
            return
        elseif isequal(answer,'OK')
            app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMin = app.EdgeMinLayerEditField.Value;
        end
        
    case 'maxLayer'
        if  app.signalFromApp
            app.signalFromApp = false;
            return
        end
        if rem(app.EdgeMaxLayerEditField.Value,1) ~= 0
            app.signalFromApp = true;
            app.EdgeMaxLayerEditField.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMax;
            return
        end
            answer = uiconfirm(app.UIFigure,'Overwrite max layer value?','Overwrite max layer','Icon','warning');
        if isequal(answer,'Cancel')
            app.signalFromApp = true;
            app.EdgeMaxLayerEditField.Value = app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMax;
            return
        elseif isequal(answer,'OK')
            app.lens.mLines{app.SamplingLineSpinner.Value}.leftPhaseMax = app.EdgeMaxLayerEditField.Value;
        end
end
        app.uppdateLinePlot()
end