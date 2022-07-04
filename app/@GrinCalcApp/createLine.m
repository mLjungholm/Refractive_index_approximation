function createLine(app,style)

switch style
    case 'sampling'
        app.addTwoConsole('Draw a sampling line. Press "enter" when done or "Q" to cansel')
        roi = drawline(app.UIAxes);
        app.listenForKey = true;
        uiwait(app.UIFigure)
        if isequal(app.keyIn,'return')
            ind = app.lens.AppCreateSamplingLine(roi);
            if ~isnan(app.lens.pixelSize)
                app.lens.mLines{ind}.pixelSize = app.lens.pixelSize;
            end
            tdata = {ind,nan,false,true};
            app.SamplingLineTable.Data = sortrows([app.SamplingLineTable.Data;tdata],1);
            if app.lens.lineNums > 1
                app.SamplingLineSpinner.Enable = true;
%                 app.SamplingLineSpinner.Limits = [1, app.lens.lineNums]; 
            elseif app.lens.lineNums == 1
                app.SamplingLineSpinner.Enable = false;
%                 app.SamplingLineSpinner.Limits = [1, 2];
            end
            app.SamplingLineSpinner.Value = ind; 
            app.currentLineInd = ind;
            app.ShowLines.Value = true;
        else
            delete(roi)
        end
        app.uppdateImstackImage()
        app.keyIn = '';
        
    case 'support'      
        app.addTwoConsole('Draw a support line. Press "enter" when done or "Q" to cansel')
        roi = drawline(app.UIAxes);
        app.listenForKey = true;
        uiwait(app.UIFigure)
        if isequal(app.keyIn,'return')
            ind = app.lens.AppCreateSupportLine(roi);
%             lineName = strcat('Support Line','',num2str(ind));
%             uitreenode(app.SupportLinesNode,'Text',lineName,"NodeData",ind);
%             expand(app.LineTree)
%             app.lineInd = ind;
%             app.lineType = 'support';
%             app.SamplingLineSpinner.Enable = false;
            tdata = {ind,''};
            app.SupportLineTable.Data = sortrows([app.SupportLineTable.Data; tdata],1);
            app.ShowLines.Value = true;
        else
            delete(roi)
        end
        app.uppdateImstackImage()
        app.keyIn = '';
end
end