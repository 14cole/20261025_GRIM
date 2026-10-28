classdef GrimPlotApp < handle
    % GRIMPLOTAPP Simplified GRIM plotting GUI for MATLAB R2024b.
    % Base MATLAB only. The included GrimBackend connects all visible actions.
    %
    % RUN (put this file on the MATLAB path):
    %   app = GrimPlotApp;
    %   app = GrimPlotApp('gui');  % Optional empty GUI without a backend.
    %
    % REPLACE THE INCLUDED BACKEND (OPTIONAL):
    %   app.OnAction = @myBackend;
    % The handler signature is:
    %   function myBackend(app, action, state, payload)
    %       switch action
    %           case 'load'
    %               % payload.Files is a cell array of full file paths.
    %               % Read files, then update the GUI with setDatasets.
    %           case 'datasets_selected'
    %               % Read state.SelectedDatasetIDs, then setParameters.
    %           case 'plot'
    %               ax = app.getAxes();
    %               % plot(ax, x, y) or polarplot(ax, thetaRadians, rho)
    %           case 'crop'
    %               % Open your dialog / call your crop routine.
    %       end
    %   end
    %
    % POPULATE THE GUI (illustrative metadata, no data are loaded):
    %   app.setDatasets({'Run A','A.mat'; 'Run B','B.mat'}, {'a','b'});
    %   app.selectDatasets(1);  % Does not emit a callback.
    %   app.setParameters({'HH','VV'}, [8 9 10], [0 5], -180:10:180);
    %   app.UserData = yourBackendState;
    %
    % API:
    %   setDatasets(rows, ids) : N-by-2 cell/string array; unique IDs optional.
    %   selectDatasets(rows)  : Select rows in underlying table Data order.
    %   setParameters(pol, frequency, elevation, azimuth)
    %                         : Values in YOUR native units; no conversion.
    %   getState()            : Snapshot of IDs, rows, parameters and options.
    %   getAxes()             : Currently visible Cartesian or polar axes.
    %   setStatus(text)       : Set the footer message.
    %   requestPlot()         : Emit a 'plot' request using current controls.
    %
    % ACTIONS (all use the same OnAction function):
    %   load, save, save_all, export_datasets, delete_datasets,
    %   rename_dataset, datasets_selected, parameters_changed,
    %   plot, plot_options_changed, export_plot, clear_plot,
    %   plus the operation IDs listed in buildOperations below.
    % Rename payload: Row, ID, PreviousName, NewName. The GUI edits the
    % displayed name; your backend owns persistence and validation.
    % Parameter-change payload: Parameter. Values/indices are in state.
    % Other actions supply an empty struct unless documented above.
    % Multi-dataset operation order is TABLE ROW ORDER, not click order.
    % Updating GUI values through the public methods emits no callbacks.
    % Auto Plot requests a plot after user dataset/parameter/option changes.
    % Hold, grid and auto-scale are also applied to the visible axes.
    % Scale, phase, markers, and legend are implemented by GrimPlots.
    %
    % File-drop support is intentionally not included: Load uses MATLAB's
    % native multi-file dialog. No third-party drag/drop helper is required.
    % No ISAR, beta operations, provenance/history, or axis-units controls.

    properties
        OnAction = []       % @(app, action, state, payload)
        UserData = []       % Reserved for your backend.
    end

    properties (SetAccess = private)
        UIFigure
        DatasetTable
        ParameterLists = struct()
        Buttons = struct()
        Controls = struct()
        CartesianAxes
        PolarAxes
        Backend
    end

    properties (Access = private)
        MainGrid
        OperationsPanel
        SettingsPanel
        PlotGrid
        CartesianPanel
        PolarPanel
        StatusLabel
        DatasetIDs = cell(0,1)
        ParameterValues = struct('Polarization', {{}}, ...
            'Frequency', [], 'Elevation', [], 'Azimuth', [])
    end

    methods
        function app = GrimPlotApp(mode)
            if nargin<1, mode='full'; end
            app.buildUI();
            if ~strcmpi(mode,'gui')
                app.Backend=GrimBackend(app);
                app.setStatus('Ready. Load datasets to begin. Internal axes: degrees and GHz.');
            end
            app.UIFigure.Visible = 'on';
        end

        function delete(app)
            if ~isempty(app.UIFigure) && isgraphics(app.UIFigure)
                app.UIFigure.CloseRequestFcn = [];
                delete(app.UIFigure);
            end
        end

        function setDatasets(app, rows, ids)
            % Table contains display metadata only; keep arrays in UserData.
            if isempty(rows)
                rows = cell(0,2);
            end
            assert((iscell(rows) || isstring(rows)) && size(rows,2) == 2, ...
                'rows must be an N-by-2 cell or string array: Name, Source.');
            rows = cellstr(string(rows));
            n = size(rows,1);
            if nargin < 3
                ids = cellstr("dataset_" + string((1:n)'));
            end
            ids = cellstr(string(ids(:)));
            assert(numel(ids) == n && numel(unique(ids)) == n, ...
                'Supply one unique ID per dataset.');
            app.DatasetTable.Selection = [];
            app.DatasetTable.Data = rows;
            app.DatasetIDs = ids;
            app.setParameters({}, [], [], []);
            app.setStatus(sprintf('%d dataset(s)', n));
        end

        function selectDatasets(app, rows)
            if isempty(rows)
                app.DatasetTable.Selection = [];
                return
            end
            validateattributes(rows, {'numeric'}, ...
                {'vector','integer','positive','<=',size(app.DatasetTable.Data,1)});
            app.DatasetTable.Selection = rows(:)';
        end

        function setParameters(app, pol, frequency, elevation, azimuth)
            % Numeric values are retained exactly; labels are display only.
            values = {pol, frequency, elevation, azimuth};
            names = {'Polarization','Frequency','Elevation','Azimuth'};
            for k = 1:4
                v = values{k};
                if ischar(v)
                    v = {v};
                end
                v = reshape(v, 1, []);
                app.ParameterValues.(names{k}) = v;
                list = app.ParameterLists.(names{k});
                list.Value = {};
                list.ItemsData = [];
                list.Items = cellstr(string(v));
                list.ItemsData = 1:numel(v);
                if ~isempty(v)
                    list.Value = 1;
                end
            end
        end

        function state = getState(app)
            rows = sort(app.DatasetTable.Selection);
            state.SelectedRows = rows;
            state.SelectedDatasetIDs = app.DatasetIDs(rows);
            names = fieldnames(app.ParameterLists);
            state.ParameterIndices = struct();
            state.Parameters = struct();
            for k = 1:numel(names)
                name = names{k};
                indices = app.ParameterLists.(name).Value;
                if isempty(indices)
                    indices = [];
                end
                state.ParameterIndices.(name) = indices;
                values = app.ParameterValues.(name);
                state.Parameters.(name) = values(indices);
            end
            state.PlotType = app.Controls.PlotType.Value;
            names = {'AutoPlot','AutoScale','Hold','Legend','Grid', ...
                'Phase','Markers','Scale','XMin','XMax','YMin','YMax', ...
                'ColorMin','ColorMax','Colormap','Colorbar'};
            state.Options = struct();
            for k = 1:numel(names)
                state.Options.(names{k}) = app.Controls.(names{k}).Value;
            end
        end

        function ax = getAxes(app)
            if strcmp(app.Controls.PlotType.Value, 'azimuth_polar')
                ax = app.PolarAxes;
            else
                ax = app.CartesianAxes;
            end
        end

        function setStatus(app, message)
            app.StatusLabel.Text = char(string(message));
        end

        function requestPlot(app)
            app.showAxes();
            if ~app.validLimits()
                return
            end
            app.applyAxesOptions();
            if app.emit('plot', struct()) && isgraphics(app.UIFigure)
                % Plot functions may reset axes properties unless held.
                app.applyAxesOptions();
            end
        end
    end

    methods (Access = private)
        function buildUI(app)
            app.UIFigure = uifigure('Name','GRIM | Plotting', ...
                'Position',[80 80 1400 820], 'Visible','off', ...
                'Color',[0.95 0.96 0.97]);
            app.UIFigure.CloseRequestFcn = @(~,~) delete(app);
            outer = uigridlayout(app.UIFigure,[2 1]);
            outer.RowHeight = {'1x',24};
            outer.Padding = [10 10 10 10];
            app.MainGrid = uigridlayout(outer,[1 3]);
            app.MainGrid.Layout.Row = 1;
            app.MainGrid.ColumnWidth = {350,230,'1x'};
            app.MainGrid.Padding = [0 0 0 0];
            app.MainGrid.Scrollable = 'on';
            app.buildSidebar();
            app.buildOperations();
            app.buildPlotArea();
            app.StatusLabel = uilabel(outer, 'Text','Ready', ...
                'FontColor',[0.3 0.35 0.4]);
            app.StatusLabel.Layout.Row = 2;
            app.UIFigure.SizeChangedFcn = @(~,~) app.updateLayout();
            app.updateLayout();
        end

        function buildSidebar(app)
            sidebar = uigridlayout(app.MainGrid,[2 1]);
            sidebar.Layout.Column = 1;
            sidebar.RowHeight = {'3x','2x'};
            sidebar.Padding = [0 0 0 0];
            datasets = uipanel(sidebar,'Title','Datasets','FontWeight','bold');
            datasets.Layout.Row = 1;
            g = uigridlayout(datasets,[3 3]);
            g.RowHeight = {30,30,'1x'};
            specs = {'Load...','load'; 'Save','save'; 'Save All','save_all'; ...
                'Export...','export_datasets'; 'Delete','delete_datasets'};
            for k = 1:size(specs,1)
                b = app.actionButton(g,specs{k,1},specs{k,2});
                b.Layout.Row = ceil(k/3);
                b.Layout.Column = mod(k-1,3)+1;
            end
            app.Buttons.load.ButtonPushedFcn = @(~,~) app.chooseFiles();
            app.Buttons.delete_datasets.Layout.Column = [2 3];
            app.DatasetTable = uitable(g, 'Data',cell(0,2), ...
                'ColumnName',{'Name','Source / Output'}, ...
                'ColumnEditable',[true false], 'ColumnSortable',false, ...
                'ColumnWidth',{'1x','1x'}, 'RowName',{}, ...
                'SelectionType','row','Multiselect','on');
            app.DatasetTable.Layout.Row = 3;
            app.DatasetTable.Layout.Column = [1 3];
            app.DatasetTable.SelectionChangedFcn = @(~,~) app.selectionChanged();
            app.DatasetTable.CellEditCallback = @(~,e) app.nameEdited(e);

            p = uipanel(sidebar,'Title','Parameters','FontWeight','bold');
            p.Layout.Row = 2;
            g = uigridlayout(p,[3 4]);
            g.RowHeight = {24,'1x',22};
            names = {'Polarization','Frequency','Elevation','Azimuth'};
            for k = 1:4
                label = uilabel(g,'Text',names{k},'FontSize',11);
                label.Layout.Row = 1;
                label.Layout.Column = k;
                list = uilistbox(g,'Items',{},'Multiselect','on');
                list.Layout.Row = 2;
                list.Layout.Column = k;
                name = names{k};
                list.ValueChangedFcn = @(~,~) app.parametersChanged(name);
                app.ParameterLists.(name) = list;
            end
            hint = uilabel(g,'Text','Ctrl / Shift: select multiple values', ...
                'FontSize',11,'FontColor',[0.4 0.4 0.4]);
            hint.Layout.Row = 3;
            hint.Layout.Column = [1 4];
        end

        function buildOperations(app)
            app.OperationsPanel = uipanel(app.MainGrid, ...
                'Title','Dataset Operations','FontWeight','bold');
            app.OperationsPanel.Layout.Column = 2;
            g = uigridlayout(app.OperationsPanel,[16 2]);
            g.RowHeight = [repmat({30},1,15), {'1x'}];
            g.RowSpacing = 6;
            g.Scrollable = 'on';
            app.operationGroup(g,1,'Combine',{
                'Coherent +','coherent_add'; 'Coherent -','coherent_subtract';
                'Coherent /','coherent_divide'; 'Incoherent +','incoherent_add';
                'Incoherent -','incoherent_subtract'; 'Delta dB','db_difference';
                'Join / Merge','join'; 'Overlap','overlap'});
            app.operationGroup(g,6,'Transform',{
                'Crop / Slice','crop'; 'Statistics','statistics';
                'Percentile','percentile'; 'Align','align';
                'Regrid','regrid'; 'Mirror','mirror';
                'Wrap','wrap'; 'Shift','shift';
                'Round','round'; 'Offset','offset';
                'dBke / dBsm','extrusion'; 'Flip Phase','flip_phase';
                'Medianize','medianize'; 'Duplicate','duplicate'});
            b = app.actionButton(g,'Swap Elevation / Azimuth','swap_el_az');
            b.Layout.Row = 14;
            b.Layout.Column = [1 2];
        end

        function operationGroup(app, g, row, titleText, specs)
            label = uilabel(g,'Text',titleText,'FontWeight','bold');
            label.Layout.Row = row;
            label.Layout.Column = [1 2];
            for k = 1:size(specs,1)
                b = app.actionButton(g,specs{k,1},specs{k,2});
                b.Layout.Row = row + ceil(k/2);
                b.Layout.Column = mod(k-1,2)+1;
            end
        end

        function buildPlotArea(app)
            app.PlotGrid = uigridlayout(app.MainGrid,[5 1]);
            app.PlotGrid.Layout.Column = 3;
            app.PlotGrid.RowHeight = {32,32,30,0,'1x'};
            app.PlotGrid.Padding = [0 0 0 0];

            top = uigridlayout(app.PlotGrid,[1 4]);
            top.Layout.Row = 1;
            top.ColumnWidth = {155,'1x',110,110};
            top.Padding = [0 0 0 0];
            app.Controls.Operations = uibutton(top,'state', ...
                'Text','Dataset Operations','Value',true, ...
                'ValueChangedFcn',@(~,~) app.toggleOperations());
            label = uilabel(top,'Text','Plotting','FontWeight','bold','FontSize',16);
            label.Layout.Column = 2;
            b = app.actionButton(top,'Export Plot','export_plot');
            b.Layout.Column = 3;
            app.Controls.Settings = uibutton(top,'state','Text','Plot Settings', ...
                'ValueChangedFcn',@(~,~) app.toggleSettings());
            app.Controls.Settings.Layout.Column = 4;

            mode = uigridlayout(app.PlotGrid,[1 4]);
            mode.Layout.Row = 2;
            mode.ColumnWidth = {'1x',80,80,80};
            mode.Padding = [0 0 0 0];
            app.Controls.PlotType = uidropdown(mode, ...
                'Items',{'Azimuth (Rectangular)','Azimuth (Polar)','Frequency', ...
                'Elevation Sweep','Waterfall','RF Compare','Delta Map','Sector Statistics'}, ...
                'ItemsData',{'azimuth_rect','azimuth_polar','frequency', ...
                'elevation_sweep','waterfall','compare','delta_map','sector_stats'}, ...
                'ValueChangedFcn',@(~,~) app.plotTypeChanged());
            app.Buttons.plot = uibutton(mode,'Text','Plot', ...
                'BackgroundColor',[0.13 0.38 0.62],'FontColor',[1 1 1], ...
                'ButtonPushedFcn',@(~,~) app.requestPlot());
            app.Buttons.clear_plot = uibutton(mode,'Text','Clear', ...
                'ButtonPushedFcn',@(~,~) app.clearPlot());
            app.Buttons.fit = uibutton(mode,'Text','Fit Axes', ...
                'ButtonPushedFcn',@(~,~) app.fitAxes());

            options = uigridlayout(app.PlotGrid,[1 5]);
            options.Layout.Row = 3;
            options.Padding = [0 0 0 0];
            labels = {'Auto Plot','Auto Scale','Hold','Legend','Grid'};
            names = {'AutoPlot','AutoScale','Hold','Legend','Grid'};
            defaults = [false true false true true];
            for k = 1:numel(names)
                app.Controls.(names{k}) = uicheckbox(options,'Text',labels{k}, ...
                    'Value',defaults(k),'ValueChangedFcn',@(~,~) app.optionsChanged());
            end
            app.buildSettings();
            app.CartesianPanel = uipanel(app.PlotGrid,'BorderType','none');
            app.CartesianPanel.Layout.Row = 5;
            g = uigridlayout(app.CartesianPanel,[1 1]);
            g.Padding = [0 0 0 0];
            app.CartesianAxes = uiaxes(g);
            title(app.CartesianAxes,'');
            xlabel(app.CartesianAxes,'Azimuth');
            ylabel(app.CartesianAxes,'Response');
            grid(app.CartesianAxes,'on');
            app.PolarPanel = uipanel(app.PlotGrid,'BorderType','none','Visible','off');
            app.PolarPanel.Layout.Row = 5;
            app.PolarAxes = polaraxes(app.PolarPanel,'Units','normalized', ...
                'Position',[0.12 0.12 0.76 0.76]);
            app.PolarAxes.ThetaZeroLocation = 'top';
            app.PolarAxes.ThetaDir = 'clockwise';
        end

        function buildSettings(app)
            app.SettingsPanel = uipanel(app.PlotGrid,'Title','Plot Settings','Visible','off');
            app.SettingsPanel.Layout.Row = 4;
            g = uigridlayout(app.SettingsPanel,[4 4]);
            g.RowHeight = {26,26,26,26};
            g.ColumnWidth = {100,'1x',100,'1x'};
            specs = {'X minimum','XMin',-180; 'X maximum','XMax',180; ...
                'Y minimum','YMin',-80; 'Y maximum','YMax',0; ...
                'Color minimum','ColorMin',-80; 'Color maximum','ColorMax',0};
            for k = 1:size(specs,1)
                row = ceil(k/2);
                col = 1 + 2*mod(k-1,2);
                label = uilabel(g,'Text',specs{k,1});
                label.Layout.Row = row;
                label.Layout.Column = col;
                field = uieditfield(g,'numeric','Value',specs{k,3}, ...
                    'ValueChangedFcn',@(~,~) app.optionsChanged());
                field.Layout.Row = row;
                field.Layout.Column = col+1;
                app.Controls.(specs{k,2}) = field;
            end
            app.Controls.Scale = uidropdown(g,'Items',{'Dataset default','dB','Linear'}, ...
                'ValueChangedFcn',@(~,~) app.optionsChanged());
            app.Controls.Scale.Layout.Row = 4;
            app.Controls.Scale.Layout.Column = 1;
            app.Controls.Scale.Tooltip = 'Response scaling requested from your plot backend';
            app.Controls.Colormap = uidropdown(g,'Items',{'turbo','parula','gray','jet','hot'}, ...
                'ValueChangedFcn',@(~,~) app.optionsChanged());
            app.Controls.Colormap.Layout.Row = 4;
            app.Controls.Colormap.Layout.Column = 2;
            app.Controls.Colormap.Tooltip = 'Image colormap';
            extra = uigridlayout(g,[1 3]);
            extra.Layout.Row = 4;
            extra.Layout.Column = [3 4];
            extra.Padding = [0 0 0 0];
            extra.ColumnWidth = {62,75,80};
            for name = {'Phase','Markers','Colorbar'}
                key = name{1};
                app.Controls.(key) = uicheckbox(extra,'Text',key, ...
                    'Value',strcmp(key,'Colorbar'), ...
                    'ValueChangedFcn',@(~,~) app.optionsChanged());
            end
        end

        function b = actionButton(app, parent, label, action)
            b = uibutton(parent,'Text',label, ...
                'ButtonPushedFcn',@(~,~) app.emit(action,struct()));
            app.Buttons.(action) = b;
        end

        function chooseFiles(app)
            [files, folder] = uigetfile({ ...
                '*.grim;*.mat;*.csv;*.tsv;*.pio;*.cmplx_di;*.ptm;*.ss;*.out;*.cst_data;*.txt;*.dat;*.asc;*.ascii','Supported datasets'; ...
                '*.*','All files'},'Load datasets','MultiSelect','on');
            if isequal(files,0)
                return
            end
            if ischar(files)
                files = {files};
            end
            payload.Files = cellfun(@(f) fullfile(folder,f),files,'UniformOutput',false);
            app.emit('load',payload);
        end

        function selectionChanged(app)
            % Clear stale coordinates before your handler supplies new ones.
            app.setParameters({},[],[],[]);
            if app.emit('datasets_selected',struct())
                app.maybePlot();
            end
        end

        function nameEdited(app, event)
            row = event.Indices(1);
            payload = struct('Row',row,'ID',app.DatasetIDs{row}, ...
                'PreviousName',event.PreviousData,'NewName',event.NewData);
            if ~app.emit('rename_dataset',payload)
                app.DatasetTable.Data{row,1}=event.PreviousData;
            end
        end

        function parametersChanged(app, name)
            if app.emit('parameters_changed',struct('Parameter',name))
                app.maybePlot();
            end
        end

        function plotTypeChanged(app)
            app.showAxes();
            app.optionsChanged();
        end

        function optionsChanged(app)
            if app.emit('plot_options_changed',struct())
                if app.validLimits()
                    app.applyAxesOptions();
                    app.maybePlot();
                end
            end
        end

        function maybePlot(app)
            if isgraphics(app.UIFigure) && app.Controls.AutoPlot.Value
                app.requestPlot();
            end
        end

        function showAxes(app)
            polar = strcmp(app.Controls.PlotType.Value,'azimuth_polar');
            app.CartesianPanel.Visible = app.onOff(~polar);
            app.PolarPanel.Visible = app.onOff(polar);
        end

        function ok = validLimits(app)
            c = app.Controls;
            ok = c.AutoScale.Value || (c.XMin.Value < c.XMax.Value && ...
                c.YMin.Value < c.YMax.Value && c.ColorMin.Value < c.ColorMax.Value);
            if ~ok
                app.setStatus('Each minimum must be below its maximum, or enable Auto Scale.');
            end
        end

        function applyAxesOptions(app)
            ax = app.getAxes();
            c = app.Controls;
            grid(ax,app.onOff(c.Grid.Value));
            hold(ax,app.onOff(c.Hold.Value));
            if isequal(ax,app.PolarAxes)
                if c.AutoScale.Value
                    rlim(ax,'auto');
                else
                    rlim(ax,[c.YMin.Value c.YMax.Value]);
                end
            else
                if c.AutoScale.Value
                    xlim(ax,'auto'); ylim(ax,'auto'); clim(ax,'auto');
                else
                    xlim(ax,[c.XMin.Value c.XMax.Value]);
                    ylim(ax,[c.YMin.Value c.YMax.Value]);
                    clim(ax,[c.ColorMin.Value c.ColorMax.Value]);
                end
                colormap(ax,c.Colormap.Value);
            end
        end

        function clearPlot(app)
            cla(app.CartesianAxes);
            cla(app.PolarAxes);
            legend(app.CartesianAxes,'off');
            legend(app.PolarAxes,'off');
            colorbar(app.CartesianAxes,'off');
            app.emit('clear_plot',struct());
        end

        function fitAxes(app)
            app.Controls.AutoScale.Value = true;
            app.applyAxesOptions();
            app.emit('plot_options_changed',struct());
        end

        function toggleOperations(app)
            visible = app.Controls.Operations.Value;
            app.OperationsPanel.Visible = app.onOff(visible);
            app.updateLayout();
        end

        function updateLayout(app)
            % Preserve usable control sizes; small windows scroll instead
            % of crushing the controls. Large windows grow the plot area.
            position = app.UIFigure.Position;
            operationsWidth = 230 * double(app.Controls.Operations.Value);
            plotWidth = max(540, position(3) - 20 - 20 - 350 - operationsWidth);
            app.MainGrid.ColumnWidth = {350,operationsWidth,plotWidth};
            app.MainGrid.RowHeight = {max(500,position(4)-54)};
        end

        function toggleSettings(app)
            visible = app.Controls.Settings.Value;
            app.SettingsPanel.Visible = app.onOff(visible);
            if visible
                app.PlotGrid.RowHeight = {32,32,30,170,'1x'};
            else
                app.PlotGrid.RowHeight = {32,32,30,0,'1x'};
            end
        end

        function ok = emit(app, action, payload)
            ok = true;
            if isempty(app.OnAction)
                app.setStatus(sprintf('%s: no backend callback connected.',action));
                return
            end
            try
                app.OnAction(app,action,app.getState(),payload);
            catch err
                ok = false;
                if isgraphics(app.UIFigure)
                    app.setStatus(sprintf('%s failed: %s',action,err.message));
                    uialert(app.UIFigure,err.message,'Backend callback error');
                end
                warning('GrimPlotApp:BackendError','%s',getReport(err,'extended','hyperlinks','off'));
            end
        end
    end

    methods (Static, Access = private)
        function value = onOff(condition)
            if condition
                value = 'on';
            else
                value = 'off';
            end
        end
    end
end
