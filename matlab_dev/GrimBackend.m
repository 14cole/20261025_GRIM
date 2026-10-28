classdef GrimBackend < handle
    % Owns datasets and connects every action in GrimPlotApp to MATLAB code.
    % No Python bridge, subprocess, third-party package, or paid toolbox.
    properties (SetAccess=private)
        Datasets = {}
    end
    properties (Access=private)
        App
        NextID = 1
        Busy = false
    end
    methods
        function backend = GrimBackend(app)
            backend.App=app;
            app.OnAction=@(ui,action,state,payload) backend.handle(action,state,payload);
        end

        function addDatasets(backend,ds)
            % Public entry point for callers with in-memory GrimData structs.
            if isstruct(ds), ds=num2cell(ds); end
            for k=1:numel(ds), ds{k}=GrimData.validate(ds{k}); end
            rows=[];
            for k=1:numel(ds)
                d=GrimData.validate(ds{k}); d.ID=sprintf('dataset_%d',backend.NextID);
                backend.NextID=backend.NextID+1;
                if isempty(d.Name), d.Name=d.ID; end
                backend.Datasets{end+1}=d; rows(end+1)=numel(backend.Datasets); %#ok<AGROW>
            end
            backend.refresh(rows);
        end

        function handle(backend,action,state,payload)
            if backend.Busy, return; end
            backend.Busy=true; release=onCleanup(@() backend.unlock()); %#ok<NASGU>
            app=backend.App; selected=state.SelectedRows;
            switch action
                case 'load'
                    loaded={}; errors={};
                    for k=1:numel(payload.Files)
                        path=payload.Files{k}; app.setStatus(['Loading ' path]); drawnow limitrate nocallbacks
                        try, loaded{end+1}=GrimIO.load(path); %#ok<AGROW>
                        catch err, errors{end+1}=sprintf('%s\n%s',path,err.message); end %#ok<AGROW>
                    end
                    if ~isempty(loaded), backend.addDatasets(loaded); end
                    app.setStatus(sprintf('Loaded %d dataset(s).',numel(loaded)));
                    if ~isempty(errors), uialert(app.UIFigure,strjoin(errors,[newline newline]),'Some files could not be loaded'); end
                case 'datasets_selected'
                    backend.populateParameters(selected);
                case 'parameters_changed'
                    app.setStatus('Parameter selection updated.');
                case 'rename_dataset'
                    assert(~isempty(strtrim(char(string(payload.NewName)))),'Dataset name must not be empty.');
                    backend.Datasets{payload.Row}.Name=char(string(payload.NewName));
                    app.setStatus('Dataset renamed.');
                case 'delete_datasets'
                    backend.requireSelected(selected);
                    backend.Datasets(selected)=[];
                    if isempty(backend.Datasets), backend.refresh([]);
                    else, backend.refresh(min(selected(1),numel(backend.Datasets))); end
                case {'save','save_all','export_datasets'}
                    if strcmp(action,'save_all'), selected=1:numel(backend.Datasets); end
                    backend.requireSelected(selected);
                    backend.saveDatasets(selected,strcmp(action,'export_datasets'),state);
                case 'plot'
                    backend.requireSelected(selected);
                    GrimPlots.render(app,backend.Datasets(selected),state);
                case 'export_plot'
                    [file,folder]=uiputfile({'*.png','PNG image';'*.pdf','Vector PDF';'*.jpg','JPEG image';'*.tif','TIFF image'},'Export plot','plot.png');
                    if isequal(file,0), return; end
                    path=fullfile(folder,file);
                    exportgraphics(app.getAxes(),path,'Resolution',300);
                    app.setStatus(['Exported ' path]);
                case 'clear_plot'
                    app.setStatus('Plot cleared.');
                case 'plot_options_changed'
                    app.setStatus('Plot options updated. Click Plot to redraw, or enable Auto Plot.');
                otherwise
                    backend.requireSelected(selected);
                    results=backend.operation(action,backend.Datasets(selected),state);
                    if ~isempty(results)
                        backend.addDatasets(results);
                        app.setStatus(sprintf('%s: created %d dataset(s).',strrep(action,'_',' '),numel(results)));
                    end
            end
        end
    end

    methods (Access=private)
        function unlock(backend)
            backend.Busy=false;
        end

        function requireSelected(~,rows)
            assert(~isempty(rows),'Select one or more datasets first.');
        end

        function refresh(backend,selection)
            n=numel(backend.Datasets); rows=cell(n,2); ids=cell(n,1);
            for k=1:n
                d=backend.Datasets{k}; rows(k,:)={d.Name,d.Source}; ids{k}=d.ID;
            end
            backend.App.setDatasets(rows,ids);
            backend.App.selectDatasets(selection);
            backend.populateParameters(selection);
        end

        function populateParameters(backend,selection)
            app=backend.App;
            if isempty(selection), app.setParameters({},[],[],[]); return; end
            % First selected row defines parameter coordinates; every other
            % dataset is matched by coordinate value, never by array index.
            d=backend.Datasets{selection(1)};
            app.setParameters(d.Polarization,d.Frequency,d.Elevation,d.Azimuth);
            app.ParameterLists.Azimuth.Value=1:numel(d.Azimuth);
            % One cut initially; frequency/elevation sweeps can be selected.
            app.setStatus(sprintf('%d selected | %s | %s | deg, GHz', ...
                numel(selection),d.Name,GrimData.logUnit(d)));
        end

        function saveDatasets(backend,selection,isExport,state)
            app=backend.App; filter={'*.grim','GRIM native';'*.mat','MATLAB dataset'};
            ext='.grim';
            if isExport
                [choice,ok]=listdlg('PromptString','Export format','SelectionMode','single', ...
                    'ListString',{'CSV','TSV','Pioneer PIO (one cut)','PTM (one cut)'});
                if ~ok, return; end
                extensions={'.csv','.tsv','.pio','.ptm'}; ext=extensions{choice};
                filter={['*' ext],upper(ext(2:end))};
            end
            if numel(selection)==1
                d=backend.Datasets{selection};
                [file,folder]=uiputfile(filter,'Save dataset',[GrimBackend.fileName(d.Name) ext]);
                if isequal(file,0), return; end
                paths={fullfile(folder,file)};
            else
                folder=uigetdir(pwd,'Save selected datasets'); if isequal(folder,0), return; end
                paths=cell(size(selection));
                for k=1:numel(selection)
                    stem=GrimBackend.fileName(backend.Datasets{selection(k)}.Name);
                    path=fullfile(folder,[stem ext]); suffix=2;
                    while isfile(path)||any(strcmp(paths,path))
                        path=fullfile(folder,sprintf('%s_%d%s',stem,suffix,ext)); suffix=suffix+1;
                    end
                    paths{k}=path;
                end
            end
            for k=1:numel(selection)
                row=selection(k); d=backend.Datasets{row}; path=paths{k};
                [~,~,extension]=fileparts(path);
                if ismember(lower(extension),{'.pio','.ptm'})
                    assert(numel(state.Parameters.Elevation)==1 && numel(state.Parameters.Polarization)==1, ...
                        'Choose exactly one elevation and polarization for PIO/PTM export.');
                    % Export full azimuth/frequency cut at selected El/Pol.
                    values=struct('Azimuth',d.Azimuth,'Elevation',state.Parameters.Elevation, ...
                        'Frequency',d.Frequency,'Polarization',{state.Parameters.Polarization});
                    d=GrimData.select(d,values);
                end
                app.setStatus(['Saving ' path]); drawnow limitrate nocallbacks
                GrimIO.save(path,d);
                if ~isExport, backend.Datasets{row}.Source=path; end
            end
            if ~isExport, backend.refresh(selection); end
            app.setStatus(sprintf('Saved %d dataset(s).',numel(selection)));
        end

        function results = operation(~,action,ds,state)
            results={};
            switch action
                case {'coherent_add','coherent_subtract','coherent_divide','incoherent_add','incoherent_subtract','db_difference'}
                    results={GrimData.combine(ds,action)};
                case 'join', results={GrimData.join(ds)};
                case 'overlap', results=GrimData.overlap(ds);
                case 'align'
                    assert(numel(ds)>=2,'Select a reference dataset first, followed by datasets to align.');
                    names=GrimData.axes(); values=struct();
                    for k=1:4, values.(names{k})=ds{1}.(names{k}); end
                    for k=2:numel(ds)
                        results{end+1}=GrimData.derived(GrimData.select(ds{k},values),'align'); %#ok<AGROW>
                    end
                case 'crop'
                    for k=1:numel(ds)
                        results{end+1}=GrimData.derived(GrimData.select(ds{k},state.Parameters),'crop'); %#ok<AGROW>
                    end
                case 'regrid'
                    [dim,ok]=listdlg('PromptString','Axis to interpolate (complex samples where phase is known)', ...
                        'SelectionMode','single','ListString',{'Azimuth (deg)','Elevation (deg)','Frequency (GHz)'});
                    if ~ok, return; end
                    answer=inputdlg({'New coordinates: start:step:stop or comma-separated values'},'Regrid',[1 65],{'-180:1:180'});
                    if isempty(answer), return; end
                    values=GrimBackend.parseNumbers(answer{1}); names=GrimData.axes();
                    for k=1:numel(ds), results{end+1}=GrimData.regrid(ds{k},names{dim},values); end %#ok<AGROW>
                case {'statistics','percentile'}
                    percentile=90; interpolate=false;
                    if strcmp(action,'statistics')
                        methods={'mean','median','min','max','std'};
                        [choice,ok]=listdlg('PromptString','Statistic on linear power (repeated on original grid)', ...
                            'SelectionMode','single','ListString',methods); if ~ok, return; end
                        method=methods{choice};
                        [dims,ok]=listdlg('PromptString','Axes to reduce','SelectionMode','multiple', ...
                            'ListString',GrimData.axes(),'InitialValue',1); if ~ok, return; end
                    else
                        method='percentile'; dims=1;
                        answer=inputdlg({'Azimuth percentile (0 to 100)'},'Percentile',1,{'90'});
                        if isempty(answer), return; end
                        percentile=str2double(answer{1});
                        choice=questdlg('Percentile method','Percentile','Nearest rank','Linear interpolation','Cancel','Nearest rank');
                        if isempty(choice)||strcmp(choice,'Cancel'), return; end
                        interpolate=strcmp(choice,'Linear interpolation');
                    end
                    for k=1:numel(ds), results{end+1}=GrimData.reduce(ds{k},dims,method,percentile,interpolate); end %#ok<AGROW>
                case 'medianize'
                    answer=inputdlg({'Azimuth window width (deg)','Slide between output centers (deg)'},'Medianize',1,{'10','1'});
                    if isempty(answer), return; end
                    values=GrimBackend.scalars(answer);
                    for k=1:numel(ds), results{end+1}=GrimData.medianize(ds{k},values(1),values(2)); end %#ok<AGROW>
                case 'wrap'
                    choices={'Azimuth [-180,180)','Azimuth [0,360)','Phase [-180,180)','Phase [0,360)'};
                    [choice,ok]=listdlg('PromptString','Wrap coordinates or stored phase','SelectionMode','single','ListString',choices);
                    if ~ok, return; end
                    for k=1:numel(ds)
                        if choice<=2
                            origin=-180; if choice==2, origin=0; end
                            out=GrimData.transform(ds{k},'wrap',origin);
                        else
                            out=ds{k}; origin=-pi; if choice==4, origin=0; end
                            out.Phase=mod(out.Phase-origin,2*pi)+origin;
                            out=GrimData.derived(out,'wrap_phase');
                        end
                        results{end+1}=out; %#ok<AGROW>
                    end
                otherwise
                    switch action
                        case {'duplicate','flip_phase','swap_el_az'}, values=[];
                        case 'offset'
                            answer=inputdlg({'Level offset (dB)'},'Offset',1,{'0'});
                        case 'mirror'
                            answer=inputdlg({'Mirror about azimuth (deg)'},'Mirror',1,{'0'});
                        case 'shift'
                            answer=inputdlg({'Azimuth shift (deg)','Elevation shift (deg)','Phase shift (deg)'},'Shift',1,{'0','0','0'});
                        case 'round'
                            answer=inputdlg({'Azimuth decimal places','Elevation decimal places','Frequency decimal places'},'Round',1,{'3','3','6'});
                        case 'extrusion'
                            answer=inputdlg({'Uniform extrusion length (m); converts dBsm <-> dBke'},'Extrusion conversion',1,{'1'});
                        otherwise, error('Unknown action: %s',action);
                    end
                    if ~ismember(action,{'duplicate','flip_phase','swap_el_az'})
                        if isempty(answer), return; end
                        values=GrimBackend.scalars(answer);
                        if strcmp(action,'round'), assert(all(values==fix(values)),'Decimal places must be integers.'); end
                    end
                    for k=1:numel(ds), results{end+1}=GrimData.transform(ds{k},action,values); end %#ok<AGROW>
            end
        end
    end

    methods (Static)
        function values = scalars(answers)
            values=cellfun(@str2double,answers);
            assert(all(isfinite(values)),'Enter finite numeric values.');
        end

        function values = parseNumbers(text)
            % Parse coordinate expressions without eval/str2num.
            text=strtrim(regexprep(char(text),'[\[\]]',''));
            if contains(text,':')
                parts=strsplit(text,':'); assert(ismember(numel(parts),[2 3]),'Use start:step:stop.');
                nums=cellfun(@str2double,parts); assert(all(isfinite(nums)),'Invalid range.');
                step=1; if numel(nums)==3, step=nums(2); end
                assert(step~=0&&(nums(end)-nums(1))*step>=0,'Step does not lead toward the stop value.');
                assert(abs((nums(end)-nums(1))/step)<1e7,'Requested coordinate grid is too large.');
                values=nums(1):step:nums(end);
            else
                parts=regexp(text,'[,\s]+','split'); values=cellfun(@str2double,parts);
            end
            assert(~isempty(values)&&all(isfinite(values)),'Enter numeric coordinates.');
            values=unique(values);
        end

        function name = fileName(text)
            name=regexprep(char(text),'[<>:"/\\|?*\x00-\x1F]','_');
            name=regexprep(strtrim(name),'[. ]+$','');
            if isempty(name), name='dataset'; end
        end
    end
end
