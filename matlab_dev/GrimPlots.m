classdef GrimPlots
    % Plot preparation/rendering ported from GRIM_Backend/plotting/modes.
    % Simplified controls: individual sweeps; azimuth/frequency image plane;
    % sector mean over selected azimuths. No ISAR or beta processing.
    methods (Static)
        function render(app,ds,state)
            assert(~isempty(ds),'Select a dataset before plotting.');
            for k=1:numel(ds), ds{k}=GrimData.select(ds{k},state.Parameters); end
            for k=2:numel(ds), GrimData.compatible(ds{1},ds{k}); end
            ax=app.getAxes(); options=state.Options; mode=state.PlotType;
            if ~options.Hold
                cla(ax,'reset');
                if isequal(ax,app.PolarAxes), ax.ThetaZeroLocation='top'; ax.ThetaDir='clockwise'; end
                if isequal(ax,app.CartesianAxes), colorbar(ax,'off'); end
            end
            hold(ax,'on'); marker='none'; if options.Markers, marker='o'; end
            switch mode
                case {'azimuth_rect','azimuth_polar','frequency','elevation_sweep'}
                    dim=1; if strcmp(mode,'frequency'), dim=3;
                    elseif strcmp(mode,'elevation_sweep'), dim=2; end
                    GrimPlots.lines(ax,ds,options,dim,strcmp(mode,'azimuth_polar'),marker);
                case 'waterfall'
                    assert(numel(ds)==1,'Waterfall displays one dataset at a time.');
                    d=ds{1}; GrimPlots.singlePlane(d);
                    values=GrimPlots.values(d,options);
                    c=reshape(values,numel(d.Azimuth),numel(d.Frequency)).';
                    GrimPlots.image(ax,d.Azimuth,d.Frequency,c);
                    title(ax,d.Name,'Interpreter','none');
                    GrimPlots.imageLabels(ax,options,GrimPlots.valueLabel(d,options));
                case 'delta_map'
                    assert(numel(ds)==2,'Delta Map requires exactly two datasets.');
                    assert(~options.Phase,'Delta Map compares magnitude levels; turn off Phase.');
                    d=ds{1}; GrimPlots.singlePlane(d);
                    c=GrimData.db(ds{1})-GrimData.db(ds{2}); c(~isfinite(c))=NaN;
                    c=reshape(c,numel(d.Azimuth),numel(d.Frequency)).';
                    GrimPlots.image(ax,d.Azimuth,d.Frequency,c);
                    title(ax,[ds{1}.Name ' minus ' ds{2}.Name],'Interpreter','none');
                    GrimPlots.imageLabels(ax,options,'Difference (dB)');
                case 'sector_stats'
                    assert(~options.Phase,'Sector statistics use linear power; turn off Phase.');
                    for k=1:numel(ds)
                        d=ds{k};
                        % Count a closed-turn direction only once in statistics.
                        if numel(d.Azimuth)>1 && abs(d.Azimuth(end)-d.Azimuth(1)-360)<1e-6
                            d=GrimData.slice(d,{1:numel(d.Azimuth)-1,1:numel(d.Elevation),1:numel(d.Frequency),1:numel(d.Polarization)});
                        end
                        reduced=GrimData.reduce(d,1,'mean');
                        GrimPlots.lines(ax,{reduced},options,1,false,marker);
                    end
                    title(ax,'Mean linear power over selected azimuth samples');
                case 'compare'
                    assert(numel(ds)==2,'RF Compare requires exactly two datasets.');
                    shape=GrimData.shape(ds{1}); dim=find(shape(1:3)>1,1);
                    assert(~isempty(dim)&&nnz(shape(1:3)>1)==1&&shape(4)==1, ...
                        'RF Compare: select one sweep axis, one value on other axes, and one polarization.');
                    GrimPlots.lines(ax,ds,options,dim,false,marker);
                    if options.Phase
                        result=GrimPlots.phaseAgreement(ds{1}.Phase(:),ds{2}.Phase(:));
                        title(ax,sprintf('Phase agreement %.1f / 100 | circular RMS %.3g deg',result.Score,result.RMS));
                    else
                        result=GrimPlots.agreement(ds{1}.Power(:),ds{2}.Power(:));
                        title(ax,sprintf('RF agreement %.1f / 100 | mean absolute difference %.3g dB',result.Score,result.MAE));
                    end
                otherwise, error('Unknown plot mode: %s',mode);
            end
            if ~ismember(mode,{'waterfall','delta_map'})
                if options.Legend, legend(ax,'show','Interpreter','none','Location','best');
                else, legend(ax,'off'); end
            end
            grid(ax,GrimPlots.onOff(options.Grid)); hold(ax,GrimPlots.onOff(options.Hold));
            app.setStatus(sprintf('%s | %d dataset(s)',strrep(mode,'_',' '),numel(ds)));
        end

        function y = values(d,options)
            if options.Phase, y=rad2deg(d.Phase);
            elseif strcmp(options.Scale,'Linear'), y=d.Power;
            else, y=GrimData.db(d); end
        end

        function label = valueLabel(d,options)
            if options.Phase, label='Phase (deg)';
            elseif strcmp(options.Scale,'Linear')
                switch d.Quantity
                    case 'sigma_3d', label='RCS (m^2)';
                    case 'sigma_2d', label='Scattering width (m)';
                    otherwise, label='Power ratio';
                end
            else, label=['Response (' GrimData.logUnit(d) ')']; end
        end

        function lines(ax,ds,options,dim,isPolar,marker)
            names=GrimData.axes(); count=0;
            for k=1:numel(ds), dims=GrimData.shape(ds{k}); count=count+prod(dims)/dims(dim); end
            assert(count<=400,'More than 400 curves requested. Select fewer parameter combinations.');
            for k=1:numel(ds)
                d=ds{k}; x=d.(names{dim}); dims=GrimData.shape(d);
                order=[dim,setdiff(1:4,dim,'stable')]; other=dims(order(2:end));
                y=reshape(permute(GrimPlots.values(d,options),order),numel(x),[]);
                for j=1:size(y,2)
                    [a,b,c]=ind2sub(other,j); indices=[a b c]; parts={d.Name};
                    for m=1:3
                        axis=order(m+1); v=d.(names{axis}); index=indices(m);
                        if axis==4, label=v{index}; else, label=sprintf('%.6g',v(index)); end
                        parts{end+1}=[names{axis} '=' label]; %#ok<AGROW>
                    end
                    label=strjoin(parts,' | ');
                    if isPolar
                        polarplot(ax,deg2rad(x),y(:,j),'LineWidth',1.2,'Marker',marker,'DisplayName',label);
                    else
                        plot(ax,x,y(:,j),'LineWidth',1.2,'Marker',marker,'DisplayName',label);
                    end
                end
            end
            if ~isPolar
                units={'deg','deg','GHz'}; xlabel(ax,[names{dim} ' (' units{dim} ')']);
                ylabel(ax,GrimPlots.valueLabel(ds{1},options));
            end
        end

        function singlePlane(d)
            assert(numel(d.Elevation)==1&&numel(d.Polarization)==1, ...
                'For an image, select one elevation and one polarization; select multiple azimuths/frequencies.');
        end

        function image(ax,x,y,c)
            % Edges retain physical spacing for nonuniform coordinate grids.
            xe=GrimPlots.edges(x); ye=GrimPlots.edges(y); [X,Y]=meshgrid(xe,ye);
            padded=c([1:end end],[1:end end]);
            surface(ax,X,Y,zeros(size(X)),padded,'EdgeColor','none','FaceColor','flat', ...
                'AlphaData',double(isfinite(padded)),'FaceAlpha','flat','AlphaDataMapping','none');
            view(ax,2); ax.YDir='normal'; xlabel(ax,'Azimuth (deg)'); ylabel(ax,'Frequency (GHz)');
        end

        function e = edges(x)
            x=x(:)';
            if numel(x)==1, e=[x-0.5,x+0.5];
            else, middle=(x(1:end-1)+x(2:end))/2; e=[x(1)-(middle(1)-x(1)),middle,x(end)+(x(end)-middle(end))]; end
        end

        function imageLabels(ax,options,label)
            colormap(ax,options.Colormap);
            if options.Colorbar, bar=colorbar(ax); bar.Label.String=label;
            else, colorbar(ax,'off'); end
            legend(ax,'off');
        end

        function result = agreement(a,b)
            valid=isfinite(a)&isfinite(b)&a>0&b>0; a=a(valid); b=b(valid);
            assert(numel(a)>=2,'RF agreement requires at least two positive finite matched samples.');
            globalScore=GrimPlots.componentScore(a,b); groups=GrimPlots.sectors(numel(a)); scores=zeros(1,numel(groups));
            for k=1:numel(groups), scores(k)=GrimPlots.componentScore(a(groups{k}),b(groups{k})); end
            result.Score=.75*globalScore+.25*GrimData.percentile(scores(:),25,true);
            result.MAE=mean(abs(10*log10(a)-10*log10(b)));
        end

        function score = componentScore(a,b)
            scale=max([abs(a(:));abs(b(:))]); if scale==0, score=100; return; end
            aa=a/scale; bb=b/scale;
            nrmse=sqrt(mean((aa-bb).^2)/mean((aa.^2+bb.^2)/2));
            score=100*(.45*max(GrimPlots.concordance(a,b),0)+ ...
                .30*max(GrimPlots.concordance(10*log10(a),10*log10(b)),0)+ ...
                .25*(1-min(nrmse/sqrt(2),1)));
        end

        function value = concordance(a,b)
            scale=max([abs(a(:));abs(b(:))]); if scale==0, value=1; return; end
            a=a/scale; b=b/scale; ma=mean(a); mb=mean(b); da=a-ma; db=b-mb;
            den=mean(da.^2)+mean(db.^2)+(ma-mb)^2;
            if den==0, value=double(isequal(a,b)); else, value=2*mean(da.*db)/den; end
        end

        function groups = sectors(n)
            count=min(6,max(1,floor(n/3))); widths=repmat(floor(n/count),1,count);
            widths(1:mod(n,count))=widths(1:mod(n,count))+1;
            offsets=[0 cumsum(widths)]; groups=cell(1,count);
            for k=1:count, groups{k}=(offsets(k)+1):offsets(k+1); end
        end

        function result = phaseAgreement(a,b)
            valid=isfinite(a)&isfinite(b); residual=angle(exp(1i*(a(valid)-b(valid))));
            assert(numel(residual)>=2,'Phase agreement requires two finite matched samples.');
            groups=GrimPlots.sectors(numel(residual)); scores=zeros(1,numel(groups));
            for k=1:numel(groups), scores(k)=50*(1+mean(cos(residual(groups{k})))); end
            result.Score=.75*50*(1+mean(cos(residual)))+.25*GrimData.percentile(scores(:),25,true);
            result.RMS=sqrt(mean(rad2deg(residual).^2));
        end

        function value = onOff(state)
            value='off'; if state, value='on'; end
        end
    end
end
