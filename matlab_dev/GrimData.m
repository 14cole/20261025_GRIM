classdef GrimData
    % Numerical core translated from GRIM_Backend/datasets/{grid,axes,
    % arithmetic,combine,transforms}. All arrays use [azimuth,elevation,
    % frequency,polarization] order. Internal units: degrees, GHz, radians.
    % Power is physical sigma (m^2), scattering width (m), or power ratio.
    % NaN power = missing sample; NaN phase = unknown phase, never zero phase.
    methods (Static)
        function d = create(az,el,f,pol,power,phase,quantity)
            if nargin < 7 || isempty(phase), phase = nan(size(power)); end
            if nargin < 8, quantity = 'sigma_3d'; end
            d = struct('Name','','Source','','ID','','Azimuth',double(az(:)'), ...
                'Elevation',double(el(:)'),'Frequency',double(f(:)'), ...
                'Polarization',{cellstr(upper(strtrim(string(reshape(pol,1,[])))))}, ...
                'Power',double(power),'Phase',double(phase),'Quantity',char(quantity), ...
                'Metadata',struct(),'History','');
            d = GrimData.validate(d);
        end

        function d = validate(d)
            names = GrimData.axes();
            assert(all(isfield(d,[names,{'Name','Source','ID','Power','Phase','Quantity','Metadata','History'}])), ...
                'Incomplete MATLAB dataset structure.');
            for k=1:3, d.(names{k})=double(reshape(d.(names{k}),1,[])); end
            d.Polarization=cellstr(upper(strtrim(string(reshape(d.Polarization,1,[])))));
            dims = GrimData.shape(d);
            assert(all(dims>0),'Dataset axes must not be empty.');
            assert(numel(d.Power)==prod(dims) && numel(d.Phase)==prod(dims), ...
                'Sample shape does not match the four axes.');
            d.Power = reshape(double(d.Power),dims);
            d.Phase = reshape(double(d.Phase),dims);
            assert(isreal(d.Power) && isreal(d.Phase),'Power and phase must be real arrays.');
            assert(~any(isinf(d.Power(:)) | d.Power(:)<0),'Power must be nonnegative finite or NaN.');
            assert(~any(isinf(d.Phase(:))),'Phase must be finite or NaN.');
            assert(ismember(d.Quantity,{'sigma_3d','sigma_2d','power_ratio'}),'Unknown physical quantity.');
            for k=1:4
                v = d.(names{k});
                if k<4, assert(all(isfinite(v)),'Coordinates must be finite.'); end
                assert(numel(unique(v))==numel(v),'Duplicate coordinates are not allowed.');
                if k<4
                    [d.(names{k}),order] = sort(v);
                    ix = {':',':',':',':'}; ix{k}=order;
                    d.Power=d.Power(ix{:}); d.Phase=d.Phase(ix{:});
                end
            end
            assert(all(d.Frequency>0),'Frequencies must be positive.');
            d.Phase(~isfinite(d.Power))=NaN;
        end

        function names = axes()
            names = {'Azimuth','Elevation','Frequency','Polarization'};
        end

        function n = shape(d)
            n = [numel(d.Azimuth),numel(d.Elevation),numel(d.Frequency),numel(d.Polarization)];
        end

        function unit = logUnit(d)
            switch d.Quantity
                case 'sigma_3d', unit='dBsm';
                case 'sigma_2d', unit='dBke';
                otherwise, unit='dB';
            end
        end

        function scale = logFactor(d)
            % GRIM dBke is 10*log10(k0 * physical scattering width).
            if strcmp(d.Quantity,'sigma_2d')
                scale=reshape(2*pi*d.Frequency*1e9/299792458,[1 1 numel(d.Frequency) 1]);
            else
                scale=1;
            end
        end

        function y = db(d)
            y=10*log10(d.Power.*GrimData.logFactor(d));
        end

        function z = field(d)
            assert(~any(isfinite(d.Power(:)) & d.Power(:)>0 & isnan(d.Phase(:))), ...
                'This operation requires known phase for every nonzero sample.');
            z=sqrt(d.Power).*exp(1i*d.Phase);
            z(d.Power==0)=0;
        end

        function d = fromField(d,z)
            d.Power=abs(z).^2; d.Phase=angle(z);
            d.Phase(~isfinite(d.Power))=NaN;
        end

        function ix = match(source, requested, tolerance)
            if nargin<3, tolerance=1e-6; end
            if iscell(source) || isstring(source)
                [found,ix]=ismember(string(requested),string(source));
                assert(all(found),'A selected polarization is not present in every dataset.');
                return
            end
            ix=zeros(1,numel(requested));
            for k=1:numel(requested)
                candidates=find(abs(source-requested(k))<=tolerance);
                assert(numel(candidates)==1,'A coordinate is missing or ambiguous within tolerance.');
                ix(k)=candidates;
            end
            assert(numel(unique(ix))==numel(ix),'Multiple selected coordinates match the same sample.');
        end

        function d = select(d,values)
            names=GrimData.axes(); ix=cell(1,4);
            for k=1:4
                v=values.(names{k});
                assert(~isempty(v),'Select at least one value on every parameter axis.');
                ix{k}=GrimData.match(d.(names{k}),v);
            end
            d=GrimData.slice(d,ix);
        end

        function d = slice(d,ix)
            names=GrimData.axes();
            d.Power=d.Power(ix{:}); d.Phase=d.Phase(ix{:});
            for k=1:4, d.(names{k})=d.(names{k})(ix{k}); end
        end

        function compatible(a,b)
            assert(strcmp(a.Quantity,b.Quantity),'Datasets have different physical quantities.');
            names=GrimData.axes();
            for k=1:4
                assert(numel(a.(names{k}))==numel(b.(names{k})),'Dataset axes differ; use Align or Regrid.');
                ix=GrimData.match(b.(names{k}),a.(names{k}));
                assert(isequal(ix(:)',1:numel(ix)),'Dataset axis order differs; use Align.');
            end
            ca=GrimData.meta(a,'angular_coordinate_system','conic');
            cb=GrimData.meta(b,'angular_coordinate_system','conic');
            assert(strcmp(ca,cb),'Datasets use different angular coordinate systems.');
            assert(strcmp(GrimData.meta(a,'elevation_coordinate_convention','signed_elevation'), ...
                GrimData.meta(b,'elevation_coordinate_convention','signed_elevation')), ...
                'Datasets use different elevation conventions (for example native SENTRi theta).');
        end

        function v = meta(d,key,fallback)
            v=fallback;
            if isfield(d.Metadata,key) && ~isempty(d.Metadata.(key))
                v=char(string(d.Metadata.(key)));
            end
        end

        function out = combine(ds,op)
            assert(numel(ds)>=2,'Select at least two datasets.');
            out=ds{1};
            for k=2:numel(ds), GrimData.compatible(out,ds{k}); end
            if ismember(op,{'coherent_divide','db_difference'})
                assert(numel(ds)==2,'This operation requires exactly two datasets.');
            end
            if startsWith(op,'coherent_')
                for k=2:numel(ds)
                    left=GrimData.meta(ds{1},'time_convention',''); right=GrimData.meta(ds{k},'time_convention','');
                    assert(isempty(left)||isempty(right)||strcmp(left,right),'Coherent time conventions differ.');
                end
                z=GrimData.field(out);
                for k=2:numel(ds)
                    next=GrimData.field(ds{k});
                    switch op
                        case 'coherent_add', z=z+next;
                        case 'coherent_subtract', z=z-next;
                        case 'coherent_divide'
                            valid=isfinite(next) & abs(next)>0;
                            result=complex(nan(size(z)),nan(size(z)));
                            result(valid)=z(valid)./next(valid); z=result;
                    end
                end
                z(~isfinite(z))=complex(NaN,NaN);
                if strcmp(op,'coherent_divide'), assert(any(isfinite(z(:))),'No usable nonzero denominator samples.'); end
                out=GrimData.fromField(out,z);
                if strcmp(op,'coherent_divide'), out.Quantity='power_ratio'; end
            else
                p=out.Power;
                for k=2:numel(ds)
                    q=ds{k}.Power;
                    switch op
                        case 'incoherent_add', p=p+q;
                        case 'incoherent_subtract'
                            difference=p-q;
                            tol=8*eps.*max(abs(p),abs(q));
                            assert(~any(difference(:)<-tol(:)), ...
                                'Incoherent subtraction would produce negative physical power.');
                            difference(difference<0)=0; p=difference;
                        case 'db_difference'
                            valid=isfinite(p) & isfinite(q) & q>0;
                            result=nan(size(p)); result(valid)=p(valid)./q(valid); p=result;
                        otherwise, error('Unknown combination operation: %s',op);
                    end
                end
                p(~isfinite(p))=NaN;
                out.Power=p; out.Phase=nan(size(p));
                if strcmp(op,'db_difference'), out.Quantity='power_ratio'; end
            end
            out=GrimData.derived(out,op);
        end

        function d = derived(d,operation)
            % Drop raw sample arrays and stale normalization claims after math.
            stale={'rcs_amp_real','rcs_amp_imag','raw_complex_amplitude_preserved', ...
                'storage_kind','storage_domain','rcs_domain','power_domain'};
            for k=1:numel(stale)
                if isfield(d.Metadata,stale{k}), d.Metadata=rmfield(d.Metadata,stale{k}); end
            end
            d.ID=''; d.Source='';
            d.Name=[d.Name ' | ' strrep(operation,'_',' ')];
            d.History=sprintf('%s\nMATLAB: %s',d.History,operation);
        end

        function d = regrid(d,axisName,values)
            names=GrimData.axes(); dim=find(strcmp(names,axisName));
            assert(~isempty(dim) && dim<4,'Choose a numeric axis.');
            values=unique(values(:)'); old=d.(axisName);
            assert(numel(old)>=2,'Interpolation needs at least two source coordinates.');
            assert(all(isfinite(values)) && ~isempty(values),'New coordinates must be finite.');
            order=[dim,setdiff(1:4,dim,'stable')];
            dims=GrimData.shape(d); shape=dims(order); shape(1)=numel(values);
            if any(isfinite(d.Phase(:)))
                z=sqrt(d.Power).*exp(1i*d.Phase); z(d.Power==0)=0;
                m=reshape(permute(z,order),numel(old),[]);
                z=interp1(old,m,values,'linear',NaN);
                z=ipermute(reshape(z,shape),order);
                % Unknown phase is not invented: interpolate its power only.
                p=reshape(permute(d.Power,order),numel(old),[]);
                p=ipermute(reshape(interp1(old,p,values,'linear',NaN),shape),order);
                phase=angle(z); valid=isfinite(z); p(valid)=abs(z(valid)).^2;
                d.Power=p; d.Phase=phase;
            else
                m=reshape(permute(d.Power,order),numel(old),[]);
                d.Power=ipermute(reshape(interp1(old,m,values,'linear',NaN),shape),order);
                d.Phase=nan(size(d.Power));
            end
            d.(axisName)=values;
            d=GrimData.derived(d,['regrid_' axisName]);
        end

        function d = reduce(d,dims,method,percentile,interpolate)
            if nargin<4, percentile=50; end
            if nargin<5, interpolate=false; end
            dims=unique(dims); other=setdiff(1:4,dims,'stable'); order=[dims,other];
            shape=GrimData.shape(d); n=prod(shape(dims));
            a=reshape(permute(d.Power,order),n,[]);
            switch method
                case 'mean', r=mean(a,1,'omitnan');
                case 'median', r=median(a,1,'omitnan');
                case 'min', r=min(a,[],1,'omitnan');
                case 'max', r=max(a,[],1,'omitnan');
                case 'std', r=std(a,1,1,'omitnan');
                case 'percentile', r=GrimData.percentile(a,percentile,interpolate);
                otherwise, error('Unknown statistic.');
            end
            r(all(isnan(a),1))=NaN;
            reducedShape=shape(order); reducedShape(1:numel(dims))=1;
            r=ipermute(reshape(r,reducedShape),order);
            repeats=ones(1,4); repeats(dims)=shape(dims);
            d.Power=repmat(r,repeats); d.Phase=nan(size(d.Power));
            d=GrimData.derived(d,method);
        end

        function r = percentile(a,p,interpolate)
            assert(isscalar(p) && isfinite(p) && p>=0 && p<=100,'Percentile must be from 0 to 100.');
            r=nan(1,size(a,2));
            for k=1:size(a,2)
                v=sort(a(isfinite(a(:,k)),k)); n=numel(v);
                if n==0, continue; end
                if interpolate
                    pos=1+(n-1)*p/100; lo=floor(pos); hi=ceil(pos);
                    r(k)=v(lo)+(pos-lo)*(v(hi)-v(lo));
                else
                    r(k)=v(max(1,ceil(n*p/100)));
                end
            end
        end

        function d = medianize(d,window,slide)
            assert(isfinite(window) && window>0 && isfinite(slide) && slide>0,'Window and slide must be positive.');
            x=d.Azimuth; assert(numel(x)>=2,'Medianize requires two or more azimuths.');
            periodic=x(end)-x(1)>=360-max(1.5*median(diff(x)),1e-6);
            if periodic
                assert(window<=360,'Periodic median window cannot exceed 360 degrees.');
                if abs(x(end)-x(1)-360)<1e-7
                    first=d.Power(1,:,:,:); last=d.Power(end,:,:,:);
                    first(isnan(first))=last(isnan(first)); d.Power(1,:,:,:)=first;
                    x=x(1:end-1); d.Power=d.Power(1:end-1,:,:,:);
                end
                centers=x(1):slide:(x(1)+360-slide*1e-9);
            else
                first=x(1)+window/2; last=x(end)-window/2;
                if last<first, centers=(x(1)+x(end))/2;
                else, centers=first:slide:last; end
            end
            shape=GrimData.shape(d); shape(1)=numel(centers); result=nan(shape);
            for k=1:numel(centers)
                distance=x-centers(k);
                if periodic, distance=mod(distance+180,360)-180; end
                use=abs(distance)<=window/2+1e-10;
                if any(use), result(k,:,:,:)=median(d.Power(use,:,:,:),1,'omitnan'); end
            end
            d.Azimuth=centers; d.Power=result; d.Phase=nan(size(result));
            d=GrimData.derived(d,'medianize');
        end

        function outputs = overlap(ds)
            assert(numel(ds)>=2,'Select at least two datasets.');
            names=GrimData.axes(); requested=struct();
            for dim=1:4
                v=ds{1}.(names{dim}); keep=true(size(v));
                for k=2:numel(ds)
                    source=ds{k}.(names{dim});
                    for j=1:numel(v)
                        if dim==4, found=sum(strcmp(source,v{j}));
                        else, found=sum(abs(source-v(j))<=1e-6); end
                        keep(j)=keep(j) && found==1;
                    end
                end
                requested.(names{dim})=v(keep);
            end
            outputs=cell(size(ds)); common=[];
            for k=1:numel(ds)
                outputs{k}=GrimData.select(ds{k},requested);
                if k==1, common=isfinite(outputs{k}.Power);
                else, common=common & isfinite(outputs{k}.Power); end
            end
            for k=1:numel(ds)
                outputs{k}.Power(~common)=NaN; outputs{k}.Phase(~common)=NaN;
                outputs{k}=GrimData.derived(outputs{k},'overlap');
            end
        end

        function d = join(ds)
            assert(numel(ds)>=2,'Select at least two datasets.');
            d=ds{1}; names=GrimData.axes();
            for k=2:numel(ds)
                assert(strcmp(d.Quantity,ds{k}.Quantity),'Join requires matching physical quantities.');
                assert(strcmp(GrimData.meta(d,'angular_coordinate_system','conic'), ...
                    GrimData.meta(ds{k},'angular_coordinate_system','conic')),'Coordinate conventions differ.');
                assert(strcmp(GrimData.meta(d,'elevation_coordinate_convention','signed_elevation'), ...
                    GrimData.meta(ds{k},'elevation_coordinate_convention','signed_elevation')),'Elevation conventions differ.');
                for dim=1:4
                    a=d.(names{dim}); b=ds{k}.(names{dim});
                    for j=1:numel(b)
                        if dim==4
                            if ~any(strcmp(a,b{j})), a{end+1}=b{j}; end
                        elseif ~any(abs(a-b(j))<=1e-6), a(end+1)=b(j);
                        end
                    end
                    if dim<4, a=sort(a); end
                    d.(names{dim})=a;
                end
            end
            d.Power=nan(GrimData.shape(d)); d.Phase=nan(size(d.Power));
            for k=1:numel(ds)
                next=ds{k}; ix=cell(1,4);
                for dim=1:4, ix{dim}=GrimData.match(d.(names{dim}),next.(names{dim})); end
                p=d.Power(ix{:}); ph=d.Phase(ix{:});
                both=isfinite(p)&isfinite(next.Power);
                tol=1e-10+1e-6*max(abs(p),abs(next.Power));
                assert(~any(both(:) & abs(p(:)-next.Power(:))>tol(:)), ...
                    'Conflicting overlapping samples. Crop or align datasets before joining.');
                phaseBoth=both & p>0 & isfinite(ph) & isfinite(next.Phase);
                assert(~any(abs(angle(exp(1i*(ph(phaseBoth)-next.Phase(phaseBoth)))))>1e-6), ...
                    'Conflicting phase in overlapping samples.');
                put=isnan(p); p(put)=next.Power(put); ph(put)=next.Phase(put);
                put=isnan(ph) & isfinite(next.Phase); ph(put)=next.Phase(put);
                d.Power(ix{:})=p; d.Phase(ix{:})=ph;
            end
            d=GrimData.derived(d,'join');
        end

        function d = transform(d,op,arg)
            if nargin<3, arg=[]; end
            switch op
                case 'duplicate'
                case 'flip_phase', d.Phase=-d.Phase;
                case 'offset', d.Power=d.Power*10^(arg/10);
                case 'shift'
                    d.Azimuth=d.Azimuth+arg(1); d.Elevation=d.Elevation+arg(2);
                    d.Phase=mod(d.Phase+deg2rad(arg(3))+pi,2*pi)-pi;
                case 'mirror', d.Azimuth=2*arg-d.Azimuth;
                case 'round'
                    d.Azimuth=GrimData.roundEven(d.Azimuth,arg(1));
                    d.Elevation=GrimData.roundEven(d.Elevation,arg(2));
                    d.Frequency=GrimData.roundEven(d.Frequency,arg(3));
                case 'wrap'
                    d.Azimuth=mod(d.Azimuth-arg,360)+arg;
                    seam=abs(d.Azimuth-arg)<1e-9 | abs(d.Azimuth-(arg+360))<1e-9;
                    d.Azimuth(seam)=arg;
                    [~,indices]=unique(d.Azimuth,'stable');
                    for j=reshape(indices,1,[])
                        repeats=find(abs(d.Azimuth-d.Azimuth(j))<=1e-9); repeats(repeats==j)=[];
                        for source=reshape(repeats,1,[])
                            p=d.Power(j,:,:,:); ph=d.Phase(j,:,:,:);
                            q=d.Power(source,:,:,:); qh=d.Phase(source,:,:,:);
                            missing=isnan(p); p(missing)=q(missing); ph(missing)=qh(missing);
                            same=isfinite(p)&isfinite(q)&abs(p-q)<=1e-10+1e-6*abs(q);
                            fill=same&isnan(ph)&isfinite(qh); ph(fill)=qh(fill);
                            d.Power(j,:,:,:)=p; d.Phase(j,:,:,:)=ph;
                        end
                    end
                    d=GrimData.slice(d,{indices,1:numel(d.Elevation),1:numel(d.Frequency),1:numel(d.Polarization)});
                case 'swap_el_az'
                    temp=d.Azimuth; d.Azimuth=d.Elevation; d.Elevation=temp;
                    d.Power=permute(d.Power,[2 1 3 4]); d.Phase=permute(d.Phase,[2 1 3 4]);
                case 'extrusion'
                    assert(arg>0 && isfinite(arg),'Extrusion length must be positive.');
                    scale=reshape(2*arg^2*d.Frequency*1e9/299792458,[1 1 numel(d.Frequency) 1]);
                    if strcmp(d.Quantity,'sigma_2d')
                        d.Power=d.Power.*scale; d.Quantity='sigma_3d';
                    elseif strcmp(d.Quantity,'sigma_3d')
                        d.Power=d.Power./scale; d.Quantity='sigma_2d';
                    else, error('Cannot convert a dimensionless ratio.'); end
                otherwise, error('Unknown transform %s',op);
            end
            d=GrimData.validate(d); d=GrimData.derived(d,op);
        end

        function result = roundEven(values,decimals)
            % NumPy rounds halfway values to even (MATLAB defaults away).
            scale=10^decimals; assert(isfinite(scale)&&scale>0,'Decimal count is out of range.');
            scaled=abs(values)*scale; integer=floor(scaled); fraction=scaled-integer;
            up=fraction>.5 | (fraction==.5 & mod(integer,2)==1);
            result=sign(values).*(integer+up)/scale;
        end
    end
end
