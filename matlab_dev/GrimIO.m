classdef GrimIO
    % MATLAB-native file dispatch and GRIM native/flat-table interchange.
    % Ported from GRIM_Backend/io/native.py, csv.py, loaders.py.
    methods (Static)
        function d = load(path)
            [~,name,ext]=fileparts(path);
            switch lower(ext)
                case {'.grim','.npz'}, d=GrimIO.fromNative(GrimNpy.readArchive(path));
                case '.mat'
                    content=builtin('load',path);
                    assert(isfield(content,'dataset'),'MAT file must contain a dataset saved by this app.');
                    d=GrimData.validate(content.dataset);
                case {'.pio','.cmplx_di'}, d=GrimVendorIO.readPio(path);
                case '.ptm', d=GrimVendorIO.readPtm(path);
                case '.out', d=GrimVendorIO.readOut(path);
                case '.ss', d=GrimVendorIO.readSS(path);
                otherwise, d=GrimIO.readTable(path);
            end
            d.Name=name; d.Source=path; d.ID='';
        end

        function save(path,d)
            d=GrimData.validate(d); [folder,~,ext]=fileparts(path);
            if isempty(folder), folder=pwd; end
            assert(isfolder(folder),'Output directory does not exist.');
            % Complete the file before replacing the chosen destination.
            temporary=[tempname(folder) ext]; cleanup=onCleanup(@() GrimIO.removeTemporary(temporary)); %#ok<NASGU>
            switch lower(ext)
                case '.mat'
                    dataset=d; builtin('save',temporary,'dataset','-v7.3');
                case '.grim', GrimNpy.writeArchive(temporary,GrimIO.toNative(d));
                case {'.csv','.tsv'}, GrimIO.writeTable(temporary,d);
                case '.pio', GrimVendorIO.writePio(temporary,d);
                case '.ptm', GrimVendorIO.writePtm(temporary,d);
                otherwise, error('Use .grim, .mat, .csv, .tsv, .pio or .ptm.');
            end
            [ok,msg]=movefile(temporary,path,'f'); assert(ok,'%s',msg);
        end

        function removeTemporary(path)
            if isfile(path), delete(path); end
        end

        function d = fromNative(a)
            required={'azimuths','elevations','frequencies','polarizations'};
            assert(all(isfield(a,required)),'Missing native dataset axes.');
            units=struct();
            if isfield(a,'units') && strlength(string(a.units))>0, units=jsondecode(char(a.units)); end
            az=GrimIO.angleDegrees(a.azimuths,GrimIO.get(units,'azimuth','deg'));
            el=GrimIO.angleDegrees(a.elevations,GrimIO.get(units,'elevation','deg'));
            f=double(a.frequencies)*GrimIO.frequencyGHz(GrimIO.get(units,'frequency','GHz'));
            quantity=GrimIO.get(units,'rcs_linear_quantity','');
            if isempty(quantity)
                logUnit=GrimIO.get(units,'rcs_log_unit','dBsm');
                switch lower(logUnit)
                    case 'dbke', quantity='sigma_2d';
                    case 'db', quantity='power_ratio';
                    otherwise, quantity='sigma_3d';
                end
            end
            version=double(GrimIO.get(a,'grim_format_version',1));
            assert(ismember(version,[1 2]),'Unsupported GRIM archive version.');
            if version==2
                kind=GrimIO.get(a,'storage_kind',''); domain=GrimIO.get(a,'storage_domain','');
                assert(ismember(kind,{'complex','power','power_phase'}),'Invalid native storage kind.');
                assert(ismember(domain,{'physical','raw_field'}) && ...
                    (~strcmp(domain,'raw_field')||strcmp(kind,'complex')),'Invalid native storage domain.');
                actual=isfield(a,{'rcs_complex','rcs_power','rcs_phase'});
                expected=[strcmp(kind,'complex'),~strcmp(kind,'complex'),strcmp(kind,'power_phase')];
                assert(isequal(actual,expected),'Native sample arrays do not match storage_kind.');
            end
            legacy=strcmp(GrimIO.get(a,'power_domain',''),'delta_amp_sq') && ...
                strcmp(GrimIO.get(a,'rcs_domain',''),'delta');
            if isfield(a,'rcs_complex')
                z=double(a.rcs_complex);
                valid=(isfinite(real(z))&isfinite(imag(z)))|(isnan(real(z))&isnan(imag(z)));
                assert(all(valid(:)),'Invalid complex samples in archive.');
                if strcmp(GrimIO.get(a,'storage_domain','physical'),'raw_field') && ~legacy
                    if strcmp(quantity,'sigma_2d')
                        z=z./reshape(2*sqrt(2*pi*f*1e9/299792458),[1 1 numel(f) 1]);
                    elseif strcmp(quantity,'sigma_3d'), z=z*sqrt(4*pi);
                    else, error('Raw field requires a physical RCS quantity.'); end
                end
                p=abs(z).^2; phase=angle(z);
            elseif legacy && all(isfield(a,{'rcs_amp_real','rcs_amp_imag'}))
                z=complex(double(a.rcs_amp_real),double(a.rcs_amp_imag)); p=abs(z).^2; phase=angle(z);
            else
                assert(isfield(a,'rcs_power'),'Missing native samples.');
                p=double(a.rcs_power); phase=GrimIO.get(a,'rcs_phase',nan(size(p)));
            end
            d=GrimData.create(az,el,f,a.polarizations,p,phase,quantity);
            reserved=[required,{'units','rcs_power','rcs_phase','rcs_complex','source_path','history', ...
                'grim_format_version','storage_kind','storage_domain','rcs_amp_real','rcs_amp_imag', ...
                'raw_complex_amplitude_preserved','rcs_domain','power_domain'}];
            names=setdiff(fieldnames(a),reserved);
            for k=1:numel(names), d.Metadata.(names{k})=a.(names{k}); end
            names=setdiff(fieldnames(units),{'azimuth','elevation','frequency','rcs_log_unit','rcs_linear_quantity'});
            for k=1:numel(names), d.Metadata.(names{k})=units.(names{k}); end
            d.History=char(GrimIO.get(a,'history',''));
        end

        function a = toNative(d)
            a=d.Metadata;
            names=fieldnames(a);
            for k=1:numel(names)
                value=a.(names{k});
                if isstruct(value), a.(names{k})=jsonencode(value); end
            end
            a.azimuths=d.Azimuth; a.elevations=d.Elevation; a.frequencies=d.Frequency;
            a.polarizations=d.Polarization; a.rcs_power=d.Power;
            a.grim_format_version=int16(2); a.storage_domain='physical';
            if any(isfinite(d.Phase(:)))
                a.rcs_phase=d.Phase; a.storage_kind='power_phase'; a.rcs_domain='power_phase';
            else
                a.storage_kind='power'; a.rcs_domain='linear_rcs';
            end
            a.power_domain='linear_rcs'; a.source_path=d.Source; a.history=d.History;
            units=struct('azimuth','deg','elevation','deg','frequency','GHz', ...
                'rcs_log_unit',GrimData.logUnit(d),'rcs_linear_quantity',d.Quantity);
            keys={'angular_coordinate_system','elevation_coordinate_convention'};
            for k=1:numel(keys)
                if isfield(d.Metadata,keys{k}), units.(keys{k})=d.Metadata.(keys{k}); end
            end
            a.units=jsonencode(units);
        end

        function d = readTable(path)
            % Vendor tables need their own header interpretation.
            text=fileread(path); compact=lower(regexprep(text(1:min(end,32000)),'[^a-zA-Z0-9]',''));
            if contains(compact,'rcsppdbsm') || contains(compact,'rcsphiscatphiinc') || contains(compact,'rcsphiincphiscat')
                d=GrimVendorIO.readSentri(path); return
            end
            if contains(compact,'thetaphi') || contains(compact,'thetadegphideg') || endsWith(lower(path),'.cst_data')
                d=GrimVendorIO.readCst(path); return
            end
            opts=detectImportOptions(path,'FileType','text','VariableNamingRule','preserve');
            t=readtable(path,opts); names=lower(string(t.Properties.VariableNames));
            required={'azimuth','elevation','frequency','polarization'};
            if ~all(ismember(required,names))
                d=GrimVendorIO.readCst(path); return
            end
            col=@(name) t.(t.Properties.VariableNames{find(names==name,1)});
            az=double(col('azimuth')); el=double(col('elevation')); f=double(col('frequency'));
            pol=cellstr(string(col('polarization')));
            if ismember('frequency_unit',names), unit=GrimIO.singleText(col('frequency_unit'));
            else
                unit='GHz'; if median(f,'omitnan')>=1e6, unit='Hz'; end
                warning('GrimIO:LegacyUnits','No frequency_unit column; interpreting frequency as %s.',unit);
            end
            f=f*GrimIO.frequencyGHz(unit);
            if ismember('azimuth_unit',names), az=GrimIO.angleDegrees(az,GrimIO.singleText(col('azimuth_unit'))); end
            if ismember('elevation_unit',names), el=GrimIO.angleDegrees(el,GrimIO.singleText(col('elevation_unit'))); end
            quantity='sigma_3d';
            if ismember('rcs_linear_quantity',names), quantity=GrimIO.singleText(col('rcs_linear_quantity'));
            elseif ismember('magnitude_unit',names)
                unit=lower(GrimIO.singleText(col('magnitude_unit')));
                if ismember(unit,{'dbke','m'}), quantity='sigma_2d';
                elseif ismember(unit,{'db','1'}), quantity='power_ratio'; end
            elseif ismember('magnitude_dbke',names), quantity='sigma_2d'; end
            candidates={};
            magnitudeNames={'magnitude_power_linear','magnitude_linear','magnitude_dbsm','magnitude_dbke','magnitude_db','magnitude'};
            for k=1:numel(magnitudeNames)
                name=magnitudeNames{k}; if ~ismember(name,names), continue; end
                p=double(col(name));
                if strcmp(name,'magnitude')
                    assert(ismember('magnitude_unit',names),'magnitude requires a magnitude_unit column.');
                    unit=lower(GrimIO.singleText(col('magnitude_unit')));
                else, unit=name; end
                if contains(unit,'db')
                    p=10.^(p/10);
                    if contains(unit,'dbke'), p=p./(2*pi*f*1e9/299792458); end
                end
                candidates{end+1}=p; %#ok<AGROW>
            end
            assert(~isempty(candidates),'No supported magnitude column.'); p=candidates{1};
            for k=2:numel(candidates)
                q=candidates{k}; both=isfinite(p)&isfinite(q);
                assert(~any(abs(p(both)-q(both))>5e-6*abs(q(both))),'Redundant magnitude columns disagree.');
                missing=isnan(p); p(missing)=q(missing);
            end
            phase=nan(size(p));
            if ismember('phase_rad',names), phase=double(col('phase_rad'));
            elseif ismember('phase_deg',names), phase=deg2rad(double(col('phase_deg')));
            elseif ismember('phase',names)
                phase=double(col('phase'));
                unit='deg'; if ismember('phase_unit',names), unit=GrimIO.singleText(col('phase_unit')); end
                phase=deg2rad(GrimIO.angleDegrees(phase,unit));
            end
            d=GrimIO.rows(az,el,f,pol,p,phase,quantity);
            if ismember('angular_coordinate_system',names)
                d.Metadata.angular_coordinate_system=GrimIO.singleText(col('angular_coordinate_system'));
            end
        end

        function d = rows(az,el,f,pol,p,phase,quantity)
            assert(~isempty(p),'File contains no samples.');
            assert(all(isfinite(az)&isfinite(el)&isfinite(f)),'Invalid table coordinates.');
            [a,~,ia]=unique(az); [e,~,ie]=unique(el); [freq,~,ifr]=unique(f);
            labels=unique(upper(strtrim(string(pol))),'stable');
            [~,ip]=ismember(upper(strtrim(string(pol))),labels);
            shape=[numel(a),numel(e),numel(freq),numel(labels)];
            P=nan(shape); Ph=nan(shape); indices=sub2ind(shape,ia,ie,ifr,ip);
            % Merge duplicate rows only if finite values are equivalent.
            for k=1:numel(p)
                j=indices(k);
                if isfinite(P(j)) && isfinite(p(k))
                    assert(abs(P(j)-p(k))<=1e-10+5e-6*abs(p(k)),'Conflicting duplicate table rows.');
                    if isfinite(Ph(j)) && isfinite(phase(k)) && p(k)>0
                        assert(abs(angle(exp(1i*(Ph(j)-phase(k)))))<1e-6,'Duplicate phase values disagree.');
                    end
                elseif isnan(P(j)), P(j)=p(k); end
                if isnan(Ph(j)), Ph(j)=phase(k); end
            end
            d=GrimData.create(a,e,freq,labels,P,Ph,quantity);
        end

        function writeTable(path,d)
            [a,e,f,p]=ndgrid(d.Azimuth,d.Elevation,d.Frequency,1:numel(d.Polarization));
            n=numel(a); repeat=@(s) repmat(string(s),n,1);
            t=table(repeat('grim.flat-rcs.v1'),a(:),repeat('deg'),e(:),repeat('deg'),f(:),repeat('GHz'), ...
                reshape(string(d.Polarization(p(:))),[],1),repeat(d.Quantity),repeat(GrimData.logUnit(d)), ...
                repeat(GrimData.meta(d,'angular_coordinate_system','conic')),d.Power(:),d.Phase(:), ...
                'VariableNames',{'grim_csv_schema','azimuth','azimuth_unit','elevation','elevation_unit', ...
                'frequency','frequency_unit','polarization','rcs_linear_quantity','rcs_log_unit', ...
                'angular_coordinate_system','magnitude_power_linear','phase_rad'});
            delimiter=','; if endsWith(lower(path),'.tsv'), delimiter='\t'; end
            writetable(t,path,'FileType','text','Delimiter',delimiter,'Encoding','UTF-8');
        end

        function value = get(s,key,fallback)
            value=fallback; if isfield(s,key), value=s.(key); end
        end

        function value = singleText(column)
            values=unique(strtrim(string(column))); values=values(~ismissing(values)&strlength(values)>0);
            assert(numel(values)==1,'A metadata/unit column has inconsistent values.'); value=char(values);
        end

        function scale = frequencyGHz(unit)
            switch lower(strtrim(char(unit)))
                case 'hz', scale=1e-9;
                case 'khz', scale=1e-6;
                case 'mhz', scale=1e-3;
                case 'ghz', scale=1;
                otherwise, error('Unsupported frequency unit: %s',unit);
            end
        end

        function degrees = angleDegrees(values,unit)
            switch lower(strtrim(char(unit)))
                case {'deg','degree','degrees'}, degrees=double(values);
                case {'rad','radian','radians'}, degrees=rad2deg(double(values));
                otherwise, error('Unsupported angle unit: %s',unit);
            end
        end
    end
end
