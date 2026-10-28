classdef GrimVendorIO
    % MATLAB ports of GRIM's Pioneer, PTM, OUT, SENTRi, CST and Xpatch readers.
    methods (Static)
        function d = readOut(path)
            lines=splitlines(string(fileread(path))); data=zeros(0,4);
            pol=GrimVendorIO.polarization(path);
            for k=1:numel(lines)
                line=strtrim(lines(k));
                if startsWith(lower(line),'# polarization:')
                    pol=char(strtrim(extractAfter(line,':'))); pol=strrep(pol,'"','');
                end
                line=regexprep(line,'#.*$',''); if strlength(strtrim(line))==0, continue; end
                row=sscanf(line,'%f')'; assert(numel(row)==4,'OUT line %d must have four numeric columns.',k);
                data(end+1,:)=row; %#ok<AGROW>
            end
            f=data(:,1); power=10.^(data(:,3)/10)./(2*pi*f*1e9/299792458);
            d=GrimIO.rows(data(:,2),zeros(size(f)),f,repmat({pol},numel(f),1),power,deg2rad(data(:,4)),'sigma_2d');
        end

        function pol = polarization(path)
            [~,stem]=fileparts(path); token=regexp(upper(stem),'(?<![A-Z0-9])(HH|VV|HV|VH)(?![A-Z0-9])','match','once');
            pol='NA'; if ~isempty(token), pol=token; end
        end

        function d = readPio(path)
            fid=fopen(path,'r'); assert(fid>=0,'Cannot open PIO file.'); cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
            h=struct();
            while true
                line=fgetl(fid); assert(ischar(line),'PIO has no Offset header.');
                token=regexp(line,'^\s*([^=]+)=(.*)$','tokens','once');
                if isempty(token), continue; end
                key=matlab.lang.makeValidName(lower(strtrim(token{1}))); h.(key)=strtrim(token{2});
                if strcmp(key,'offset'), break; end
            end
            offset=str2double(h.offset); nx=str2double(h.xsize); ny=str2double(h.ysize);
            assert(all(isfinite([offset nx ny])) && nx>0 && ny>0 && ...
                all([offset nx ny]==fix([offset nx ny])) && offset>=ftell(fid),'Invalid PIO framing.');
            endian='ieee-le'; if contains(lower(GrimIO.get(h,'order','little')),'big'), endian='ieee-be'; end
            precision=lower(h.precision); assert(ismember(precision,{'single','double'}),'Invalid PIO precision.');
            type=lower(GrimIO.get(h,'type','complex')); assert(ismember(type,{'real','complex'}),'Invalid PIO Type.');
            assert(strcmpi(GrimIO.get(h,'dataformat','binary'),'binary'),'PIO data must be binary.');
            parts=1+strcmp(type,'complex'); fseek(fid,offset,'bof');
            raw=fread(fid,parts*nx*ny,['*' precision],0,endian);
            assert(numel(raw)==parts*nx*ny,'Truncated PIO payload.');
            if parts==2, z=complex(double(raw(1:2:end)),double(raw(2:2:end)));
            else, z=complex(double(raw),zeros(size(raw))); z(isnan(raw))=complex(NaN,NaN); end
            valid=(isfinite(real(z))&isfinite(imag(z)))|(isnan(real(z))&isnan(imag(z)));
            assert(all(valid),'Invalid PIO complex sample.');
            z=reshape(z,[nx 1 ny 1]);
            while ~feof(fid)
                line=fgetl(fid); if ~ischar(line)||isempty(strtrim(line)), continue; end
                token=regexp(line,'^\s*([^=]+)=(.*)$','tokens','once'); assert(~isempty(token),'Invalid PIO footer.');
                key=matlab.lang.makeValidName(lower(strtrim(token{1})));
                if ~isfield(h,key), h.(key)=strtrim(token{2}); end
            end
            assert(ismember(lower(h.xname),{'azimuth','position'})&&strcmpi(h.yname,'frequency'),'Unsupported PIO axes.');
            az=GrimIO.angleDegrees(GrimVendorIO.pioAxis(h,'x',nx),h.xunits);
            f=GrimVendorIO.pioAxis(h,'y',ny)*GrimIO.frequencyGHz(h.yunits);
            el=GrimIO.angleDegrees(str2double(GrimIO.get(h,'elevation','0')),GrimIO.get(h,'elevationunits','deg'));
            pol=GrimIO.get(h,'polarity',GrimVendorIO.polarization(path));
            % A wrapped closing sample can repeat the initial coordinate.
            if numel(az)>1
                steps=mod(diff(az)+180,360)-180;
                steps(steps==-180 & diff(az)>0)=180;
                span=sum(steps); rawSpan=az(end)-az(1);
                tol=max(1e-7,8*double(eps('single'))*max(360,max(abs(az))));
                if abs(mod(rawSpan+180,360)-180)<=tol && ...
                        (abs(abs(span)-360)<=tol || abs(abs(rawSpan)-360)<=tol)
                    unwrapped=[az(1),az(1)+cumsum(steps)];
                    if any(diff(az(1:end-1))<=0)&&~all(diff(az(1:end-1))<0), az=unwrapped; end
                    az=az(1:end-1); z=z(1:end-1,:,:,:);
                end
            end
            d=GrimData.create(az,el,f,{pol},abs(z).^2,angle(z),'sigma_3d');
            d.History=GrimIO.get(h,'log','');
        end

        function values = pioAxis(h,prefix,count)
            key=[prefix 'vals'];
            if isfield(h,key), values=sscanf(regexprep(h.(key),'[:,]',' '),'%f')';
            else
                start=str2double(h.([prefix 'start']));
                if isfield(h,[prefix 'step'])
                    values=start+(0:count-1)*str2double(h.([prefix 'step']));
                elseif count==1, values=start;
                else, values=linspace(start,str2double(h.([prefix 'stop'])),count); end
            end
            assert(numel(values)==count && all(isfinite(values)),'Invalid PIO axis.');
        end

        function writePio(path,d)
            GrimVendorIO.checkCut(d); z=GrimData.field(d);
            az=sprintf('%.17g:',d.Azimuth); freq=sprintf('%.17g:',d.Frequency);
            header=sprintf(['Pioneer\nType=Complex\nDataFormat=Binary\nPrecision=Double\nOrder=Little Endian\n' ...
                'XName=Azimuth\nXUnits=deg\nXSize=%d\nXVals=%s\nYName=Frequency\nYUnits=GHz\nYSize=%d\nYVals=%s\n' ...
                'Elevation=%.17g\nElevationUnits=deg\nPolarity=%s\n'], ...
                numel(d.Azimuth),az(1:end-1),numel(d.Frequency),freq(1:end-1),d.Elevation,d.Polarization{1});
            offset=numel(header)+numel(sprintf('Offset=%012d\n',0));
            header=[header sprintf('Offset=%012d\n',offset)];
            fid=fopen(path,'w','ieee-le'); assert(fid>=0,'Cannot create PIO.'); cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid,uint8(header),'uint8'); fwrite(fid,[real(z(:))';imag(z(:))'],'double');
        end

        function checkCut(d)
            assert(numel(d.Elevation)==1 && numel(d.Polarization)==1,'Export one elevation and one polarization per PIO/PTM file.');
            assert(strcmp(d.Quantity,'sigma_3d'),'PIO/PTM export requires 3-D RCS (dBsm).');
        end

        function d = readPtm(path)
            info=dir(path); fid=fopen(path,'r'); assert(fid>=0,'Cannot open PTM.'); cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
            assert(info.bytes>=291,'Truncated PTM header.'); valid=false;
            for machine={'ieee-le','ieee-be'}
                fseek(fid,0,'bof'); integers=fread(fid,9,'int32',0,machine{1}); values=fread(fid,9,'single',0,machine{1});
                n=integers(1); nf=info.bytes/(8*(n+1));
                valid=n>=1 && nf>=37 && nf==fix(nf) && all(isfinite(values)) && ...
                    values(3)>0 && values(4)>=0 && values(3)-values(4)/2>0;
                if valid, endian=machine{1}; break; end
            end
            assert(valid,'Invalid PTM header/framing.');
            fseek(fid,72,'bof'); pol=strtrim(char(fread(fid,2,'*uint8')'));
            fseek(fid,124,'bof'); configuration=char(fread(fid,50,'*uint8')');
            fseek(fid,8*nf,'bof'); raw=fread(fid,2*n*nf,'single',0,endian);
            assert(numel(raw)==2*n*nf && all(isfinite(raw)),'Invalid PTM samples.');
            z=reshape(complex(raw(1:2:end),raw(2:2:end)),nf,n).';
            az=mod(values(1)+(0:n-1)*values(2)+180,360)-180;
            f=linspace(values(3)-values(4)/2,values(3)+values(4)/2,nf);
            d=GrimData.create(az,values(6),f,{pol},reshape(abs(z).^2,[n 1 nf 1]), ...
                reshape(angle(z),[n 1 nf 1]),'sigma_3d');
            d.Metadata.time_convention='exp(+jwt)';
            d.Metadata.ptm_configuration=strtrim(strrep(configuration,'*',''));
        end

        function writePtm(path,d)
            GrimVendorIO.checkCut(d); nf=numel(d.Frequency); n=numel(d.Azimuth);
            assert(nf>=37,'PTM requires at least 37 frequencies for its header block.');
            assert(ismember(d.Polarization{1},{'HH','VV','HV','VH'}),'Unsupported PTM polarization.');
            delta=0; if n>1, delta=mean(diff(d.Azimuth)); end
            assert(n<3 || max(abs(diff(d.Azimuth)-delta))<1e-6,'PTM needs a uniformly spaced azimuth axis.');
            assert(max(abs(diff(d.Frequency)-mean(diff(d.Frequency))))<1e-6,'PTM needs a uniformly spaced frequency axis.');
            assert(n<2 || d.Azimuth(end)-d.Azimuth(1)<360-1e-6,'Remove the repeated full-turn endpoint before PTM export.');
            z=GrimData.field(d);
            convention=lower(GrimData.meta(d,'time_convention','exp(+jwt)'));
            if contains(convention,'-jwt'), z=conj(z);
            else, assert(contains(convention,'+jwt'),'Unrecognized PTM time convention.'); end
            assert(all(isfinite(z(:))) && all(abs(real(z(:)))<=realmax('single')) && ...
                all(abs(imag(z(:)))<=realmax('single')),'PTM cannot store missing or overflowing samples.');
            fid=fopen(path,'w','ieee-le'); assert(fid>=0,'Cannot create PTM.'); cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid,[n zeros(1,8)],'int32');
            fwrite(fid,[d.Azimuth(1),delta,mean(d.Frequency([1 end])),diff(d.Frequency([1 end])),0,d.Elevation,0,0,0],'single');
            fwrite(fid,uint8(d.Polarization{1}),'uint8');
            GrimVendorIO.padded(fid,d.Name,50); GrimVendorIO.padded(fid,'GRIM_GC_V1',50);
            GrimVendorIO.padded(fid,'MATLAB GRIM',46); fwrite(fid,nf,'int32');
            [~,name,ext]=fileparts(path); GrimVendorIO.padded(fid,[name ext],50);
            GrimVendorIO.padded(fid,upper(char(datetime('now','Format','ddMMMyyyy'))),9);
            GrimVendorIO.padded(fid,char(datetime('now','Format','HH:mm:ss')),8);
            fwrite(fid,zeros(8*nf-291,1,'uint8'),'uint8');
            z=reshape(z,n,nf).'; fwrite(fid,[real(z(:))';imag(z(:))'],'single');
        end

        function padded(fid,text,width)
            bytes=uint8(char(text)); bytes(bytes>127)=uint8('?'); bytes=bytes(1:min(width,end));
            fwrite(fid,[bytes repmat(uint8('*'),1,width-numel(bytes))],'uint8');
        end

        function cells = textCells(path)
            text=fileread(path); lines=splitlines(string(text)); delimiter=',';
            for k=1:min(100,numel(lines))
                if contains(lower(lines(k)),{'freq','theta','azimuth'})
                    if contains(lines(k),sprintf('\t')), delimiter=sprintf('\t');
                    elseif contains(lines(k),';'), delimiter=';'; end
                    break
                end
            end
            cells=readcell(path,'FileType','text','Delimiter',delimiter);
        end

        function value = cellText(c)
            if isempty(c)||(isnumeric(c)&&isscalar(c)&&isnan(c)), value='';
            else, value=char(string(c)); end
        end

        function d = readSentri(path)
            cells=GrimVendorIO.textCells(path); header=0; map=struct(); scale=1e-9;
            for r=1:size(cells,1)
                map=struct();
                for j=1:size(cells,2)
                    token=lower(regexprep(GrimVendorIO.cellText(cells{r,j}),'[^a-zA-Z0-9]',''));
                    if ismember(token,{'freqmhz','frequency'}), map.f=j; if strcmp(token,'freqmhz'), scale=1e-3; end
                    elseif ismember(token,{'theta','thetadeg'}), map.el=j;
                    elseif ismember(token,{'phi','phideg'}), map.az=j;
                    else
                        suffix='';
                        if startsWith(token,{'rcspp','efieldphasepp'}) || contains(token,{'phiscatphiinc','phiincphiscat','phiphi'}), suffix='HH';
                        elseif startsWith(token,{'rcstt','efieldphasett'}) || contains(token,{'thetascatthetainc','thetaincthetascat','thetatheta'}), suffix='VV';
                        elseif startsWith(token,{'rcspt','efieldphasept'}) || contains(token,{'phiscatthetainc','thetaincphiscat','phitheta'}), suffix='HV';
                        elseif startsWith(token,{'rcstp','efieldphasetp'}) || contains(token,{'thetascatphiinc','phiincthetascat','thetaphi'}), suffix='VH'; end
                        if ~isempty(suffix)
                            if contains(token,'phase'), map.(['phase' suffix])=j;
                            elseif startsWith(token,'rcs'), map.(['power' suffix])=j; end
                        end
                    end
                end
                if all(isfield(map,{'f','el','az','powerHH','powerVV','powerHV','powerVH'})), header=r; break; end
            end
            assert(header>0,'Unrecognized SENTRi header.');
            [az,el,f,pol,p,ph]=GrimVendorIO.wideRows(cells,header,map,scale);
            if ~(all(az>=0)&&all(az<=180)), az=mod(az+180,360)-180; end
            d=GrimIO.rows(az,el,f,pol,p,ph,'sigma_3d');
            d.Metadata.angular_coordinate_system='conic';
            d.Metadata.elevation_coordinate_convention='sentri_theta_top_zero';
            d.Metadata.sentri_elevation_convention='sentri_theta_top_zero';
            d.Metadata.time_convention='exp(+jwt)';
        end

        function d = readCst(path)
            cells=GrimVendorIO.textCells(path); header=0; map=struct(); wide=false; scale=1;
            for r=1:size(cells,1)
                map=struct(); wide=false;
                for j=1:size(cells,2)
                    token=lower(regexprep(GrimVendorIO.cellText(cells{r,j}),'[^a-zA-Z0-9]',''));
                    if startsWith(token,{'frequency','freq'})
                        map.f=j; unit=regexp(token,'(ghz|mhz|khz|hz)','match','once');
                        if isempty(unit), unit='ghz'; end
                        scale=GrimIO.frequencyGHz(unit);
                    elseif startsWith(token,'elevation'), map.el=j;
                    elseif startsWith(token,'azimuth'), map.az=j;
                    elseif startsWith(token,'theta'), map.el=j; wide=true;
                    elseif startsWith(token,'phi'), map.az=j;
                    elseif ismember(token,{'pol','polarity','polarization'}), map.pol=j;
                    elseif contains(token,{'rcs','magnitude','phase'})
                        suffix='';
                        if contains(token,'thetatheta'), suffix='VV';
                        elseif contains(token,'phitheta'), suffix='HV';
                        elseif contains(token,'thetaphi'), suffix='VH';
                        elseif contains(token,'phiphi'), suffix='HH'; end
                        if contains(token,'phase')
                            assert(~contains(token,'rad'),'CST phase must be declared in degrees.');
                            map.(['phase' suffix])=j;
                        elseif contains(token,{'dbsm','dbm2'})
                            map.(['power' suffix])=j;
                        end
                    elseif ismember(token,{'iq','complexiq','complexsample','complexamplitude'}), map.iq=j;
                    end
                end
                if all(isfield(map,{'f','el','az'})) && ...
                        ((wide&&isfield(map,'powerVV')) || (isfield(map,'pol')&&(isfield(map,'power')||isfield(map,'iq'))))
                    header=r; break
                end
            end
            assert(header>0,'Unrecognized CST table header.');
            if wide
                [az,theta,f,pol,p,ph]=GrimVendorIO.wideRows(cells,header,map,scale); el=90-theta;
            else
                az=[]; el=[]; f=[]; pol={}; p=[]; ph=[];
                for r=header+1:size(cells,1)
                    frequency=str2double(string(cells{r,map.f})); if isnan(frequency), continue; end
                    f(end+1,1)=frequency*scale; az(end+1,1)=str2double(string(cells{r,map.az})); %#ok<AGROW>
                    el(end+1,1)=str2double(string(cells{r,map.el})); pol{end+1,1}=GrimVendorIO.cellText(cells{r,map.pol}); %#ok<AGROW>
                    power=NaN; phase=NaN;
                    if isfield(map,'power'), power=10^(str2double(string(cells{r,map.power}))/10); end
                    if isfield(map,'phase'), phase=deg2rad(str2double(string(cells{r,map.phase}))); end
                    if isfield(map,'iq')
                        token=lower(GrimVendorIO.cellText(cells{r,map.iq}));
                        token=regexprep(token,'[\s()\[\]{}*]',''); token=strrep(token,'j','i');
                        token=regexprep(token,'(?<=\d)d(?=[+-]?\d)','e');
                        z=str2double(token);
                        if isfinite(z)
                            if isfinite(power), assert(abs(power-abs(z)^2)<=5e-6*power,'CST IQ/magnitude mismatch.'); end
                            if isfinite(phase)&&abs(z)>0, assert(abs(angle(exp(1i*(phase-angle(z)))))<=1e-6,'CST IQ/phase mismatch.'); end
                            power=abs(z)^2; phase=angle(z);
                        else
                            assert(~isnan(power),'CST row has neither a usable magnitude nor complex IQ.');
                        end
                    end
                    p(end+1,1)=power; ph(end+1,1)=phase; %#ok<AGROW>
                end
            end
            az=mod(az+180,360)-180;
            d=GrimIO.rows(az,el,f,pol,p,ph,'sigma_3d');
        end

        function [az,el,f,pol,p,ph] = wideRows(cells,header,map,scale)
            az=[]; el=[]; f=[]; pol={}; p=[]; ph=[];
            labels={'VV','HV','VH','HH'};
            for r=header+1:size(cells,1)
                frequency=str2double(string(cells{r,map.f})); if isnan(frequency), continue; end
                for k=1:4
                    label=labels{k}; if ~isfield(map,['power' label]), continue; end
                    az(end+1,1)=str2double(string(cells{r,map.az})); %#ok<AGROW>
                    el(end+1,1)=str2double(string(cells{r,map.el})); f(end+1,1)=frequency*scale; %#ok<AGROW>
                    pol{end+1,1}=label; p(end+1,1)=10^(str2double(string(cells{r,map.(['power' label])}))/10); %#ok<AGROW>
                    phase=NaN;
                    if isfield(map,['phase' label]), phase=deg2rad(str2double(string(cells{r,map.(['phase' label])}))); end
                    ph(end+1,1)=phase; %#ok<AGROW>
                end
            end
        end

        function d = readSS(path)
            fid=fopen(path,'r','ieee-be'); assert(fid>=0,'Cannot open SS file.'); cleaner=onCleanup(@() fclose(fid)); %#ok<NASGU>
            info=dir(path); nb=fread(fid,1,'int32'); nd=fread(fid,1,'int32'); nf=(nd-408)/32;
            assert(nb>=648&&nf>0&&nf==fix(nf)&&nb+nd<=info.bytes,'Invalid SS framing.');
            candidates=[];
            for offset=648+(0:3)*256
                if offset+588<=nb && GrimVendorIO.readAt(fid,offset+80,'int32')==nf && ...
                        ismember(GrimVendorIO.readAt(fid,offset+424,'int32'),[1 2]) && ...
                        ismember(GrimVendorIO.readAt(fid,offset+456,'int32'),[1 2])
                    candidates(end+1)=offset; %#ok<AGROW>
                end
            end
            assert(numel(candidates)==1,'SS header-C is missing or ambiguous.'); c=candidates(1);
            ifreq=GrimVendorIO.readAt(fid,c+424,'int32'); imono=GrimVendorIO.readAt(fid,c+456,'int32');
            f1=GrimVendorIO.readAt(fid,c+428,'single'); f2=GrimVendorIO.readAt(fid,c+432,'single');
            if ifreq==2, fseek(fid,nb-4*nf,'bof'); freq=fread(fid,nf,'single')';
            else, freq=linspace(f1,f2,nf); end
            count=info.bytes/(nb+nd); assert(count==fix(count),'Truncated SS record.');
            angles=zeros(count,4); samples=complex(zeros(count,nf,4));
            for k=1:count
                pos=(k-1)*(nb+nd);
                assert(GrimVendorIO.readAt(fid,pos,'int32')==nb && ...
                    GrimVendorIO.readAt(fid,pos+4,'int32')==nd,'Inconsistent SS record framing.');
                assert(GrimVendorIO.readAt(fid,pos+c+80,'int32')==nf && ...
                    GrimVendorIO.readAt(fid,pos+c+424,'int32')==ifreq && ...
                    GrimVendorIO.readAt(fid,pos+c+456,'int32')==imono,'Inconsistent SS header-C.');
                fseek(fid,pos+nb+280,'bof'); angles(k,:)=fread(fid,4,'single')';
                fseek(fid,pos+nb+408,'bof'); values=fread(fid,8*nf,'single');
                assert(numel(values)==8*nf,'Truncated SS samples.');
                good=(isfinite(values(1:2:end))&isfinite(values(2:2:end))) | ...
                    (isnan(values(1:2:end))&isnan(values(2:2:end)));
                assert(all(good),'SS contains invalid complex samples.');
                z=reshape(complex(values(1:2:end),values(2:2:end)),4,nf).'; samples(k,:,:)=reshape(z,[1 nf 4]);
            end
            variation=@(v) numel(unique(round(mod(v,360),4)));
            ni=max(variation(angles(:,1)),variation(angles(:,2)));
            no=max(variation(angles(:,3)),variation(angles(:,4)));
            assert(~(imono==2&&ni>1&&no>1),'Bistatic SS varying both angular pairs cannot fit a four-axis grid.');
            cols=[1 2]; seamOffsets=[472 476 480];
            if no>ni, cols=[3 4]; seamOffsets=[496 500 504]; end
            az=angles(:,cols(1)); el=angles(:,cols(2));
            start=GrimVendorIO.readAt(fid,c+seamOffsets(1),'single'); stop=GrimVendorIO.readAt(fid,c+seamOffsets(2),'single');
            steps=GrimVendorIO.readAt(fid,c+seamOffsets(3),'int32');
            if abs(stop-start-360)<5e-5 && steps>1
                % Preserve a declared closing +180 rather than colliding with -180.
                for k=2:numel(az)
                    if abs(az(k)-start)<5e-5 && any(abs(az(1:k-1)-start)<5e-5 & abs(el(1:k-1)-el(k))<5e-5)
                        az(k)=stop;
                    end
                end
            end
            az=round(az,4); el=round(el,4);
            assert(size(unique([az el],'rows'),1)==count,'SS angular coordinate collision.');
            [rr,ff,pp]=ndgrid(1:count,1:nf,1:4); labels={'VV','VH','HV','HH'};
            d=GrimIO.rows(az(rr(:)),el(rr(:)),reshape(freq(ff(:)),[],1), ...
                reshape(labels(pp(:)),[],1),abs(samples(:)).^2,angle(samples(:)),'sigma_3d');
            d.Metadata.measurement_geometry='monostatic';
            if imono==2, d.Metadata.measurement_geometry='bistatic'; end
        end

        function value = readAt(fid,offset,precision)
            fseek(fid,offset,'bof'); value=fread(fid,1,precision);
            assert(~isempty(value),'Truncated binary header.');
        end
    end
end
