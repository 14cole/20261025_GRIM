classdef GrimNpy
    % NumPy NPY/NPZ serialization implemented in MATLAB, without Python.
    % Supports numeric, complex, logical, byte-string and Unicode arrays.
    % Object/pickle and structured dtypes are deliberately rejected.
    methods (Static)
        function data = readArchive(path)
            folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
            % Check names before extraction so archive members stay in temp.
            archive=java.util.zip.ZipFile(java.io.File(path));
            closer=onCleanup(@() archive.close()); %#ok<NASGU>
            entries=archive.entries();
            while entries.hasMoreElements()
                entry=entries.nextElement(); name=char(entry.getName());
                assert(~entry.isDirectory() && ~contains(name,{'/','\','..'}) && endsWith(name,'.npy'), ...
                    'Native archive must contain flat .npy members only.');
            end
            unzip(path,folder); files=dir(fullfile(folder,'*.npy')); data=struct();
            for k=1:numel(files)
                [~,name]=fileparts(files(k).name);
                assert(isvarname(name),'Invalid native archive member name.');
                data.(name)=GrimNpy.read(fullfile(folder,files(k).name));
            end
        end

        function writeArchive(path,data)
            folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
            names=fieldnames(data); files=cell(size(names));
            for k=1:numel(names)
                name=names{k}; value=data.(name); files{k}=[name '.npy'];
                shape=[];
                if ismember(name,{'rcs_power','rcs_phase','rcs_complex'})
                    shape=[numel(data.azimuths),numel(data.elevations),numel(data.frequencies),numel(data.polarizations)];
                end
                GrimNpy.write(fullfile(folder,files{k}),value,shape);
            end
            zipPath=fullfile(folder,'dataset.zip'); zip(zipPath,files,folder);
            [ok,msg]=copyfile(zipPath,path,'f'); assert(ok,'%s',msg);
        end

        function value = read(path)
            fid=fopen(path,'r','ieee-le'); assert(fid>=0,'Cannot open NPY file.');
            cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
            magic=fread(fid,6,'*uint8')';
            assert(isequal(magic,uint8([147 double('NUMPY')])),'Invalid NPY signature.');
            version=fread(fid,2,'*uint8');
            assert(ismember(version(1),[1 2 3]),'Unsupported NPY version.');
            if version(1)==1, length=fread(fid,1,'uint16');
            else, length=fread(fid,1,'uint32'); end
            assert(~isempty(length) && length<=16*1024*1024,'Invalid NPY header length.');
            header=char(fread(fid,length,'*uint8')');
            token=regexp(header,'''descr''\s*:\s*''([^'']+)''','tokens','once');
            assert(~isempty(token),'Unsupported NPY dtype descriptor.'); dtype=token{1};
            token=regexp(header,'''shape''\s*:\s*\(([^)]*)\)','tokens','once');
            assert(~isempty(token),'Missing NPY shape.');
            dims=sscanf(strrep(token{1},',',' '),'%f')';
            assert(all(dims>=0 & dims==fix(dims)),'Invalid NPY dimensions.');
            count=prod(dims); fortran=contains(header,'''fortran_order'': True') || ...
                ~isempty(regexp(header,'''fortran_order''\s*:\s*True','once'));
            parsed=regexp(dtype,'^([<>=|])([biufcSU?])(\d+)$','tokens','once');
            assert(~isempty(parsed),'Unsupported/object NPY dtype: %s',dtype);
            endian=parsed{1}; kind=parsed{2}; width=str2double(parsed{3});
            byteWidth=width; if strcmp(kind,'U'), byteWidth=4*width; end
            offset=ftell(fid); fseek(fid,0,'eof'); available=ftell(fid)-offset; fseek(fid,offset,'bof');
            assert(count*byteWidth<=available,'Truncated NPY array.');
            machine='ieee-le'; if strcmp(endian,'>'), machine='ieee-be'; end
            if ismember(kind,{'S','U'})
                value=cell(count,1);
                for k=1:count
                    if strcmp(kind,'U')
                        codes=fread(fid,width,'*uint32',0,machine); codes=codes(codes~=0);
                        value{k}=GrimNpy.decodeUnicode(codes);
                    else
                        bytes=fread(fid,width,'*uint8')'; bytes=bytes(bytes~=0);
                        value{k}=native2unicode(bytes,'UTF-8');
                    end
                end
                if isempty(dims), value=value{1}; return; end
            else
                switch kind
                    case {'b','?'}, precision='uint8';
                    case 'i', precision=sprintf('int%d',8*width);
                    case 'u', precision=sprintf('uint%d',8*width);
                    case {'f','c'}
                        realWidth=width; if strcmp(kind,'c'), realWidth=width/2; end
                        assert(ismember(realWidth,[4 8]),'Unsupported floating-point precision.');
                        if realWidth==4, precision='single'; else, precision='double'; end
                end
                n=count; if strcmp(kind,'c'), n=2*n; end
                value=fread(fid,n,['*' precision],0,machine);
                assert(numel(value)==n,'Truncated NPY array.');
                if strcmp(kind,'c'), value=complex(value(1:2:end),value(2:2:end)); end
                if ismember(kind,{'b','?'}), value=logical(value); end
            end
            if isempty(dims), return; end
            if numel(dims)==1, value=reshape(value,1,dims(1));
            elseif fortran, value=reshape(value,dims);
            else, value=permute(reshape(value,fliplr(dims)),numel(dims):-1:1); end
        end

        function write(path,value,shape)
            if nargin<3, shape=[]; end
            text=ischar(value)||isstring(value)||iscellstr(value);
            if text
                if ischar(value), strings={value}; shape=[];
                else, strings=cellstr(value); shape=numel(strings); end
                codes=cellfun(@GrimNpy.encodeUnicode,strings,'UniformOutput',false);
                width=max([1;cellfun(@numel,codes(:))]); dtype=sprintf('<U%d',width);
            else
                if isempty(shape)
                    if isscalar(value), shape=[];
                    elseif isvector(value), shape=numel(value);
                    else, shape=size(value); end
                end
                switch class(value)
                    case 'double', dtype='<f8';
                    case 'single', dtype='<f4';
                    case 'logical', dtype='|b1';
                    otherwise
                        c=class(value); token=regexp(c,'^(u?int)(8|16|32|64)$','tokens','once');
                        assert(~isempty(token),'Unsupported metadata type: %s',c);
                        kind='i'; if startsWith(c,'u'), kind='u'; end
                        dtype=sprintf('<%s%d',kind,str2double(token{2})/8);
                end
                if ~isreal(value)
                    if isa(value,'single'), dtype='<c8'; else, dtype='<c16'; end
                end
            end
            tuple='';
            for k=1:numel(shape), tuple=[tuple sprintf('%d, ',shape(k))]; end %#ok<AGROW>
            header=sprintf('{''descr'': ''%s'', ''fortran_order'': True, ''shape'': (%s), }',dtype,tuple);
            padding=mod(64-mod(10+numel(header)+1,64),64);
            header=[header repmat(' ',1,padding) newline];
            assert(numel(header)<65536,'NPY header too long.');
            fid=fopen(path,'w','ieee-le'); assert(fid>=0,'Cannot create NPY file.');
            cleanup=onCleanup(@() fclose(fid)); %#ok<NASGU>
            fwrite(fid,uint8([147 double('NUMPY') 1 0]),'uint8');
            fwrite(fid,numel(header),'uint16'); fwrite(fid,uint8(header),'uint8');
            if text
                for k=1:numel(codes)
                    row=zeros(width,1,'uint32'); row(1:numel(codes{k}))=codes{k}; fwrite(fid,row,'uint32');
                end
            elseif ~isreal(value)
                pair=[real(value(:))';imag(value(:))']; fwrite(fid,pair,class(value));
            elseif islogical(value), fwrite(fid,value,'uint8');
            else, fwrite(fid,value,class(value)); end
        end

        function codes = encodeUnicode(text)
            units=uint32(char(text));
            if ~any(units>=55296 & units<=57343), codes=units; return; end
            codes=zeros(1,numel(units),'uint32'); n=0; k=1;
            while k<=numel(units)
                code=units(k);
                if code>=55296 && code<=56319 && k<numel(units)
                    next=units(k+1);
                    if next>=56320 && next<=57343
                        code=65536+(code-55296)*1024+(next-56320); k=k+1;
                    end
                end
                n=n+1; codes(n)=code; k=k+1;
            end
            codes=codes(1:n);
        end

        function text = decodeUnicode(codes)
            if all(codes<65536), text=char(uint16(codes(:)')); return; end
            units=zeros(1,2*numel(codes),'uint16'); n=0;
            for k=1:numel(codes)
                code=double(codes(k)); assert(code<=1114111,'Invalid Unicode codepoint.');
                if code>65535
                    code=code-65536; units(n+1)=55296+floor(code/1024);
                    units(n+2)=56320+mod(code,1024); n=n+2;
                else, n=n+1; units(n)=code; end
            end
            text=char(units(1:n));
        end
    end
end
