function tests = test_grim
% Run from matlab_dev: results = runtests('tests'); assertSuccess(results)
% Fixed Python-reference fixtures are data only; these tests use MATLAB only.
    tests=functiontests(localfunctions);
end

function setupOnce(tc)
    folder=fileparts(mfilename('fullpath')); addpath(fileparts(folder));
    tc.TestData.Fixtures=fullfile(folder,'fixtures');
    tc.TestData.Reference=load(fullfile(folder,'fixtures','python_reference.mat'));
    tc.TestData.A=GrimIO.load(fullfile(folder,'fixtures','base.grim'));
    tc.TestData.B=GrimIO.load(fullfile(folder,'fixtures','other.grim'));
end

function compare(tc,actual,expected,tolerance)
    if nargin<4, tolerance=2e-6; end
    verifyEqual(tc,actual.Azimuth(:),expected.Azimuth(:),'AbsTol',1e-6);
    verifyEqual(tc,actual.Elevation(:),expected.Elevation(:),'AbsTol',1e-6);
    verifyEqual(tc,actual.Frequency(:),expected.Frequency(:),'AbsTol',1e-6);
    verifyEqual(tc,strjoin(actual.Polarization,','),strtrim(expected.Polarization));
    verifyEqual(tc,actual.Quantity,strtrim(expected.Quantity));
    verifyEqual(tc,size(actual.Power),size(expected.Power));
    verifyEqual(tc,isnan(actual.Power),isnan(expected.Power));
    valid=isfinite(actual.Power)&isfinite(expected.Power);
    verifyEqual(tc,actual.Power(valid),expected.Power(valid),'RelTol',tolerance,'AbsTol',1e-11);
    verifyEqual(tc,isnan(actual.Phase),isnan(expected.Phase));
    valid=isfinite(actual.Phase)&isfinite(expected.Phase)&expected.Power>0;
    difference=abs(angle(exp(1i*(actual.Phase(valid)-expected.Phase(valid)))));
    verifyLessThanOrEqual(tc,difference,tolerance*ones(size(difference)));
end

function testNativeInputs(tc)
    compare(tc,tc.TestData.A,tc.TestData.Reference.base);
    compare(tc,tc.TestData.B,tc.TestData.Reference.other);
    for key={'raw_2d','power_only'}
        d=GrimIO.load(fullfile(tc.TestData.Fixtures,[key{1} '.grim']));
        compare(tc,d,tc.TestData.Reference.(key{1}));
    end
    raw=GrimIO.load(fullfile(tc.TestData.Fixtures,'raw_2d.grim'));
    verifyEqual(tc,GrimData.db(raw),tc.TestData.Reference.raw_2d_db,'AbsTol',1e-10);
end

function testVendorReaders(tc)
    files=tc.TestData.Reference.file_names;
    for k=1:numel(files)
        name=strtrim(files{k}); d=GrimIO.load(fullfile(tc.TestData.Fixtures,name));
        key=['file_' strrep(name,'.','_')]; compare(tc,d,tc.TestData.Reference.(key));
    end
end

function testArithmeticAgainstPython(tc)
    ops={'coherent_add','coherent_subtract','coherent_divide','incoherent_add','incoherent_subtract','db_difference'};
    for k=1:numel(ops)
        d=GrimData.combine({tc.TestData.A,tc.TestData.B},ops{k});
        compare(tc,d,tc.TestData.Reference.(ops{k}),1e-10);
    end
end

function testTransformsAgainstPython(tc)
    specs={'flip_phase',[];'offset',3;'mirror',7;'swap_el_az',[];'shift',[3 -2 45];'extrusion',2};
    for k=1:size(specs,1)
        d=GrimData.transform(tc.TestData.A,specs{k,1},specs{k,2});
        compare(tc,d,tc.TestData.Reference.(specs{k,1}),1e-10);
    end
    d=GrimData.regrid(tc.TestData.A,'Azimuth',[-30 -15 0 15 30]);
    compare(tc,d,tc.TestData.Reference.regrid,1e-10);
    d=GrimData.medianize(tc.TestData.A,30,15); compare(tc,d,tc.TestData.Reference.medianize,1e-10);
end

function testStatisticsAgainstPython(tc)
    compare(tc,GrimData.reduce(tc.TestData.A,[1 3],'mean'),tc.TestData.Reference.mean,1e-10);
    compare(tc,GrimData.reduce(tc.TestData.A,1,'percentile',90,false),tc.TestData.Reference.percentile,1e-10);
    compare(tc,GrimData.reduce(tc.TestData.A,1,'percentile',90,true),tc.TestData.Reference.percentile_linear,1e-10);
end

function testAgreementAgainstPython(tc)
    ref=tc.TestData.Reference.agreement; result=GrimPlots.agreement(ref.Left(:),ref.Right(:));
    verifyEqual(tc,result.Score,ref.Score,'AbsTol',1e-10);
    verifyEqual(tc,result.MAE,ref.MAE,'AbsTol',1e-10);
end

function testInterchangeRoundTrips(tc)
    folder=tempname; mkdir(folder); cleanup=onCleanup(@() rmdir(folder,'s')); %#ok<NASGU>
    for ext={'.grim','.mat','.csv','.tsv'}
        path=fullfile(folder,['roundtrip' ext{1}]); GrimIO.save(path,tc.TestData.A);
        compare(tc,GrimIO.load(path),tc.TestData.Reference.base,1e-8);
    end
    d=GrimIO.load(fullfile(tc.TestData.Fixtures,'cut.pio'));
    path=fullfile(folder,'roundtrip.pio'); GrimIO.save(path,d);
    compare(tc,GrimIO.load(path),tc.TestData.Reference.file_cut_pio,1e-8);
    % PTM needs a uniformly spaced non-wrapping cut for this writer fixture.
    d=GrimData.create([-10 0 10],0,linspace(8,12,40),{'HH'},ones(3,1,40),zeros(3,1,40));
    path=fullfile(folder,'roundtrip.ptm'); GrimIO.save(path,d); out=GrimIO.load(path);
    verifyEqual(tc,out.Power,d.Power,'AbsTol',1e-6); verifyEqual(tc,out.Azimuth,d.Azimuth);
end

function testMissingPhaseAndInvalidOperations(tc)
    d=tc.TestData.A; d.Phase(:)=NaN;
    verifyTrue(tc,throws(@() GrimData.field(d)));
    verifyTrue(tc,throws(@() GrimData.combine({tc.TestData.B,tc.TestData.A},'incoherent_subtract')));
    verifyTrue(tc,throws(@() GrimData.combine({tc.TestData.A},'coherent_add')));
    d=tc.TestData.A; d.Azimuth=[0 .01 .02];
    verifyTrue(tc,throws(@() GrimData.transform(d,'round',[0 3 3])));
    verifyEqual(tc,GrimData.roundEven([-2.5 -1.5 .5 1.5 2.5],0),[-2 -2 0 2 2]);
end

function testJoinCropAndOverlap(tc)
    d=tc.TestData.A;
    a=GrimData.slice(d,{1:2,1:2,1:3,1:2}); b=GrimData.slice(d,{2:3,1:2,1:3,1:2});
    compare(tc,GrimData.join({a,b}),tc.TestData.Reference.base);
    common=GrimData.overlap({a,b});
    verifyEqual(tc,common{1}.Azimuth,0); verifyEqual(tc,common{1}.Power,common{2}.Power);
    b.Power(1)=b.Power(1)+1; verifyTrue(tc,throws(@() GrimData.join({a,b})));
    partial=GrimData.regrid(d,'Azimuth',[-40 -30 0 30 40]);
    verifyTrue(tc,all(isnan(partial.Power([1 end],:,:,:)),'all'));
end

function testSafeCoordinateParser(tc)
    verifyEqual(tc,GrimBackend.parseNumbers('-2:0.5:0'),[-2 -1.5 -1 -.5 0]);
    verifyEqual(tc,GrimBackend.parseNumbers('[3, 1, 2]'),[1 2 3]);
    verifyTrue(tc,throws(@() GrimBackend.parseNumbers('system("anything")')));
    verifyTrue(tc,throws(@() GrimBackend.parseNumbers('1:0:5')));
end

function testGuiCallbacksAndPlotModes(tc)
    % Desktop-only UI smoke coverage. Numeric tests work in MATLAB -batch.
    assumeTrue(tc,usejava('desktop'),'GUI smoke test requires the MATLAB desktop.');
    app=GrimPlotApp; cleanup=onCleanup(@() delete(app)); %#ok<NASGU>
    app.UIFigure.Visible='off'; app.Backend.addDatasets({tc.TestData.A,tc.TestData.B});
    modes={'azimuth_rect','azimuth_polar','frequency','elevation_sweep','waterfall','delta_map','compare','sector_stats'};
    for k=1:numel(modes)
        mode=modes{k}; rows=1;
        if ismember(mode,{'delta_map','compare'}), rows=[1 2]; end
        app.selectDatasets(rows); app.Controls.PlotType.Value=mode;
        app.ParameterLists.Frequency.Value=1;
        if ismember(mode,{'waterfall','delta_map','frequency'}), app.ParameterLists.Frequency.Value=1:3; end
        state=app.getState();
        GrimPlots.render(app,app.Backend.Datasets(rows),state);
        verifyNotEmpty(tc,app.getAxes().Children);
    end
    app.Backend.handle('duplicate',app.getState(),struct());
    verifyEqual(tc,numel(app.Backend.Datasets),3);
end

function result = throws(callback)
    result=false; try, callback(); catch, result=true; end
end
