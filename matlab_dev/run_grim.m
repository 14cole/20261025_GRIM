function app = run_grim()
%RUN_GRIM Start the self-contained MATLAB GRIM plotting tool (R2024b).
%   app = run_grim;
% All runtime dependencies are .m files in this folder. No Python is used.
    addpath(fileparts(mfilename('fullpath')));
    app=GrimPlotApp;
end
