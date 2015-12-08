% describe fluxes
if ~exist('ew34', 'var') | ~strcmpi(ew.name, 'ew-34')
    ew34 = runs('../topoeddy/runew-34/');
end
factor = 2;
isobath = 4;

%% flux diagnostics
handles = ew34.plot_fluxts(factor, isobath, isobath);
handles.htitle.String = ['Flux of water above ' num2str(factor) 'H_{sb}'];
maximize; drawnow;
export_fig images/paper2/flux-diags.png

%% x-z sections
handles = ew34.plot_xzsection(isobath, 225);
correct_ticks('y', [], 3, handles.hax([1 3 4]));
correct_ticks('y', [], 4, handles.hax([2]));
drawnow;
export_fig images/paper2/ew-34-xzsection.png