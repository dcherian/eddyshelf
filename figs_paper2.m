% describe fluxes
if ~exist('ew', 'var') | ~strcmpi(ew.name, 'ew-34')
    ew = runs('../topoeddy/runew-34/');
end

factor = 2;
handles = ew.plot_fluxts(factor, 3, 3);
handles.htitle.String = ['Flux of water above ' num2str(factor) 'H_{sb}'];
maximize; drawnow;
export_fig images/paper2/flux-diags.png
