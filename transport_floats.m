dir = ['runs/topoeddy/runteb-04-hires-7/'];
run = runs(dir);

%% masks 

% first bin floats?
xvec = run.rgrid.xr(:,1);
yvec = run.rgrid.yr(1,:)';
sz = size(run.rgrid.xr);
floats = run.roms;
eddy = run.eddy;
bathy = run.bathy;

shelfmask = ones([1 size(floats.x,2)]);
shelfmask(floats.init(:,1) < bathy.xsb) = NaN;

for i=1:size(eddy.mask,3)
    fltloc = zeros(sz);
    ff = ceil(i*floats.fac);
    if ff > size(floats.x,1), break; end
    % determine grid indices & make float locations a matrix
    xg = vecfind( xvec, cut_nan(floats.x(ff,:) .* shelfmask) );
    yg = vecfind( yvec, cut_nan(floats.y(ff,:) .* shelfmask) );
    fltloc(sub2ind(sz,xg,yg)) = 1;
    
    % remove floats inside eddy contour and those that started offshore of
    % shelbreak
    netmask = ~eddy.mask(:,:,i);
    
    % truncate because of eddy.mask
    fltloc = fltloc(2:end-1,2:end-1);    
    
    clf;
    [cc,hh] = contour(eddy.xr/1000,eddy.yr/1000,eddy.mask(:,:,i),1,'w');
    hold on
    pcolorcen(eddy.xr/1000,eddy.yr/1000,fltloc .*netmask);
    colormap(gray);    
    set(hh,'LineWidth',2);
    if run.bathy.axis == 'y'
        liney(bathy.xsb/1000,'shelfbreak','w');
    else
        linex(bathy.xsb/1000,'shelfbreak','w');
    end
    %plot(floats.x(i,:)/1000,floats.y(i,:)/1000,'k.');
    title(num2str(i));
    pause(0.01);
end
