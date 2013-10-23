load sponge % get xr,yr

width = 30*1000;
maxvisc = 50;
bdiff = 1;

XX = max(xr(:)); YY = max(yr(:));

visc2 = ones(size(xr));

for j=1:size(xr,2)
    for i=1:size(xr,1)
        if yr(i,j) > YY-width
           visc2(i,j) =  visc2(i,j) + (yr(i,j)+width-YY)*maxvisc/width;
        end
    end
end
for j=1:size(xr,2)
    for i=1:size(xr,1)
        if xr(i,j) > XX-width
           visc2(i,j) =  visc2(i,j) + (xr(i,j)+width-XX)*maxvisc/width;
        end
    end
end
visc2(visc2>maxvisc) = maxvisc;
pcolor(xr/1000,yr/1000,visc2)
set(gca,'ydir','normal');
axis image;shading flat;colorbar