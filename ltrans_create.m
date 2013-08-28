function [] = ltrans_create(rgrid,zeta,eddy)
    fname = 'runs/ltrans_init.txt';
    % open file
    fid = fopen(fname,'w');
    if fid == -1
        error('Couldn''t open file to write.');
    end
    
    contourf(rgrid.x_rho/1000,rgrid.y_rho/1000,zeta(:,:,1)');
    shading flat; axis image
    caxis([min(zeta(:)) max(zeta(:))]);
    hold on
    [C,hc] = contour(rgrid.h,[114 500 750 1100],'k');
    clabel(C,hc);
    if exist('eddy','var') && ~isempty(eddy)
       plot(eddy.mx/1000,eddy.my/1000,'k','LineWidth',2);
    end
    
    % create distribution
    xr = rgrid.x_rho';    yr = rgrid.y_rho'; zr = permute(rgrid.z_r,[3 2 1]);
    flag = input('From image? (0/1) : ');
    if flag
        str = 'select rectangle for float deployment';
        title(str);
        disp(str);
        [floatx,floaty] = select_rect();
        xlo = find_approx(xr(:,1),min(floatx(:))*1000,1); 
        xhi = find_approx(xr(:,1),max(floatx(:))*1000,1); 
        ylo = find_approx(yr(1,:),min(floaty(:))*1000,1); 
        yhi = find_approx(yr(1,:),max(floaty(:))*1000,1); 
    else
        dd = input('[xlo,xhi,ylo,yhi] in index: ');
       xlo = dd(1); xhi = dd(2);
       ylo = dd(3); yhi = dd(4);
    end
    xlo = floor(xlo); ylo = floor(ylo);
    xhi = floor(xhi); yhi = floor(yhi);
    dd = input('[zlo,zhi] in m: ');
    zlo = dd(1); zhi = dd(2);
    
    flag = 0;
    while flag == 0
        fprintf('\n\n ( %d:%d, %d:%d, %d:%d ) | ( %.2f:%.2f, %.2f:%.2f )\n\n', ...
            xlo,xhi,ylo,yhi,zlo,zhi,xr(xlo,1)/1000,xr(xhi,1)/1000, ...
            yr(1,ylo)/1000,yr(1,yhi)/1000);
        dd = input('enter [dx,dy,dz] in pts :');
        dx = dd(1); dy = dd(2); dz = dd(3); clear dd

        dd = input('enter [tlo,thi,dt] in days: ');
        tlo = dd(1); thi = dd(2); dt = dd(3); clear dd;

        xvec = xlo:dx:xhi; yvec = ylo:dy:yhi; zvec =zlo:dz:zhi;
        tvec = tlo:dt:thi;
        X = length(xvec); Y= length(yvec); Z = length(zvec); T = length(tvec);
        fprintf('\n Total Number of floats = %d | Number released at once = %d\n\n', ...
                    X*Y*Z*T, X*Y*Z);
        flag = input('\nContinue (0/1)?');
    end
    
    [p,q,r,s] = ndgrid(xvec,yvec,zvec,tvec);
    ind = [p(:) q(:) r(:)];
    zz = zr(xvec,yvec,zvec);
    zz = repmat(zz,[1 1 1 T]);
    
    init = [xr(p(:),1) yr(1,q(:))' zz(:) s(:)*86400];
    
    
    % write to file
    fprintf(fid,'%.2f,%.2f,%.2f,%d \n',init');
    fprintf('Added %d floats to %s \n\n\n',size(init,1),fname);

    fclose(fid);
    
    disp('ROMS Parameters');
    fprintf(['\n Ft0 = %d | Fx0 = %d | Fy0 = %d | Fz0 = %d' ...
             '\n Fdt = %d | Fdx = %d | Fdy = %d | Fdz = %d\n\n'],tlo,xlo,ylo,zlo, ...
             dt,dx,dy,dz);
end