function [] = ltrans_create(rgrid,zeta)
    fname = 'runs/ltrans_init.txt';
    % open file
    fid = fopen(fname,'w');
    if fid == -1
        error('Couldn''t open file to write.');
    end
    
    contourf(zeta(:,:,1)');
    shading flat; axis image
    caxis([min(zeta(:)) max(zeta(:))]);
    hold on
    [C,hc] = contour(rgrid.h,[114 500 750 1100],'k');
    clabel(C,hc);
    
    % create distribution
    xr = rgrid.x_rho';    yr = rgrid.y_rho'; zr = permute(rgrid.z_r,[3 2 1]);
    flag = input('From image? (0/1) : ');
    if flag
        title('xlo'); [xlo,~] = ginput(1);
        title('xhi'); [xhi,~] = ginput(1);
        title('ylo'); [~,ylo] = ginput(1);
        title('yhi'); [~,yhi] = ginput(1);
    else
       xlo = 5; xhi = 34;
       ylo = 30; yhi = 128;
    end
    xlo = floor(xlo); ylo = floor(ylo);
    xhi = floor(xhi); yhi = floor(yhi);
    zlo = 1; zhi = rgrid.N;
    figure(gcf); linex([xlo xhi]); liney([ylo yhi]);
    
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
        fprintf('\n Number of floats = %.0f | ',X*Y*Z*T);
        flag = input('Continue (0/1)?');
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