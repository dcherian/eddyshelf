function [S] = bathy2_x(S,B,X,Y)

    ix1 = find_approx(S.x_rho(:,1),X-B.L_entry-B.L_tilt,1);
    ix2 = find_approx(S.x_rho(:,1),X-B.L_entry,1);
    iy1 = find_approx(S.y_rho(1,:),Y-B.L_shelf,1);
    iy2 = find_approx(S.y_rho(1,:),Y-B.L_shelf2,1);
    
    x1 = S.x_rho(ix1,1); x2 = S.x_rho(ix2,1);
    y1 = S.y_rho(1,iy1); y2 = S.y_rho(1,iy2);
    
    dx = mean(diff(S.x_rho(:,1)));
    dy = mean(diff(S.y_rho(1,:)));
    
    % equation for tilted line : y=mx+c
    m = (iy2-iy1)./(ix2-ix1);
    c = (iy1*ix2 - iy2*ix1)./(ix2-ix1);
    
    % number of points for slope - t=tilt
    nsly = ceil(B.L_slope./dy);
    theta = atan(m);
    nsltx = floor(nsly * sin(theta));
    nslty = ceil(nsly * cos(theta));
    
    % equation for displaced tilted line - for slope
    c2 = ( (iy1-nsly) *ix2 - (iy2-nsly)*ix1)./(ix2-ix1);
    
    % build shelf mask
    mask_sh = zeros(size(S.x_rho));
    mask_sh(1:ix1,iy1:end) = 1; % left region
    [xt,yt] = meshgrid(ix1:ix2,iy1:iy2); % tilt
    xt = xt'; yt = yt';
    mask_sh(ix1:ix2,iy1:iy2) = ((yt - m*xt - c) > 0);
    mask_sh(1:end,iy2:end) = 1; % right region
    
    imagesc(mask_sh');set(gca,'ydir','normal');

    % build slope mask
    mask_sl = zeros(size(S.x_rho));
    [xt,yt] = meshgrid(ix1:ix2,iy1-nsly:iy2-nsly); % tilt
    xt = xt'; yt = yt';
    mask_sl(ix1:ix2,iy1-nsly:iy2-nsly) = ((yt - m*xt - c2) > 0);
    mask_sl(:,1:end-nsly) = mask_sh(:,nsly+1:end);
    mask_sl(:,iy2:end) = 1;
    
    mask_sh = logical(mask_sh);
    mask_sl = logical(mask_sl);
    
    mask_sl = xor(mask_sh,mask_sl);
    imagesc(mask_sl');set(gca,'ydir','normal');
    
    %% make bathymetry
    S.h = zeros(size(S.h));
    S.h(mask_sh) = B.H_shelf + B.sl_shelf .* (Y - S.y_rho(mask_sh));
    
    mask_west = mask_sl & (S.x_rho < x1);
    mask_east = mask_sl & (S.x_rho > x2);
    mask_tilt = mask_sl & ~(mask_east | mask_west);
    
    hmax  = repmat(nanmax(S.h,[],2),1,size(S.h,2));
    heast = nanmax(S.h(S.x_rho > x2));
    
    S.h = S.h + mask_west.*(hmax + B.sl_slope .* (y1-S.y_rho));
    
    hdeep = nanmax(S.h(:));
    
    sl_east = (hdeep-heast)./(nsly)/dy;
    S.h = S.h + mask_east.*(heast + sl_east .* (y2-S.y_rho));
    
    hmaxt = repmat(nanmax(S.h .* ~mask_east,[],1),size(S.h,1),1);
    ntilt = repmat(sum(mask_tilt,1),size(S.h,1),1);
    sl_tilt = (hdeep-hmaxt)./ntilt./dx;
    xtilt_left = repnan(nanmin(S.x_rho.*fillnan(double(mask_tilt),0),[],1),0);
    xtilt_left = repmat(xtilt_left,size(S.h,1),1);
    S.h = S.h + repnan(mask_tilt.*(hmaxt + sl_tilt .* ((S.x_rho - xtilt_left))),0);
    
    S.h = repnan(fillnan(S.h,0),hdeep);

    imagescnan(-S.h'); colorbar
    