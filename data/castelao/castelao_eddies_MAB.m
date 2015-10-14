%clear all
%close all
%pack

%cd /media/data/Work/eddyshelf/data/castelao/

% aviso data
load aviso_daily_1992_2011_7days_40km_gap_Leff_mod_MAB
% time snapshots
load aviso_data_1992_2011_eddies_cyclones_7days_40km_gap_Leff_mod_MAB
load aviso_data_1992_2011_eddies_anticyclones_7days_40km_gap_Leff_mod_MAB
% trajectory files
load 03eddy_trajectories_7days_gap40km_1992_Leff_mod_MAB_v5
load 03eddy_trajectories_combined_7days_gap40km_1992_Leff_mod_MAB_v5
load 03eddy_trajectories_7days_gap40km_1992_Leff_mod_MAB_v5_frequency
load topo_useast
   
cx=coast(:,1);
cy=coast(:,2);
XX=Xuseast;
YY=Yuseast;
ZZ=Zuseast;

%% renato's animation
close all
 
for BB=1:size(ssh30,3)
    [cs,kk]=contourf(xx,yy,ssh30(:,:,BB),50);
    set(kk,'linestyle','none')
    colorbar
    hold on
    plot(cx,cy,'k')
    caxis([-40 40])
    
    boundary=cBOUNDARY{BB};
    for kkk=1:length(boundary)
        plot(xx(boundary{kkk}),yy(boundary{kkk}),'.r')
    end
    plot(cXXmean{BB},cYYmean{BB},'or','markerfacecolor','r')
%    plot(cXXext{BB},cYYext{BB},'or','markerfacecolor','r')
    
    boundary=aBOUNDARY{BB};
    for kkk=1:length(boundary)
        plot(xx(boundary{kkk}),yy(boundary{kkk}),'.b')
    end
    plot(aXXmean{BB},aYYmean{BB},'ob','markerfacecolor','b')
%    plot(aXXext{BB},aYYext{BB},'ob','markerfacecolor','b')
    
    title([num2str(yearssh(BB)) '    ' num2str(dayssh(BB))])
    pause(0.1)
    clf
end

%%
close all
 
for BB=900:size(ssh30,3)
    [cs,kk]=contourf(xx,yy,ssh30(:,:,BB),50);
    set(kk,'linestyle','none')
    colorbar
    hold on
    % coastline
    plot(cx,cy,'k')
    % bathymetry
    [cc,hh] = contour(XX, YY, ZZ, [500 1000 2000 3000 4000 5000 6000],'k');
    %xclabel(cc,hh);
    caxis([-40 40])
    
    boundary=cBOUNDARY{BB};
    for kkk=1:length(boundary)
        plot(xx(boundary{kkk}),yy(boundary{kkk}),'.r')
    end
    plot(cXXmean{BB},cYYmean{BB},'or','markerfacecolor','r')
%    plot(cXXext{BB},cYYext{BB},'or','markerfacecolor','r')
    
    boundary=aBOUNDARY{BB};
    for kkk=1:length(boundary)
        plot(xx(boundary{kkk}),yy(boundary{kkk}),'.b')
    end
    plot(aXXmean{BB},aYYmean{BB},'ob','markerfacecolor','b')
%    plot(aXXext{BB},aYYext{BB},'ob','markerfacecolor','b')
    
    title([num2str(yearssh(BB)) '    ' num2str(dayssh(BB))])
    pause(0.1)
    clf
end

%% example of trajectory of an anticyclone
close all

mm = 5;
aaaa = 576
te=find(yearssh==eayear{aaaa}(1) & dayssh== eaday{aaaa}(1)) ;  
[cc,hh] = contour(Xuseast(1:mm:end,1:mm:end), ...
                      Yuseast(1:mm:end,1:mm:end), ...
                      Zuseast(1:mm:end,1:mm:end), ...
                      [100 200 500 1000 2000 3000 4000],'k');
                  hold on
clabel(cc,hh);

ap=eaX{aaaa};

j =1;
hssh = pcolor(xx,yy,ssh30(:,:,te-1+j));
shading interp
caxis([-30 30])
axis image;
hold on
plot(cx,cy,'k')
ylim([31 43]);

for j=2:length(ap)
    set(hssh,'CData',ssh30(:,:,te-1+j));
    plot(eaX{aaaa}(1:j),eaY{aaaa}(1:j),'ob-','markerfacecolor','b')
    title(j)
    pause
end

%% example of trajectory of a cyclone

close all

aaaa=202;
    
te=find(yearssh==ecyear{aaaa}(1) & dayssh== ecday{aaaa}(1))   
   
ap=ecX{aaaa};
for j=1:length(ap)
pcolor(xx,yy,ssh30(:,:,te-1+j))
shading interp
caxis([-30 30])
hold on
plot(cx,cy,'k')
plot(ecX{aaaa}(1:j),ecY{aaaa}(1:j),'or','markerfacecolor','r')
plot(ecX{aaaa}(1:j),ecY{aaaa}(1:j),'r')
title(j)
pause
clf
end


%% ploting frequency


xlim=[min(min(xx)) max(max(xx))];
ylim=[min(min(yy)) max(max(yy))];
xtam=diff(xlim)*cos(mean(ylim)*pi/180);
ytam=diff(ylim);

clim=[0 40];
figure
    subplot('position',[0.100    0.438    0.3847    0.3812])
    pcolor(xx,yy,cyclones);
    shading interp
    hold on
    plot(cx,cy,'k')
    caxis(clim)
    contour(XX,YY,ZZ,[60 200 800],'color',[0.7 0.7 0.7])
    axis([xlim ylim])
    set(gca,'plotboxaspectratio',[1 ytam/xtam 1])
    set(gca,'box','on')
    text(-81.4,42,'cyclonic')

    
    subplot('position',[0.500    0.438    0.3847    0.3812]) 
    pcolor(xx,yy,anticyclones);
    shading interp
    hold on
    plot(cx,cy,'k')
    caxis(clim)
    contour(XX,YY,ZZ,[60 200 800],'color',[0.7 0.7 0.7])
    axis([xlim ylim])
    set(gca,'plotboxaspectratio',[1 ytam/xtam 1])
    set(gca,'box','on')
    text(-81.4,42,'anticyclonic')
    set(gca,'yticklabel',[])
    po=get(gca,'position')
    set(gca,'position',[po(1) po(2) po(3) po(4)]);
    po=get(gca,'position')
    cb=colorbar;
    po1=get(cb,'position');
    set(gca,'position',po)
    set(cb,'position',[po1(1)+0.09 po1(2) po1(3)*0.5 po1(4)])
    