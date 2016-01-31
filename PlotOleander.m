fname = '../data/All_Oleander_3D.mat';

ole = load(fname);

%% churchill (1986) section
% Oct 19-20, 1983
% stored in reverse order, that's why I reverse axis direction.
PlotOleanderSection(114, [11 17]);
export_fig -r150 -a2 images/oleander-oct1983.png

%contourf(max(ole.Dist(tt,:)) - repmat(ole.Dist(tt,:)', [1 500]), ...
%         -1*squeeze(ole.Depth(tt,:,:)), squeeze(ole.Temp_Fix(tt,:,:)), 40);
%colorbar;
%xlim([nanmin(ole.Dist(:)) nanmax(ole.Dist(:))/8]);

%%
tt = 480;
for tt=480:500
    cla
    contourf(repmat(ole.Dist(tt,:)', [1 500]), ...
             -1*squeeze(ole.Depth(tt,:,:)), squeeze(ole.Temp_Fix(tt,:,:)), 40);
    colorbar;
    xlim([nanmin(ole.Dist(:)) nanmax(ole.Dist(:))/8]);
    pause
end

%% Look at all XBT drops
NumTotal = 0;
NumFiltered = 0;
IndFiltered = [];
debugFlag = 0;

warning off signal:findpeaks:largeMinPeakHeight

ticstart = tic;
for tt=1:size(ole.Temp_Fix, 1)
    for xx=1:size(ole.Temp_Fix, 2)
        % disp([tt xx]);
        if all(isnan(ole.Depth(tt,xx,:))), break; end
        if max(ole.Depth(tt,xx,:)) < 50, continue; end

        T = cut_nan(squeeze(ole.Temp_Fix(tt,xx,:)));
        if ole.Year(tt,1) > 2008
            dT = smooth(diff(T), 5);
        else
            dT = smooth(diff(T), 2);
        end
        if length(dT) < 3, continue; end
        NumTotal = NumTotal + 1;

        [~,Peaks] = findpeaks(dT, 'MinPeakHeight', 0.15*max(abs(dT)));
        [~,Troughs] = findpeaks(-1*dT,  'MinPeakHeight', 0.15*max(abs(dT)));

        Peaks = sort([Peaks; Troughs]);

        if debugFlag
            plot(dT); linex(Peaks);
        end

        for pp=1:length(Peaks)-1
            if (dT(Peaks(pp)) < 0 & dT(Peaks(pp+1)) > 0) ...
                    & ( (Peaks(pp+1) - Peaks(pp)) > 3)
                if debugFlag
                    linex(Peaks(pp), [], 'r');
                    % keyboard;
                end
                IndFiltered(NumFiltered+1,:) = [tt xx];
                NumFiltered = NumFiltered+1;
            end
        end
    end
end
toc(ticstart);

%% save images
failed = [];
ticstart = tic;
hfig = figure; maximize;
for tt=unique(IndFiltered(:,1))'
    inds = find(IndFiltered(:,1) == tt);
    try
        PlotOleanderSection(tt, IndFiltered(inds,2)', hfig);
    catch ME
        failed = [failed tt];
    end
    export_fig('-r150','-a2', ['images/oleander/filtered/' ...
                        datestr([ole.Year(tt,1) ole.Month(tt,1) nanmin(ole.Day(tt,:)) 0 0 0], ...
                                'yyyy-mm-dd') '.png']);
end
disp('failed = ');
disp(failed);
toc(ticstart);