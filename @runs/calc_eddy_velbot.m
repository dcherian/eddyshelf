function [] = calc_eddy_velbot(runs)

    if isempty(runs.ubot) || isempty(runs.vbot)
        runs.read_velbot;
    end

    % get bottom velocity scale
    u = runs.ubot;
    v = runs.vbot;
    %if size(vel,3) > length(runs.eddy.V)

    debug = 0;

    vel = sqrt(avg1(u(:,2:end-1,:), 1).^2 + ...
               avg1(v(2:end-1,:,:), 2).^2);

    %vel = avg1(u(:,2:end-1,:),1);
    runs.eddy.Vb = squeeze(mean(mean(vel .* runs.eddy.mask, 1), ...
                               2))';

    if debug
        eddbot = dc_roms_read_data(runs.dir, runs.eddname, [], {'z' 1 ...
                            1}, [], runs.rgrid, 'his', 'single') > runs.eddy_thresh;
        rhobot = dc_roms_read_data(runs.dir, runs.eddname, [], {'z' 1 ...
                            1}, [], runs.rgrid, 'his', 'single') > runs.eddy_thresh;

        Vb2 = squeeze(max(max(vel .* runs.eddy.mask .* eddbot(2:end-1,2:end-1,:),[], 1), ...
                          [], 2))';

        Vb3 = squeeze(max(max(vel .* runs.eddy.mask .* rhobot(2:end-1,2:end-1,:),[], 1), ...
                          [], 2))';

        figure;
        hold all;
        plot(runs.eddy.Vb);
        plot(Vb2);
        plot(Vb3);
        legend('eddy mask', 'eddy dye', 'rho mask');
    end

    runs.eddy.hash = githash([mfilename('fullpath') '.m']);

    eddy = runs.eddy;
    save([runs.dir '/eddytrack.mat'], 'eddy');
end