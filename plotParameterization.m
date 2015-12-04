function [] = plotParameterization(name)

    load(['./params/param_' name]);

    figure;
    insertAnnotation(['plotParameterization | use_hash = ' hash(1:15)]);
    slopeplot = 3;
    hax = packboth(3,3);

    for ii=1:length(slope)
        if ii >= slopeplot, hh = ii+1; else hh = ii; end
        axes(hax(hh));

        norm(ii,:) = 1; %norm(ii,:) / 1000;
        for rr=1:length(norm(ii,:))
            paramvalue = slope(ii) * plotx(ii,rr) + intercept(ii);
            herrhi = (slope(ii)+merr(ii)) * plotx(ii,rr) + (intercept(ii) + cerr(ii)) - ...
                     paramvalue;
            herrlo = abs((slope(ii)-merr(ii)) * plotx(ii,rr) + (intercept(ii) - cerr(ii)) - ...
                         paramvalue);

            plot(paramvalue./norm(ii,rr), diags(ii,rr)./norm(ii,rr), '.', ...
                 'Color', color(rr,:));
            hold on;
            h = errbar(paramvalue./norm(ii,rr), diags(ii,rr)./norm(ii,rr), ...
                       herrlo./norm(ii,rr), herrhi./norm(ii,rr), '-', 'horiz');
            h.Color = color(rr,:);

            if ~isempty(isnan(err))
                h = errbar(paramvalue./norm(ii,rr), diags(ii,rr)./norm(ii,rr), ...
                           err(ii,rr)./norm(ii,rr), '-');
                h.Color = color(rr,:);
            end
        end

        xlabel(''); ylabel(''); title('');
        beautify([14 16 18]); pbaspect([1.732 1 1]);
        ggplot;
        hax(hh).XColor = [1 1 1];
        hax(hh).YColor = [1 1 1];

        ylim([0 max(ylim)]);
        xlim([0 max(xlim)]);
        line45;
    end


    hax(1).YColor = [1 1 1]*0.3;
    hax(4).YColor = [1 1 1]*0.3;
    hax(7).YColor = [1 1 1]*0.3;

    hax(7).XColor = [1 1 1]*0.3;
    hax(8).XColor = [1 1 1]*0.3;
    hax(9).XColor = [1 1 1]*0.3;


    ylabel(['Diagnosed ' name '/ Eddy volume flux']);
    xlabel(['Prediction / Eddy volume flux']);
end