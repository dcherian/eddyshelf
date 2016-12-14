function [] = read_diagnostics()

    fname = '''../topoeddy/runew-6452/ocean_dia.nc.new02''';

    var = 'vbar';
    terms = {'accel';'cor'; 'prsgrd';'xadv'; 'yadv'; ...
              'xvisc'; 'yvisc'};

    for ii=1:length(terms)
        vname = [var '_' terms{ii}];
        disp(vname);
        tic;
        eval([ vname ' = single(ncread(' fname ',''' vname '''));']);
        toc;
    end

    figure;
    hold all;
    for ii=1:length(terms)
        vname = [var '_' terms{ii}];
        eval(['plot(integrate(' vname '));']);
    end
    legend(terms);
    title(var);
end

function [out] = integrate(in)
    if ndims(in) > 3
        in = sum(in,3);
    end
    out = squeeze(sum(sum(in,1),2));
end