function [] = calc_eddy_velbot(runs)

    if isempty(runs.ubot) || isempty(runs.vbot)
        runs.read_velbot;
    end

    % get bottom velocity scale
    u = runs.ubot;
    v = runs.vbot;

    vel = sqrt(avg1(u(:,2:end-1,:), 1).^2 + ...
               avg1(v(2:end-1,:,:), 2).^2);
    %if size(vel,3) > length(runs.eddy.V)
    runs.eddy.Vb = squeeze(max(max(vel .* runs.eddy.vormask,[], 1), [], 2))';

    runs.eddy.hash = githash;

    eddy = runs.eddy;
    save([runs.dir '/eddytrack.mat'], 'eddy');

        %else
        %Ub = squeeze(max(max(bsxfun(@times, ...
        %                            vel(:,2:end-1,:), ...
        %                            runs.eddy.vormask), [], 1), [], 2));
        %end
end