function [] = CalcNonShelfWaterCensus(runs)

    if ~isfield(runs, 'water')
        return;
    end
    runs.water.nonsh.shelf = runs.water.shvol ...
        - runs.water.sh.shelf;
end