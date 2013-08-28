% applies mask to all fields in structure data

function [data] = apply_mask(data,mask)

    names = fieldnames(data);   

    for i=1:length(names)
        data.(char(names(i,:))) = cut_nan(data.(char(names(i,:))) .* fillnan(mask,0));
    end