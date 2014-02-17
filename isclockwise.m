function [out] = isclockwise(inx,iny)
    % based on http://stackoverflow.com/questions/14505565/detect-if-a-set-of-points-in-an-array-that-are-the-vertices-of-a-complex-polygon
    area = sum(inx(1:end-1) .* iny(2:end) - inx(2:end) .* iny(1:end-1));
    
    area = area + inx(end).*iny(1) - inx(1) .* iny(end);
    
    if area < 0, out = true; else out = false; end
%     
%     area2 = 0;
%     for ii=1:length(inx)
%         jj = ii + 1;
%         if jj > length(inx), jj = 1; end
%         area2 = area2 + inx(ii) * iny(jj) - inx(jj)*iny(ii);
%     end
%     area2
%     area
end
