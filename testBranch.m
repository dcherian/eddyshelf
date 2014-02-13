% function that helps me test branches

function [] = testBranch(indices, sz, num)
    [ix,iy] = ind2sub(sz,indices);
    
    %figure;
    plot(ix,iy,'*-');
    hold on
    plot(ix(1),iy(1),'ro','MarkerSize',16);
    plot(ix(end),iy(end),'ko','MarkerSize',16);
    
    if ~exist('num','var'), num = []; end
    text(ix(1),iy(1),num2str(num));
    title(['red = start | black = end | isclockwise = ' ...
        num2str(isclockwise(ix,iy))]);
    
end