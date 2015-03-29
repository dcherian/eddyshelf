% I want to fit y = y_0 tanh((x-x0)/X) + y1 * (x-x0)/X + yc
function [y0,X,y1,x0,yc] = tanh_fit(x, y, plot_flag, test)

    if ~exist('test', 'var'), test = 0; end
    if ~exist('plot_flag', 'var'), plot_flag = 0; end

    if test
        test_fit();
        return;
    end

    mny = mean(y);
    y = y - mean(y);

    % The fit is very sensitive to initial guesses.
    % if getting crap results, play around with these values.
    initGuess(1) = y(1);
    initGuess(2) = mean(x(:))/3;
    initGuess(3) = 1e-3*initGuess(1);
    initGuess(4) = mean(x(:));
    initGuess(5) =  y(end);

    opts = optimset('MaxFunEvals',1e7, 'TolFun', 1e-11);
    [fit2,~,exitflag] = fminsearch(@(fit) fiterror(fit,x,y), ...
                                   initGuess,opts);

    y0 = fit2(1);
    X = fit2(2);
    y1 = fit2(3);
    x0 = fit2(4);
    yc = fit2(5) + mny;
    y = y + mny;

    if plot_flag
        figure;
        plot(x,y,'k*'); hold all
        plot(x, y0*tanh((x-x0)/X) + y1*((x-x0)/X) + yc);
        %linex([1 2 3]*X);
    end
end

function [E] = fiterror(fit,x,y)
% x = (T0,H,a)
    y0 = fit(1); X = fit(2); y1 = fit(3); x0 = fit(4); yc = fit(5);

    E = sum((y - y0 .* tanh((x-x0)/X) - y1*((x-x0)/X) - yc).^2);
end

function [] = test_fit()
    x = [0.02:0.05:20];
    X = 2;
    y0 = 2*1e3;
    x0 = mean(x); %mean(x(:));
    y1 = 0.2*1e3;
    yc = 14e6;
    y = y0 * tanh((x-x0)/X) + y1 * ((x-x0)/X) + yc;% + rand(size(x)).*y0/10;

    [yy,xx,yy1,xx0,yc0] = runs.tanh_fit(x,y,1);

    disp(['y0 = ' num2str(yy) ' | Original = ' num2str(y0)]);
    disp(['yy1 = ' num2str(yy1) ' | Original = ' num2str(y1)]);
    disp(['X = ' num2str(xx) ' | Original = ' num2str(X)]);
    disp(['x0 = ' num2str(xx0) ' | Original = ' num2str(x0)]);
    disp(['yc = ' num2str(yc0) ' | Original = ' num2str(yc)]);
end
