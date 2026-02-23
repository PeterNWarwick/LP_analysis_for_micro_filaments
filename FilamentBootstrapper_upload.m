%% Lp analysis from filament skeletons%%
clear all;close all;
%% PARAMETERS
inputFile   = 'examplefilename.tif'; % read tif stack with clean filament skeletons
numRuns     = 2;   % number of repeats                       !!!
N           =50;    % number of filaments per repeat          !!! 1000 or more
minPixels   = 2;    % remove very small filaments
minDistance = 2;    % min distance between skeletons
maxTries    = 1000;
rng('shuffle');
pixelSize = 0.248;      % / microns per pixel (pixel size)     !!! change
ds_interp = pixelSize*1; % / microns (interpolation step)      !!! 1 is a good start
ds_pix    = 1;        % step along filament in pixels
maxLag_um = 61;        % / microns, Lc,max
binWidth  = pixelSize;      % / microns (bin size)
%% READ inputFile
info      = imfinfo(inputFile);
numSlices = numel(info);
imgSize   = [info(1).Height, info(1).Width];
filaments = {};
fid = 1;
for s = 1:numSlices
    img = imread(inputFile, s) > 0;
    CC = bwconncomp(img, 8);
    for i = 1:CC.NumObjects
        pix = CC.PixelIdxList{i};
        if numel(pix) >= minPixels
            [y,x] = ind2sub(imgSize, pix);
            filaments{fid}.coords = [x(:), y(:)];
            filaments{fid}.length = numel(pix);
            fid = fid + 1;
        end
    end
end
fprintf('Total number of filaments found in stack: %d\n', numel(filaments));
%% Bootstrap Loop
fitresult1={};
gof1={};
fitresult2={};
gof2={};
fitresult3={};
gof3={};
for run = 1:numRuns
    fprintf('Bootstrip run %d / %d\n', run, numRuns);
    % select randomly
    nSel = min(N, numel(filaments));
    idx  = randperm(numel(filaments), nSel);
    selected = filaments(idx);
    % empty target
    targetSkel = false(imgSize);
    % cleaning
    all_s   = [];
    all_cos = [];
    maxLag_pix = round(maxLag_um / pixelSize);
    %% loop over filaments
    for i = 1:numel(selected)
        coords_raw = selected{i}.coords;
        coords=orderCoordsAlongFilament(coords_raw);
        %% interpolate pixelated filament contour
        % pixel coordinates
        dcoords = diff(coords,1,1);
        s_raw = [0; cumsum(vecnorm(dcoords,2,2))];
        % clean (branches, steplength)
        dcoords = diff(coords,1,1);
        stepLen = vecnorm(dcoords,2,2);
        % interpolate
        s_fine = 0:ds_interp:s_raw(end);
        coords_interp = interp1(s_raw, coords, s_fine, 'spline');
        % new tangents
        dcoords_interp = diff(coords_interp,1,1);
        t = dcoords_interp ./ vecnorm(dcoords_interp,2,2);
        % smoothing with moving average (3–5 reasonable)
        t = smoothdata(t,1,'movmean',5);
        t = t ./ vecnorm(t,2,2);  % normalize
        % tangent arc length
        s_tan_um = (s_fine(1:end-1) + ds_interp/2) * pixelSize;
        % cosine correlation
        nT = numel(s_tan_um);  % number of tangents
        for lag = 1:nT-1
            dots = sum(t(1:end-lag,:) .* t(1+lag:end,:), 2);
            s_sep = s_tan_um(1+lag:nT) - s_tan_um(1:nT-lag);
            % safety check
            if numel(dots) ~= numel(s_sep)
                error('mismatch of s_sep and dots lengths detected')
            end
            all_s   = [all_s;   s_sep(:)];
            all_cos = [all_cos; dots(:)];
        end
    end
    %% Binning
    edges = 0:binWidth:maxLag_um;
    centers = edges(1:end-1) + binWidth/2;
    C_mean = nan(size(centers)); % mean cosine correlation
    C_std  = nan(size(centers)); % standard deviation of cosine correlation
    C_stdm  = nan(size(centers)); % standard deviation of the mean of cosine correlation
    C_n    = zeros(size(centers));
    for b = 1:numel(centers)
        inBin = all_s >= edges(b) & all_s < edges(b+1);
        if any(inBin)
            C_mean(b) = mean(all_cos(inBin));
            C_std(b)  = std(all_cos(inBin));
            C_stdm(b)  =C_std(b)./(length(all_cos(inBin))).^0.5;
            C_n(b)    = sum(inBin);
        end
    end
    %% Plotting
    figure; hold on
    errorbar(centers, C_mean, C_std, 'o-', 'LineWidth', 1.5)
    xlabel(['\it{s}\rm{} [',char(181),'m]'],'Interpreter','tex','fontsize',18,'Fontname','Helvetica');
    ylabel('\langle cos θ(s) \rangle','fontsize',18,'Fontname','Helvetica');
    title('Mean cosine correlation')
    grid on
    box on ; hold off
    %% Fitting with three different methods
    ft1 = fittype('exp(-x/(2*Lp))', ...
        'independent','x','coefficients',{'Lp'});
    start1 = [10];
    ft2 = fittype('A*exp(-x/(2*Lp))', ...
        'independent','x','coefficients',{'A','Lp'});
    start2 = [1 10];
    ft3 = fittype('exp(-x/(2*Lp))', ...
        'independent','x','coefficients',{'Lp'});
    start3 = [10];
    x1=centers;
    y1=C_mean;
    x2=centers;
    y2=C_mean;
    x3=centers;
    y3=C_mean;
    sigma1=C_std;
    sigma2=C_std;
    sigma3=C_std;
    sigma1m=C_stdm;
    sigma2m=C_stdm;
    sigma3m=C_stdm;
    xcut=pixelSize*35;
    k1=1;
    w1 = 1 ./ (sigma1m.^2);
    w2 = 1 ./ (sigma2m.^2);
    for i=1:length(x3)
        w3(i) = 1 ./ (sigma3m(i).^2).*xcut.*(1/(1+exp(-k1*(x3(i)-xcut))));
    end
    % remove infinite or NaN weights
    invalid = isnan(w1) | isinf(w1);
    w1(invalid) = [];
    x1(invalid) = [];
    y1(invalid) = [];
    sigma1(invalid) = [];
    sigma1m(invalid) = [];
    invalid = isnan(w2) | isinf(w2);
    w2(invalid) = [];
    x2(invalid) = [];
    y2(invalid) = [];
    sigma2(invalid) = [];
    sigma2m(invalid) = [];
    invalid = isnan(w3) | isinf(w3);
    w3(invalid) = [];
    x3(invalid) = [];
    y3(invalid) = [];
    sigma3(invalid) = [];
    sigma3m(invalid) = [];
    opts1 = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', start1);
    if ~isempty(w1)
        opts1.Weights = w1;
    end
    opts2 = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', start2);
    if ~isempty(w2)
        opts2.Weights = w2;
    end
    opts3 = fitoptions('Method','NonlinearLeastSquares', 'StartPoint', start3);
    if ~isempty(w3)
        opts3.Weights = w3;
    end
    [fitresult, gof] = fit(x1', y1', ft1, opts1);
    fitresult1{run}=fitresult;
    gof1{run}=gof;
    [fitresult, gof] = fit(x2', y2', ft2, opts2);
    fitresult2{run}=fitresult;
    gof2{run}=gof;
    [fitresult, gof] = fit(x3', y3', ft3, opts3);
    fitresult3{run}=fitresult;
    gof3{run}=gof;
    %% Plot bare decay
    x_fit = linspace(0, maxLag_um, 500)';
    y_fit = feval(fitresult1{run}, x_fit);
    ci = confint(fitresult1{run},0.95);
    paramNames = coeffnames(fitresult1{run});
    paramVals  = coeffvalues(fitresult1{run});
    y_lower = y_fit;
    y_upper = y_fit;
    Lp_lower = ci(1,strcmp(paramNames,'Lp'));
    Lp_upper = ci(2,strcmp(paramNames,'Lp'));
    y_lower = exp(-x_fit/(2*Lp_upper));
    y_upper = exp(-x_fit/(2*Lp_lower));
    figure; hold on; box on;
    title('Fitting')
    O1=plot(x1, y1, 'ko', 'MarkerFaceColor','k', 'DisplayName','\langle cos θ(s) \rangle');
    fill([x_fit; flipud(x_fit)], [y_upper; flipud(y_lower)], ...
        [0.8 0.8 1], 'EdgeColor','none','FaceColor','blue', 'FaceAlpha',0.5, 'DisplayName','95%-CI');
    F1=plot(x_fit, y_fit, 'b-', 'LineWidth',2.5, 'DisplayName','bare fit');
    xlabel(['\it{s}\rm{} [',char(181),'m]'],'Interpreter','tex','fontsize',18,'Fontname','Helvetica');
    ylabel('\langle cos θ(s) \rangle','fontsize',18,'Fontname','Helvetica');
    legend('show','Location','northeast');
    legend('\langle cos θ(s) \rangle','95% confidence interval','bare fit','Fontsize',12);
    grid on;
    %% Plot amplitude decay
    x_fit = linspace(0, maxLag_um, 500)';
    y_fit = feval(fitresult2{run}, x_fit);
    ci = confint(fitresult2{run},0.95);
    paramNames = coeffnames(fitresult2{run});
    paramVals  = coeffvalues(fitresult2{run});
    y_lower = y_fit;
    y_upper = y_fit;
    A_lower  = ci(1,strcmp(paramNames,'A'));
    A_upper  = ci(2,strcmp(paramNames,'A'));
    Lp_lower = ci(1,strcmp(paramNames,'Lp'));
    Lp_upper = ci(2,strcmp(paramNames,'Lp'));
    y_lower = A_lower * exp(-x_fit/(2*Lp_upper));
    y_upper = A_upper * exp(-x_fit/(2*Lp_lower));
    fill([x_fit; flipud(x_fit)], [y_upper; flipud(y_lower)], ...
        [0.8 0.8 1], 'EdgeColor','none','FaceColor','yellow', 'FaceAlpha',0.5, 'DisplayName','95%-CI');
    F2=plot(x_fit, y_fit, 'y-', 'LineWidth',2.5, 'DisplayName','Amplitude fit');
    %% Plot midrange-weighted decay
    x_fit = linspace(0, maxLag_um, 500)';
    y_fit = feval(fitresult3{run}, x_fit);
    ci = confint(fitresult3{run},0.95);
    paramNames = coeffnames(fitresult3{run});
    paramVals  = coeffvalues(fitresult3{run});
    y_lower = y_fit;
    y_upper = y_fit;
    Lp_lower = ci(1,strcmp(paramNames,'Lp'));
    Lp_upper = ci(2,strcmp(paramNames,'Lp'));
    y_lower = exp(-x_fit/(2*Lp_upper));
    y_upper = exp(-x_fit/(2*Lp_lower));
    fill([x_fit; flipud(x_fit)], [y_upper; flipud(y_lower)], ...
        [0.8 0.8 1], 'EdgeColor','none','FaceColor','red', 'FaceAlpha',0.5, 'DisplayName','95%-CI');
    F3=plot(x_fit, y_fit, 'r-', 'LineWidth',2.5, 'DisplayName','Mid-range weighted fit');
    hold off
end
%% Output
disp('Fit results');
Lp1=mean(cellfun(@(f) f.Lp, fitresult1));
Lp1_std=std(cellfun(@(f) f.Lp, fitresult1));
fprintf('e^(s/(2Lp)), w=1/Sm %2f ± %.2f\n', Lp1, Lp1_std);
Lp2=mean(cellfun(@(f) f.Lp, fitresult2));
Lp2_std=std(cellfun(@(f) f.Lp, fitresult2));
fprintf('Ae^(s/(2Lp)), w=1/Sm %2f ± %.2f\n', Lp2, Lp2_std);
Lp3=mean(cellfun(@(f) f.Lp, fitresult3));
Lp3_std=std(cellfun(@(f) f.Lp, fitresult3));
fprintf('e^(s/(2Lp)), w=(1/Sm)·B(1/(1+e^(−1·(si−B)))) %2f ± %.2f\n', Lp3, Lp3_std);







%% Functions
function coords_sorted = orderCoordsAlongFilament(coords)
N = size(coords,1);
used = false(N,1);
coords_sorted = zeros(size(coords));
% starting point (highest distance)
center = mean(coords,1);
d2 = sum((coords - center).^2,2);
[~, idx] = max(d2);
coords_sorted(1,:) = coords(idx,:);
used(idx) = true;
for k = 2:N
    last = coords_sorted(k-1,:);
    d = vecnorm(coords(~used,:) - last, 2, 2);
    idx_free = find(~used);
    [~, minIdx] = min(d);
    nextIdx = idx_free(minIdx);
    coords_sorted(k,:) = coords(nextIdx,:);
    used(nextIdx) = true;
end
end
