clear all;
close all;
%% Parameters
imgsize    = 512;      % image size in pixel
pixelSize  = 0.248;    % pixel size in um/pixel
Lp         = 10;       % filament persistence length in um
nFil       = 35;       % number of filaments in per image
Lmean      = 7.5;      % mean filament length
Lstd       = 0.9;      % lognormal distribution of filamnents
noiseLevel = 0.05;     % intensity noise between 0 and 1
filWidth   = 3;        % filament width in pixel
nImages    = 50;       % number of images
saveImages = true;     % saving images (true or false)
outDir = fullfile(pwd, 'simulated_filaments'); % path of new folder for saved images
if saveImages && ~exist(outDir, 'dir')
    mkdir(outDir);
end
%% Simulate filaments
for i = 1:nImages
    fprintf('Simulated_image_%d of %d...\n', i, nImages);
    img = simulate_filaments_image(imgsize, pixelSize, Lp, nFil, ...
                                   Lmean, Lstd, noiseLevel, filWidth);
    figure; imshow(img, []);
    title(sprintf('simulated_filament_image_%d (L_p = %.1f µm)', i, Lp));
    if saveImages
        filename = sprintf('filaments_Lp%.1f_%03d.png', Lp, i);
        imwrite(img, fullfile(outDir, filename));
    end
end
%% Functions
%% create filaments
function img = simulate_filaments_image(imgsize, pixelSize, Lp, nFil, Lmean, Lstd, noiseLevel, filWidth)
img = zeros(imgsize);
Lc_vec = lognrnd(log(Lmean), Lstd, [nFil, 1]);  % lognormal size distribution
for i = 1:nFil
    Lc = Lc_vec(i);
    N = round(Lc / (pixelSize * 0.5)); 
    if N < 5, continue; end
    ds = Lc / N;
    x0 = rand() * imgsize; % random start
    y0 = rand() * imgsize;
    dtheta = randn(N,1) * sqrt(ds / Lp); % curvature governed by Lp
    theta = cumsum(dtheta) + 2*pi*rand();
    x = x0 + cumsum(cos(theta) * ds / pixelSize);
    y = y0 + cumsum(sin(theta) * ds / pixelSize);
    valid = (x > 1 & x < imgsize & y > 1 & y < imgsize); % limit to field of view
    x = x(valid); y = y(valid);
    for k = 2:length(x)  % draw lines
        img = drawThickLine(img, [x(k-1), y(k-1)], [x(k), y(k)], filWidth, 0.9);
    end
end
img = imgaussfilt(img, 0.5); % noise and smoothing
img = img + noiseLevel * randn(size(img));
img = mat2gray(img);
end

%% draw thick lines
function img = drawThickLine(img, p1, p2, width, intensity)
[xv, yv] = bresenham(round(p1(1)), round(p1(2)), round(p2(1)), round(p2(2))); % draw line with variable size
for k = 1:length(xv)
    if xv(k) > 0 && xv(k) <= size(img,2) && yv(k) > 0 && yv(k) <= size(img,1)
        [xx, yy] = meshgrid(1:size(img,2), 1:size(img,1));
        mask = ((xx - xv(k)).^2 + (yy - yv(k)).^2) <= (width/2)^2;
        img(mask) = min(1, img(mask) + intensity);
    end
end
end
%% Bresenham line
function [x, y] = bresenham(x1, y1, x2, y2)
x1 = round(x1); y1 = round(y1); x2 = round(x2); y2 = round(y2);
dx = abs(x2 - x1); dy = abs(y2 - y1);
sx = sign(x2 - x1); sy = sign(y2 - y1);
err = dx - dy;
x = []; y = [];
while true
    x = [x x1]; y = [y y1];
    if x1 == x2 && y1 == y2, break; end
    e2 = 2*err;
    if e2 > -dy, err = err - dy; x1 = x1 + sx; end
    if e2 < dx, err = err + dx; y1 = y1 + sy; end
end
end