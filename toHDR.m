% clear workspace
clear

% load sample image
imgs = imread('data/office_1.jpg');
imgs(:,:,:,2) = imread('data/office_2.jpg');
imgs(:,:,:,3) = imread('data/office_3.jpg');
sz = size(imgs);

% gamma decode images
g = 2.2;
linear_img = uint8( 255*((single(imgs)/255).^g) );

% load camera response curve (from HDRShop) and
% exposure times
run('data/curve.m');
run('data/exps.m');

% compute radiances for each image
rad_imgs = zeros(sz);

% loop through images
for id = 1:3
    % loop through channels
    for channel = 1:3
        exps = exp(C(linear_img(:,:,channel,id)+1, channel));
        rad_imgs(:,:,channel,id) = reshape(exps/exp_time(id), sz(1:2));
    end
end
clearvars exps

% compute the average radiance map
% discard saturated or underexposed pixels
saturated = 250; underexposed = 5;

rad_img = zeros(sz(1:3));

% loop through channels
for channel = 1:3
    % this incredbly clever piece of code computes a BOOLEAN matrix
    % whose entries are 0 in case a value is underexposed/saturated
    % and 1 if they're fine.
    c = imgs(:,:,channel,:);
    good = reshape((c > underexposed & c < saturated), [sz(1) sz(2) 1 sz(4)]);
    
    % when computing the average, weight the values using GOOD. This
    % will set underexposed/saturated values to zero; then, just divide
    % by the number of 1's in the GOOD matrix. This will compute the
    % average of the good pixels only.
    weighted = (rad_imgs(:,:,channel,:) .* good);
    weights = sum(good(:,:,channel,:), 4);
    
    rad_img(:,:,channel) = weighted ./ weights;
end

clearvars c