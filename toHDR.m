% clear workspace
clear

% load sample image       
imgs = imread('data/office_1.jpg');
imgs(:,:,:,2) = imread('data/office_2.jpg');
imgs(:,:,:,3) = imread('data/office_3.jpg');
imgs(:,:,:,4) = imread('data/office_4.jpg');
imgs(:,:,:,5) = imread('data/office_5.jpg');
imgs(:,:,:,6) = imread('data/office_6.jpg');

n_images = 6;
sz = size(imgs);

% gamma decode images and erase the original ones to spare memory
g = 2.2;
linear_imgs = uint8( round(255*((single(imgs)/255).^g)) );

% load camera response curve (from HDRShop) and
% exposure times
run('data/curve.m');
run('data/exps.m');

% compute radiances for each image
rad_imgs = zeros(sz);
for id = 1:n_images
    % loop through channels
    for channel = 1:3
        exps = exp(C(imgs(:,:,channel,id)+1, channel));
        rad_imgs(:,:,channel,id) = reshape(exps/exp_time(id), sz(1:2));
    end
end
clearvars exps

% compute the average radiance map
% discard saturated or underexposed pixels
saturated = 254; underexposed = 1;
rad_img = zeros(sz(1:3));

% ------ discard the whole pixel if one channel is invalid
good = imgs(:,:,1,:) <= saturated & imgs(:,:,1,:) >= underexposed & ...
        imgs(:,:,2,:) <= saturated & imgs(:,:,2,:) >= underexposed & ...
        imgs(:,:,3,:) <= saturated & imgs(:,:,3,:) >= underexposed;

% Avoid NaNs: sometimes, all pixels will be invalid (pixel (579,194)
% in blue channel, for example, is zero in all images). In this case,
% I believe there's nothing to be done because we are uncertain of the
% actual radiance. In order to avoid NaNs, we'll just take the average
% of all pixel values: if all of them are black, in the end we'll have
% the dimmest possible value; if all of them are white, we'll have the
% brightest one, which is fair reasonable assuming we don't know
% exactly how bright they actually are. In the (unlikely) case where
% some pixels are underexposed and some are saturated, we'll end up
% computing a weighted average and the final value will be somewhere
% between the brightest and dimmest possible value.
weights = sum(good(:,:,1,:), 4);
must_change = (weights == 0);

for id = 1:n_images
    good(:,:,1,id) = good(:,:,1,id) + must_change;
end

% We need to recompute the weights, as some pixels might have changed
weights = sum(good(:,:,1,:), 4);
weighted = zeros(sz);

% when computing the average, weight the values using GOOD. This
% will set underexposed/saturated values to zero; then, just divide
% by the number of 1's in the GOOD matrix. This will compute the
% average of the good pixels only.
for channel = 1:3
    rad_img(:,:,channel) = sum(rad_imgs(:,:,channel,:) .* good, 4) ./ weights;
end

% ------ discard invalid pixels in a per-channel basis
%{
for channel = 1:3
    % this incredbly clever piece of code computes a BOOLEAN matrix
    % whose entries are 0 in case a value is underexposed/saturated
    % and 1 if they're fine.
    c = imgs(:,:,channel,:);
    good = reshape((c >= underexposed & c <= saturated), [sz(1) sz(2) 1 sz(4)]);
    
    % ---- avoid NaNs ----
    weights = sum(good(:,:,1,:), 4);
    must_change = (weights == 0);
    
    for id = 1:n_images
        good(:,:,1,id) = good(:,:,1,id) + must_change;
    end
    
    % ---- arithmetic average ----
    weights = sum(good(:,:,1,:), 4);
    weighted = sum(rad_imgs(:,:,channel,:) .* good, 4);
    rad_img(:,:,channel) = weighted ./ weights;
    
    % ---- median ----
    %rad_img(:,:,channel) = median(rad_imgs(:,:,channel,:), 4);
end
%}

% clean up to spare some memory
clearvars c
clearvars good
clearvars weighted
clearvars weights
clearvars rad_imgs
clearvars linear_imgs
clearvars must_change

% compare with MATLAB's HDR generation algorithm
files = {'office_1.jpg', 'office_2.jpg', 'office_3.jpg', ...
         'office_4.jpg', 'office_5.jpg', 'office_6.jpg'};
expTimes = [0.0333, 0.1000, 0.3333, 0.6250, 1.3000, 4.0000];
hdr = makehdr(files, 'RelativeExposure', expTimes./expTimes(1));

%-------------------------------------
%----------- Tone mapping ------------
%-------------------------------------
avg_I = zeros(3,1);
for c = 1:3
    avg_I(c) = exp(sum(sum(log(rad_img(:,:,c)+0.00001)))/numel(rad_img(:,:,c)));
end

L = 0.299*rad_img(:,:,1) + 0.586*rad_img(:,:,2) + 0.114*rad_img(:,:,3);
avg_L = exp(sum(sum(log(L+0.00001)))/numel(L));

alpha = 0.0; beta = 0.5;
gamma = 0.18;

ldr = zeros(size(rad_img));
for channel = 1:3
    I = rad_img(:,:,channel);
    
    %semi saturation constant
    Gr = alpha*I + (1-alpha)*L;
    Gr_ = alpha*avg_I(channel) + (1-alpha)*avg_L;
    ss = beta*Gr + (1-beta)*Gr_;
    
    scaled_lum = (gamma/avg_L) * I;
    
    %ldr(:,:,channel) = scaled_lum ./ (1 + scaled_lum);
    ldr(:,:,channel) = I ./ (I + ss);
end

%ldr = ldr.^(1/g);
%imtool( ldr );

%-------------
avg_I = zeros(3,1);
for c = 1:3
    avg_I(c) = exp(sum(sum(log(hdr(:,:,c)+0.00001)))/numel(hdr(:,:,c)));
end

L = 0.299*hdr(:,:,1) + 0.586*hdr(:,:,2) + 0.114*hdr(:,:,3);
avg_L = exp(sum(sum(log(L+0.00001)))/numel(L));

alpha = 0.0; beta = 0.5;
gamma = 0.18;

ldr = zeros(size(hdr));
for channel = 1:3
    I = hdr(:,:,channel);
    
    %semi saturation constant
    Gr = alpha*I + (1-alpha)*L;
    Gr_ = alpha*avg_I(channel) + (1-alpha)*avg_L;
    ss = beta*Gr + (1-beta)*Gr_;
    
    scaled_lum = (gamma/avg_I(channel)) * I;
    
    ldr(:,:,channel) = scaled_lum ./ (1 + scaled_lum);
    %ldr(:,:,channel) = I ./ (I + ss);
end

ldr = ldr.^(1/g);
imtool(ldr);