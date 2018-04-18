clear;clc;close all;

%%% The nearest neighbor index (NNI) measures the degree of clustering
%%% of the data. If 0<NNI<1, then data are clustered. If NNI~1, then data
%%% are random. If NNI~2, then data are organized in an array. The
%%% formula for this index is NNI = distObs/(0.5*sqrt(areaObs/n)) where
%%% distObs = the shortest distance of each element (bacterium) to its
%%% nearest neighbor, areaObs = the total area occupied by the data (in our
%%% case it corresponds to the area occupied by the host nuclei), and n =
%%% the number of data points, i.e. bacteria

imageName = 'Test.tif';
bacteriaChannel = 1;
nucleiChannel = 2;

%% Step 1: calculate number of bacteria

% 1a: binarize image using Otsu's method
I = imread(imageName,bacteriaChannel);
I2 = double(I);
I3 = mat2gray(I2);
bw = imbinarize(I3);

% 1b: segment image using the watershed transform
D = -bwdist(~bw);
include = D > -50;
D_min = imextendedmin(D,1).*include;
D_imp = imimposemin(D,D_min);
D_imp(D==0)=-Inf;
L = watershed(D_imp);

% 1c: cut out the background segments and clean tiny segments
L_fore = immultiply(L,bw);
L_fore_bw_clean = bwareaopen(L_fore,5);
L_fore = immultiply(L_fore,L_fore_bw_clean);

% 1d: calculate number of bacteria and show segmented image
[labeled_image,number_bacteria] = bwlabel(L_fore);
se = strel('disk',2);
final = imdilate(labeled_image,se);
RGB = label2rgb(final,'parula','w','shuffle');
figure,imshow(RGB)


%% Step 2: get distance of each bacterium to its nearest neighbor (distObs)

% 2a: find centroids of each bacterium
s = regionprops(final,'centroid');
centroids = cat(1,s.Centroid)';
idx = nearestneighbor(centroids);
distances = zeros(1,length(idx));

% 2b: calculate each bacterium's distance to nearest neighbor
for ii = 1:length(idx)
    distances_sq = (centroids(1,ii)-centroids(1,idx(ii)))^2 + ...
        (centroids(2,ii)-centroids(2,idx(ii)))^2;
    distances(1,ii) = distances_sq^0.5;
end

% 2c: calculate average distance observed (to nearest neighbor)
distObserved = mean(distances);


%% Step 3: calculate areaObserved (number of pixels of host nuclei)

% read in image, binarize it, and count the pixels
input_nuclei = imread(imageName,nucleiChannel);
I_nuc = double(input_nuclei);
I2_nuc = mat2gray(I_nuc);
bw_nuc = imbinarize(I2_nuc);
se = strel('disk',10);
bw_nuc2 = imdilate(bw_nuc,se);
pixel_nuc = double(bw_nuc2);
areaObserved = sum(pixel_nuc(:));


%% Step 4: calculate Nearest Neighbor Index (NNI)

myNNI = distObserved/(0.5*sqrt(areaObserved/number_bacteria));
disp(myNNI)
