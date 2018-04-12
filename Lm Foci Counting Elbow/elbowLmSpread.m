% Goal is to measure the area of individual infection foci in a series of
% microscopy images. To do this, we first read in an image showing 
% L. monocytogenes infection foci. The user will also provide the maximum
% number of of clusters to be explored by the algorithm. The code will find
% the optimal number of clusters (minimizes the intercluster distances).
% Then, we measure the area of each cluster, which is proportional
% to the efficiency of Listeria spread.

% The vector convHullAreas contains the final data.

clc;clear;close all;
rng(2018)

%% Code starts here

% Step 1: read in raw image and display it
I = imread('Lm_RFP.tif');
I2 = double(I);
I3 = mat2gray(I2);
imshow(I3,[])
title('\fontsize{16}Raw Data')


% Step 2: binarize image, dilate it and clean it (these numbers have to be
% tweaked depending on the nature of the image). Also, elements at the
% edges have been removed bc they are not good for the analysis (to find
% their area)
bw = imbinarize(I3);
se = strel('disk',20);
bw2 = imdilate(bw,se);
bw3 = bwareaopen(bw2,1e4);
bw4 = imclearborder(bw3);
bw5 = logical(bw.*bw4);
figure,imshow(bw5)
title('\fontsize{16}Cleaned Up Data')


% Step 3: identify centroids of each point (individual bacteria)
s = regionprops(bw5,'centroid');
X = zeros(2,size(s,1));
for i=1:size(s)
    X(1,i) = s(i).Centroid(1);
    X(2,i) = s(i).Centroid(2);
end


% Step 4: provide number of clusters and varExplained. The algorithm will
% test clusters from 1 to numClusters, so think of this as an upper limit.
% Then, set up matrix of centroids and perform k-means clustering
% best_kmeans takes in matrix X' and number of k clusters to explore.
% Each cluster will be displayed with a different color
numClusters = 10;
varExplained = 0.95;
X = [X(1,:);-X(2,:)];
[ClusterIndex,centroids] = elbow_method(X',numClusters,varExplained);
figure,plotClass(X,ClusterIndex)


% Step 5: calculate area of convex hull of each of the clusters identified
% in the previous step
convHullAreas = zeros(1,max(ClusterIndex));
figure,imshow(imdilate(bw5,strel('disk',2)))
hold on
for ii=1:max(ClusterIndex)
    newClusterInd = find(ii==ClusterIndex);
    clusterX = X(1,newClusterInd)';
    clusterY = X(2,newClusterInd)';
    [k,convHullAreas(ii)] = convhull(clusterX,clusterY);
    plot(clusterX(k),-clusterY(k),'r-', 'LineWidth', 4)
end


% Step 6: display data
T = array2table(convHullAreas','VariableNames',{'AreaConvexHull'});
fprintf('\n')
disp(T)