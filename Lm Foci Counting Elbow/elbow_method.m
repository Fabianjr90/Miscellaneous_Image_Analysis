function [IDX,C,SUMD,K]=elbow_method(X,numClusters,varExplained)

% This code computes the distorions of having a specific number of
% clusters in the data. It will then partition data into K clusters.
% IDX = cluster indices of each point.
% C = the K cluster centroid locations.
% SUMD = sums of point-to-centroid distances.
% K = the number of cluster centroids determined using ELBOW method.
% ELBOW method: computing the distortions under different cluster number
% counting from % 1 to numClusters. K is the cluster number corresponding
% to varExplained (e.g. 95%) of the % variance expained, which is the
% ratio of the between-group variance to the total variance.

figure
test_num=numClusters;
distortion=zeros(numClusters,1);
for k_temp=1:numClusters
    [~,~,sumd]=kmeans(X,k_temp,'emptyaction','drop');
    destortion_temp=sum(sumd);
    
    % try different number of clusters to find minimun disortion
    for test_count=2:test_num
        [~,~,sumd]=kmeans(X,k_temp,'emptyaction','drop');
        destortion_temp=min(destortion_temp,sum(sumd));
    end
    distortion(k_temp,1)=destortion_temp;
end

variance=distortion(1:end-1)-distortion(2:end);
distortion_percent=cumsum(variance)/(distortion(1)-distortion(end));
plot(distortion_percent,'b*--');
xlabel('Number of clusters','FontSize', 16)
ylabel('Percent variance explained by cluster number','FontSize', 16)

% We are looking to explain varExplained% of variance
[r,~]=find(distortion_percent>varExplained);
K=r(1,1)+1;

[IDX,C,SUMD]=kmeans(X,K);

end