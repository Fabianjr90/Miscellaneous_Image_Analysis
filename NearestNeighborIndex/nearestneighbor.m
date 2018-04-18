function idx = nearestneighbor(X)
%NEARESTNEIGHBOUR    find nearest neighbours
%   IDX = NEARESTNEIGHBOUR(X) finds the nearest neighbour by Euclidean
%   distance to each point (column) in X from X. X is a matrix with points
%   as columns. IDX is a vector of indices into X, such that X(:, IDX) are
%   the nearest neighbours to X. e.g. the nearest neighbour to X(:, 2) is
%   X(:, IDX(2))


% pre-define for speed
idx = zeros(1, size(X, 2));

% Loop through the set of points X, finding the neighbors
Y = zeros(size(X));
for ii = 1:size(X, 2)
    x = X(:, ii);
    
    % This is the faster than using repmat based techniques such as
    % Y = X - repmat(x, 1, size(X, 2))
    for i = 1:size(Y, 1)
        Y(i, :) = X(i, :) - x(i);
    end
    
    % Find the closest points, and remove matches beneath a radius
    dSq = sum(abs(Y).^2, 1);
    iRad = find(dSq < inf);
    iSorted = iRad(minn(dSq(iRad), 2));
    iSorted = iSorted(2:end);
    
    % Remove any bad ones
    idx(1:length(iSorted), ii) = iSorted';
end

idx(all(idx==0,2),:)=[];
if isvector(idx)
    idx = idx(:)';
end

end


%MINN   find the n most negative elements in x, and return their indices
%  in ascending order
function I = minn(x, n)

% Make sure n is no larger than length(x)
n = min(n,length(x));

% Sort the first n
[xsn, I] = sort(x(1:n));

% Go through the rest of the entries, and insert them into the sorted block
for jj = (n+1):length(x)
    j = n;
    while j > 0 && x(jj) < xsn(j)
        j = j - 1;
    end

    if j < n
        % x(i) should go into the (j+1) position
        xsn = [xsn(1:j), x(jj), xsn((j+1):(n-1))];
        I   = [I(1:j), jj, I((j+1):(n-1))];
    end
end

end
