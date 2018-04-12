function plotClass(X, label)

% This function will create scatter plots of data X using vector label

[~,n] = size(X);
if nargin == 1
    label = ones(n,1);
end
assert(n == length(label));

color = 'brgmcyk';
m = length(color);
c = max(label);

figure(gcf);
clf;
hold on;

for i = 1:c
    idc = label==i;
    scatter(X(1,idc),X(2,idc),36,color(mod(i-1,m)+1),'filled');
end

axis equal
grid on
hold off

end