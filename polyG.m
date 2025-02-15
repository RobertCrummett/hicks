function [G] = polyG(x,y,order)
% Construct the least squares polynomial operator of a certain order
G = ones(size(x));
for i=1:order
    for j=0:i
        A = x.^(j).*y.^(i-j);
        G = [G, A];
    end
end
end