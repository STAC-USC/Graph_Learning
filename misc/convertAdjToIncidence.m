
% Input: 
%       Adj: 0/1 adjacency matrix (no self loops) given the *undirected* graph connection pattern
% Output: 
%       Incidence : Incidence matrix
%       fromNodes : list of *from* nodes
%       toNodes   : list of *to* nodes
function [Incidence,fromNodes,toNodes] = convertAdjToIncidence(Adj)

% take lower triangular part
upper = ~logical(tril(Adj));
Adj(upper) = 0;

% Find edges to design
[fromNodes,toNodes] = find(Adj==1);
nNodes = size(Adj,1);
nLinks = length(fromNodes);

% create Incidence Matrix
Incidence = zeros(nNodes,nLinks);
for l=1:nLinks
    i=fromNodes(l); j=toNodes(l);
    Incidence(i,l) = 1; Incidence(j,l) = -1;
end

end