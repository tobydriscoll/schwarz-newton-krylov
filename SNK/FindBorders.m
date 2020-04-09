function [border,boundary,interface] = FindBorders(degs,domain,bbox)
% Indices of the interface and boundary points.
%
% INPUT
%     degs: vector of grid sizes
%     domain: each row gives min,max values of domain in a dimension
%     bbox: each row gives min,max values of bounding box in a dimension
%
% OUTPUT
%     (Each is a structure with fields named all,north,south,east,west.
%      Field values are logical arrays.)
%     border: all points on the border of the domain
%     boundary: all points on the true boundary
%     interface: all points on the non-boundary border
%
%  Some care is given to not double-count points.
%

South = false(degs); South(:,1) = true;
North = false(degs); North(:,end) = true;
West = false(degs); West(1,:) = true;
East  = false(degs); East(end,:) = true;

border.all = South | North | East | West;
border.north = North;
border.south = South & ~North;
border.west = West & ~(North | South);
border.east = East & ~(North | South | West);

% Might have to make this more robust!
outer = false(degs);
if domain(1,1)==bbox(1,1)
	outer = outer | West;
end

if domain(1,2)==bbox(1,2)
	outer = outer | East;
end

if domain(2,1)==bbox(2,1)
	outer = outer | South;
end

if domain(2,2)==bbox(2,2)
	outer = outer | North;
end

boundary.all = outer(:);
boundary.west = outer(:) & West(:);
boundary.east = outer(:) & East(:) & ~boundary.west;
boundary.south = outer(:) & South(:) & ~(boundary.west | boundary.east);
boundary.north = outer(:) & North(:) & ~(boundary.west | boundary.east | boundary.south);

inner = false(degs);
inner(border.all) = ~outer(border.all);
interface.all = inner(:);

interface.west = inner(:) & West(:);
interface.east = (inner(:) & East(:)) & ~interface.west;
interface.south = inner(:) & South(:) & ~(interface.west | interface.east);
interface.north = inner(:) & North(:) & ~(interface.west | interface.east | interface.south);

end