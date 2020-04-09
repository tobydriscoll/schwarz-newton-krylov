function [boundary_s,boundary,interface,interface_s,border,border_s] = FindBoundaryIndex2DSides(degs,domain,bbox)

% Indices of the interface and boundary points. 

%  Some care is given to not double count indicies.
% 
%  out_border, in_border: logical array for points in the outer boundary
%  and inner boundary.
%
%  out_border_s, in_border_s: cell array for west,east,south,north sides (in
%  that order) for the outer boundary and interface.

boundary = false(degs);
interface = false(degs);

South = false(degs); South(:,1) = true;

North = false(degs); North(:,end) = true;

West = false(degs); West(1,:) = true;

East  = false(degs); East(end,:) = true;

border = South | North | East | West;

border_s{4} = North;
border_s{3} = South & ~ North;
border_s{1} = West & ~(North | South);
border_s{2} = East & ~(North | South + West);

%Might have to make this more robust!
if domain(1,1)==bbox(1,1)
    boundary = boundary | West;
end

if domain(1,2)==bbox(1,2)
    boundary = boundary | East;
end

if domain(2,1)==bbox(2,1)
    boundary = boundary | South;
end

if domain(2,2)==bbox(2,2)
    boundary = boundary | North;
end

interface(border) = ~boundary(border);

boundary = boundary(:);
interface = interface(:);
boundary_s = cell(4,1);

boundary_s{1} = boundary & West(:);
boundary_s{2} = boundary & East(:) & ~ boundary_s{1};

boundary_s{3} = boundary & South(:) & ~(boundary_s{1} | boundary_s{2});
boundary_s{4} = boundary & North(:) & ~(boundary_s{1} | boundary_s{2} | boundary_s{3});

interface_s{1} = interface & West(:);
interface_s{2} = (interface & East(:)) & ~ interface_s{1};

interface_s{3} = interface & South(:) & ~(interface_s{1} | interface_s{2});
interface_s{4} = interface & North(:) & ~(interface_s{1} | interface_s{2} | interface_s{3});

end