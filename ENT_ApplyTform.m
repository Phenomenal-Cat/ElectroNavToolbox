function [X,Y,Z] = ENT_ApplyTform(Tform, x,y,z)

%=========================== ENT_ApplyTform.m =============================
% Apply the 4 x 4 transformation matrix 'Tform' to the input coordinates 
% [x,y,z].
%
% ELECTRONAV TOOLBOX
% Developed by Aidan Murphy, © Copyleft 2015, GNU General Public License
%==========================================================================

if all(size(Tform)==[3,3])
    Append = 0;
elseif all(size(Tform)==[4,4])
    Append = 1;
else
    error('Input ''Tform'' must either be a 3x3 or 4x4 matrix!');
end

if nargin == 2              %========== Matrix coordinate input
    InputSize = size(x);
    if find(size(x)==3)==1
        if Append == 1
            x = [x, ones(size(x,1),1)]';
        end
        verts = Tform*x;
    elseif find(size(x)==3)==2
        if Append==1
            x = [x, ones(size(x,1),1)];
        end
        verts = Tform*x';
    elseif size(x,1)==4
        verts = Tform*x;
    elseif size(x,2)==4
        verts = Tform*x';
    end
elseif nargin == 4          %========== Vector coordinate inputs
    if size(x,1)==1
        xyz = [x', y', z', ones(numel(z),1)];
    elseif size(x,2)==1
        xyz = [x, y, z, ones(numel(z),1)];
    end
    verts = Tform*xyz';
else
    error('Number of inputs to ENT_ApplyTform.m must be either 2 or 4!');
end

if nargout == 3
    X = verts(1,:);
    Y = verts(2,:);
    Z = verts(3,:);
elseif nargout == 1
    X = verts(1:3,:);
    if size(X)~= InputSize
        X = X';
    end
end

end