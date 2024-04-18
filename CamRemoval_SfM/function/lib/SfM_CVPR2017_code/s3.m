function out = s3(varargin)
%function inds = s3(i)
%function out  = s3(R,i)
%function out  = s3(A,i,j)
%
% This function, when given an index i vector, converts it to the node
% indiecs (3*i-2,3*i-1,3*i) for each value of i. It is meant for accessing
% subblocks of block matrices.
%
% Uses:
% A(s3(1),s3(2)) retrieves the subblock A(1,2).
% A(s3([1 7]),s3(2:3)) retrieves the corresponding subblocks

if nargin == 1
    i = varargin{1};
    
    if size(i,1) == 1
        stampConst = [3 3 3];
        stampVar = [0 1 2];
    else
        stampConst = [3 3 3]';
        stampVar = [0 1 2]';
    end
    
    out = 1 + kron(ones(size(i)), stampVar) + kron((i-1), stampConst);
else
    %A = varargin{1};
    %i = varargin{2};    
    %if nargin >= 3
    %    j = varargin{3};
    %else
    %    j = ones(size(i));
    %end
    %out = A(s3(i),s3(j));
    
    if nargin >= 3
        out = varargin{1}(s3(varargin{2}),s3(varargin{3}));
    else
        out = varargin{1}(s3(varargin{2}(:)),s3(1));
    end
end

end
