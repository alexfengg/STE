function snorm = diffsnorm(varargin)
%DIFFSNORM  2-norm accuracy of an approx. to a matrix.
%
%
%   snorm = diffsnorm(A,U,S,V)  computes an estimate snorm of the spectral norm
%           (the operator norm induced by the Euclidean vector norm) of A-USV',
%           using 20 iterations of a power method started with a random vector.
%
%   snorm = diffsnorm(A,U,S,V,its)  computes an estimate snorm of the spectral
%           norm (the operator norm induced by the Euclidean vector norm)
%           of A-USV', using its iterations of a power method with a random
%           starting vector; its must be a positive integer.
%
%   snorm = diffsnorm(funcT,funcTt,params,U,S,V)
%   snorm = diffsnorm(funcT,funcTt,params,U,S,V,its)
%   These forms take the names of two functions that compute the
%   application of the operators T and T-transpose on a set of 
%   column vectors. "params" are any additional parameters required by T
%   and Tt. Use [] if no such paramters are required. The signature of T
%   and Tt must be either T(v) and Tt(v), or, T(params,v) and Tt(params,v).
%   When called without a column vector as its last parameter, T must
%   return m,n,cflag, where m is the number of rows in the matrix of the
%   operator, n is the number of columns, and cflag is 1 if the operator is
%   complex and 0 otherwise.
%
%   Increasing its improves the accuracy of the estimate of the spectral norm
%   of A-USV'.
%
%   
%
%   Note: DIFFSNORM invokes RAND. To obtain repeatable results, invoke
%         RANDN('state',j) with a fixed integer j before invoking DIFFSNORM.
%
%
%   inputs (the first four are required):
%   A -- first matrix in A-USV' whose spectral norm is being estimated
%   U -- second matrix in A-USV' whose spectral norm is being estimated
%   S -- third matrix in A-USV' whose spectral norm is being estimated
%   V -- fourth matrix in A-USV' whose spectral norm is being estimated
%   its -- number of iterations of the power method to conduct;
%          its must be a positive integer, and defaults to 20
%
%   output:
%   snorm -- an estimate of the spectral norm of A-USV' (the estimate fails to
%            be accurate with exponentially low probability as its increases;
%            see references 1 and 2 below.)
%
%
%   References:
%   [1] Jacek Kuczynski and Henryk Wozniakowski,
%       Estimating the largest eigenvalues by the power and Lanczos methods
%       with a random start, SIAM Journal on Matrix Analysis and Applications,
%       13 (4): 1094-1122, 1992.
%   [2] Edo Liberty, Franco Woolfe, Per-Gunnar Martinsson, Vladimir Rokhlin,
%       and Mark Tygert, Randomized algorithms for the low-rank approximation
%       of matrices, Proceedings of the National Academy of Sciences (USA),
%       104 (51): 20167-20172, 2007. (See the appendix.)
%   [3] Franco Woolfe, Edo Liberty, Vladimir Rokhlin, and Mark Tygert,
%       A fast randomized algorithm for the approximation of matrices,
%       Applied and Computational Harmonic Analysis, 25 (3): ????-????, 2008.
%       (See Section 3.4.)
%
%
%   See also NORMEST, NORM.
%
%   Yoel Shkolnisky, October 2008
%   Based on implementation by Mark Tygert.


global T_func  Tt_func

if(nargin < 4)
  error('MATLAB:pca:malformedInput',...
        'There must be at least 4 inputs.')
end

Op=varargin{1};

if ischar(Op)  % Decide if a matrix or function handle is given,
              % and parse the input accordingly.
    [T_func,Tt_func,params,U,S,V,its]=checkInputsT(varargin{:});
    if isempty(params) % Retrieve the dimensions of T.
        [m n,cflag] = feval(T_func); 
    else
        [m n,cflag] = feval(T_func,params); 
    end
    Aflag=0; % We are using function handles
elseif isnumeric(Op) % Matrix
    [A,U,S,V,its]=checkInputsA(varargin{:});    
    [m n] = size(A); % Retrieve the dimensions of A.
    Aflag=1; % We are using an explicit matrix
    
    cflag=0;
    if ~isreal(A)
        cflag=1;
    end
    
else
    error('MATLAB:pca:malformedInput',...
        'Input 1 must be a matrix or a name of a Matlab function.')
end

checkUSV(m,n,U,S,V); % Check that U, S, and V have the correct dimensions.

if(numel(its) ~= 1)
  error('MATLAB:diffsnorm:malformedInput',...
        'its must be a scalar.')
end

if(its <= 0)
  error('MATLAB:diffsnorm:malformedInput',...
        'its must be > 0.')
end


%
% Check the number of outputs.
%
if(nargout > 1)
  error('MATLAB:diffsnorm:malformedOutput',...
        'There must be at most 1 output.')
end

%
% Define functions for uniform syntax
%
if Aflag
    T=@(x,params) A*x;
    Tt=@(x,params) (A')*x;
    params=[];
else
    T=@T_tmp;
    Tt=@Tt_tmp;
end

%
% Construct the adjoints Ua of U, Sa of S, and Va of V.
%
Ua = U';
Sa = S';
Va = V';

%
% Generate a random vector x.
%
if(~cflag && isreal(U) && isreal(S) && isreal(V))
  x = randn(n,1);
end

if(cflag || ~isreal(U) || ~isreal(S) || ~isreal(V))
  x = randn(n,1) + i*randn(n,1);
end

x = x/norm(x);

%
% Run its iterations of the power method, starting with the random x.
%
for it = 1:its
%
% Set y = (A-USV')x.
%
  y = Va*x;
  y = S*y;
  y = U*y;
  y = T(x,params)-y;
%
% Set x = (A'-VS'U')y.
%
  x = Ua*y;
  x = Sa*x;
  x = V*x;
  x = Tt(y,params)-x;
%
% Normalize x, memorizing its Euclidean norm.
%
  snorm = norm(x);
  if(snorm == 0)
    break;
  end
  x = x/snorm;
end

snorm = sqrt(snorm);

clear x y Ua Sa Va;

end

%% checkInputsT
% First argument is a function handle. Functions for computing the operator
% and its transpose are given. Parse the input arguments accordingly.
function [T_func,Tt_func,params,U,S,V,its]=checkInputsT(varargin)

narg=numel(varargin);

% Check the number of inputs.

if(narg < 6)
    error('MATLAB:pca:malformedInput',...
        'When providing functions, there must be at least 6 inputs. Use [] if needed.')
end

if(narg > 7)
    error('MATLAB:pca:malformedInput',...
        'When providing functions, there must be at most 7 inputs.')
end

T_func=varargin{1};
Tt_func=varargin{2};
params=varargin{3};
U=varargin{4};
S=varargin{5};
V=varargin{6};

% Check arguments type
if ~ischar(T_func) || ~ischar(Tt_func)
    error('MATLAB:pca:malformedInput',...
        'First two input must be the names of the operator and its adjoint.')
end
    
% Set the inputs k, its, and l to default values, if necessary.
if(narg == 6)
    its=20;
else
    its=varargin{7};
    assertInt(its,7);
end

end


%% checkInputsA
% The matrix A is given explicitly. Parse the input arguments accordingly.
function [A,U,S,V,its]=checkInputsA(varargin)

narg=numel(varargin);

% Check the number of inputs.

if (narg<4)
        error('MATLAB:pca:malformedInput',...
        'When providing a matrix, there must be at least 4 inputs.')
end

if(narg > 5)
    error('MATLAB:pca:malformedInput',...
        'When providing a matrix, there must be at most 5 inputs.')
end

A=varargin{1};
U=varargin{2};
S=varargin{3};
V=varargin{4};

% Check argument type
if(~isnumeric(A))
    error('MATLAB:pca:malformedInput',...
        'Input 1 must be a floating-point matrix.')
end

if(isempty(A))
    error('MATLAB:pca:malformedInput',...
        'Input 1 must not be empty.')
end

% Set the inputs k, its, and l to default values, if necessary.
if(narg == 4)
    its=20;
else
    its=varargin{5};
    assertInt(its,5);
end

end

%% checkUSV
% check the type and sizes of the SVD decomposition U,S,V

function checkUSV(m,n,U,S,V)


if(~isnumeric(U))
  error('MATLAB:diffsnorm:malformedInput',...
        'Input 2 must be a floating-point matrix.')
end

if(~isnumeric(S))
  error('MATLAB:diffsnorm:malformedInput',...
        'Input 3 must be a floating-point matrix.')
end

[m2 k] = size(U);
[k2 l] = size(S);
[n2 l2] = size(V);

%
% Make sure that the dimensions of A, U, S, and V are commensurate.
%
if(m ~= m2)
  error('MATLAB:diffsnorm:malformedInput',...
        'The 1st dims. of Inputs 1 and 2 must be equal.')
end

if(k ~= k2)
  error('MATLAB:diffsnorm:malformedInput',...
        'The 2nd dim. of Input 2 must equal the 1st dim. of Input 3.')
end

if(l ~= l2)
  error('MATLAB:diffsnorm:malformedInput',...
        'The 2nd dims. of Inputs 3 and 4 must be equal.')
end

if(n ~= n2)
  error('MATLAB:diffsnorm:malformedInput',...
        'The 2nd dim. of Input 1 must equal the 1st dim. of Input 4.')
end

end


%% assertInt
% Check that the argument is an integer scalar
function assertInt(x,pos)
if ~isscalar(x) || ~isreal(x) || (floor(x)~=x)
    error('MATLAB:pca:malformedInput',...
        'Input %d must be an integer',pos);
end
end

%% T_tmp
function w=T_tmp(v,params)
global T_func
if isempty(params)
    w=feval(T_func,v);
else
    w=feval(T_func,params,v);
end
end

%% Tt_tmp
function w=Tt_tmp(v,params)
global Tt_func
if isempty(params)
    w=feval(Tt_func,v);
else
    w=feval(Tt_func,params,v);
end
end