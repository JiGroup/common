function [LR_out] = interval_merge(LR_in)
%INTERVAL_MERGE - Merge real number intervals from endpoint data
%
% Syntax:  [LR_out] = function_name(LR_in)
% Inputs:
%    LR_in  - Nx2 matrix of endpoints. Left endpoints in column 1 and right
%               endpoints in column 2.
% Outputs:
%    LR_out  - Nx2 matrix of merged endpoints.
% Example: 
%    LR = interval_merge([L R]);
% Other m-files required: none
% Subfunctions: none
% MAT-files required: none
% See also: 

% Author: Pat Flaherty
% Work address
% email: pflahert@stanford.edu
% Website: http://www.
% Jan 2009; Last revision: 09-Jan-2009

N = length(LR_in(:,1));

% Mark left and right endpoints {+1,-1} and sort
A = [LR_in(:,1), repmat(1,N,1); LR_in(:,2), repmat(-1,N,1)];
A = sortrows(A,1);

% Cumulative sum and threshold marks
A(:,2) = cumsum(double(A(:,2)))>0;

% Increment through endpoints and record start/stop
LR_out = []; 
in_flag = 0; % flag if we're inside and interval
for i = 1:length(A(:,1))
    if(in_flag == 0 && A(i,2) == 1) % interval start
        L_point = A(i,1);
        in_flag = 1;
    elseif(in_flag == 1 && A(i,2) == 0) % interval stop
        LR_out = [LR_out; L_point, A(i,1)];
        in_flag = 0;
    end % else not on interval boundary
end
