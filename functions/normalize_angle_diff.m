
function diff = normalize_angle_diff(a, b)
% Calculates the shortest angular difference between two angles (a - b).
% The result is normalized to the range [-pi, pi].
%
% INPUTS:
%   a: First angle (or vector of angles) in radians.
%   b: Second angle (or vector of angles) in radians.
%
% OUTPUT:
%   diff: The normalized angular difference.

    diff = atan2(sin(a - b), cos(a - b));
end