function [x_rotated, y_rotated] = rot2D(x, y, theta, x_center, y_center)

% create a matrix input points
v = [x;y];

%2D Rotation Matrix
R = [cos(theta) -sin(theta); sin(theta) cos(theta)];


% create a matrix for center
center = repmat([x_center; y_center], 1, length(x));

% this can be done in one line as:
vo = R*(v - center) + center;
% pick out the vectors of rotated x- and y-data
x_rotated = vo(1,:);
y_rotated = vo(2,:);