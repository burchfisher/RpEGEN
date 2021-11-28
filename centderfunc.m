function [dist,xl,yl] = centderfunc(I)

% Created 11/3/21 by G. Burch Fisher

% VARIABLES
% I = BW mask created from the ROIxy data located in data.ROIBW from the
    % RPEGen script.
    % example: I = data.ROIBW{3,1};

%% STEP ONE - CENTERLINE CREATION AND PATH DISTANCE CALCULATION

% Thin I using the Image Processing toolbox to create a skeleton
It = bwmorph(I,'thin','inf');               %Creates skeleton
[y,x] = find(bwmorph(It,'endpoints'));      %Extracts x,y coords for all endpoints

% Find the start endpoint coords - should always be in the top left quadrant and
    % the closest to the x axis
idy = find(y < size(I,1)/2);
idstart = find(x == min(x(idy)));
xs = x(idstart); 
ys = y(idstart);

% Find the distance to each pixel in the skeleton from the start point
D = bwdistgeodesic(It,xs,ys,'quasi-euclidean');

% Extract location of the farthest endpoint
[ye,xe] = find(D == max(max(D))); 

clear I idstart idy It n x y xs ys

%% STEP TWO - REMOVING SPURS

% For loop to remove all spurs along the longest path from the start point
    % to the endpoint

% Works by starting at the farthest endpoint (highest distance value from 
    % the start from bwdistgeodesic) and then goes to the lowest point 
    % in a 3x3 kernal with the point of observation in the middle 
    % until it reaches the start point where distance = 0

val = D(ye,xe);
x = xe;
y = ye;

for t = 1:(size(D,1)^2);
    if t==1;               % To record the start at the end point
        xl(t) = x;
        yl(t) = y;
        dist(t) = val;
    else
        if val > 0;        % When val = 0 you have reached the start endpoint 
            p = x-1:x+1; q = y-1:y+1;   % in the upper left part of the image
            [b,a] = find(D(q,p) == min(D(q,p),[], 'all'));
            x = p(a); y = q(b);
            xl(t) = x; yl(t) = y;
            val = D(y,x);
            dist(t) = val;
        else
            return
        end
    end
end

clear a b p q posa posb t val vala valb x xe y ye 

end



