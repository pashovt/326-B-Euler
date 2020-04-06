% Beam properties
YoungsModulus = 68.9*10.^3; % GPa
PoissonRatio = 0.3;
Density = 2720; % kg/m.^3

% To metric units
ToMeter = 1000;

% Segment lengths
L1 = 139/ToMeter;
L2 = 88.41/ToMeter;
L3 = 133/ToMeter;
L4 = 88.41/ToMeter;
L5 = 203/ToMeter;
% Segment 2 and 4 tilt angle
degreeOfTilt = 43.62;

% Beam length
Length = L1 + L2 + L3 + L4 + L5;
ShapeLength = L1 + L3 + L5 + 2*cosd(degreeOfTilt)*L2;

% ZY view parameters
thickness = 4.8;
height = 25.4;
topLength = 50.8;

% Areas of sections
Side1 = thickness*height;
Side2 = thickness*height;
Top = (topLength-2*thickness)*thickness;

% Distance to centroid of segment
ToCentroidSide1 = height/2;
ToCentroidSide2 = height/2;
ToCentroidTop = height - thickness + thickness/2;

% Distance to centroid of shape
TotalArea = Side1 + Side2 + Top;
% Centroid location
YA = (Side1*ToCentroidSide1) + (Side2*ToCentroidSide2) + (Top*ToCentroidTop);
Centroid = YA/TotalArea;

% Second moment of Inertia
Ix = ((1/12)*thickness*height.^3) + Side1*abs(ToCentroidSide1-Centroid).^2 + ...
    ((1/12)*(topLength-2*thickness)*thickness.^3) + Top*abs(ToCentroidTop-Centroid).^2 + ...
    ((1/12)*thickness*height.^3) + Side2*abs(ToCentroidSide2-Centroid).^2;

% K values - circular eigen frequency
syms B1
syms B2
syms B3
syms B4
syms k
syms x
syms W(x)
% free - free beam => BC : w'''=0 & w''=0
C1 = B1*cosh(k*x);
C2 = B2*sinh(k*x);
C3 = B3*cos(k*x);
C4 = B4*sin(k*x);

% Equation
W(x) = C1 + C2 + C3 + C4; % B1*cosh(k*x) + B2*sinh(k*x) + B3*cosh(k*x) + B4*sinh(k*x);

% Boundary conditions
FirstCondition = diff(W, x, 2);
SecondCondition = diff(W, x, 3);
BC1 = FirstCondition(0);
BC2 = SecondCondition(0);
BC3 = FirstCondition(ShapeLength);
BC4 = SecondCondition(ShapeLength);

% BC1 = B1*k^2 - B3*k^2 => B1 - B3 = 0 => B1 = B3
% BC2 = B2*k^3 - B4*k^3 => B2 - B4 = 0 => B2 = B4
% BC3 = B1*k^2*cosh((32591*k)/50000) - B3*k^2*cos((32591*k)/50000) - B4*k^2*sin((32591*k)/50000) + B2*k^2*sinh((32591*k)/50000)
% BC4 = B2*k^3*cosh((32591*k)/50000) - B4*k^3*cos((32591*k)/50000) + B3*k^3*sin((32591*k)/50000) + B1*k^3*sinh((32591*k)/50000)

% Using BC1 & BC2, BC3 and BC4 can be rearranged in the following metrix
% method. Where BC3/k.^2 and BC4/k.^3
% BC3
    % B1*cosh((32591*k)/50000) - B1*cos((32591*k)/50000) + ...
    % B2*sinh((32591*k)/50000) - B2*sin((32591*k)/50000) 
% BC4
    % B1*sin((32591*k)/50000) + B1*sinh((32591*k)/50000) + ...
    % B2*cosh((32591*k)/50000) - B2*cos((32591*k)/50000)
%
% | cosh(k*Length) - cos(k*Length)  sinh(k*Length) - sin(k*Length)| |B1|
% | sin(k*Length) + sinh(k*Length)  cosh(k*Length) - cos(k*Length)| |B2|
%
% (cosh(k*Length) - cos(k*Length))*(cosh(k*Length) - cos(k*Length))
% cosh.^2 + cos.^2 - 2cosh*cos
% (sinh(k*Length) - sin(k*Length))*(sin(k*Length) + sinh(k*Length))
% sinh.^2 - sin.^2
% cosh.^2 + cos.^2 - 2cosh*cos - sinh.^2 + sin.^2
% Where
% cos.^2 + sin.^2 = 1
% cosh.^2 - sinh.^2 = 1
%
% 1 - 2cosh*cos + 1 = 0
% - 2cosh*cos = - 2
% Characteristic equation
% cosh*cos = 1
% 
% n = 0:5;
% kL = (n*pi)/(2);

kL(1) = 4.694;
kL(2) = 7.855;
kL(3) = 10.996;
kL(4) = 14.137;
kL(5) = 17.279;

k = kL/ShapeLength

BC1
BC2
BC3
BC4

% Euler formula
omega = sqrt(((YoungsModulus*Ix)/(Density*TotalArea))*k.^4)
NaturalFrequency = omega/(2*pi)
