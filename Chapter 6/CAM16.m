% CAM16
% Using algorithmic improvements suggested by Schlomer 2018: arXiv:1802.06067
%
% Function call:
% appearance = CAM16(xyz,xyzw)
%  or any combination of these name, value pairs:
% appearance = CAM16(xyz,xyzw,'adaptingLuminance',200,'whiteLuminance',1000,...
% 'relativeReferenceWhite',100,'relativeBackgroundLuminance',20,...
% 'Condition','average','D',1);
%
% Input Parameters:
% xyz: CIE XYZ values for stimuli (absolute or relative)
% xyzw: CIE XYZ values for white point (absolute or relative)
% Name, Value Pairs: 
% 'adaptingLuminance' (optional): Absolute luminance of background (La) in cd/m2, typically 20% of white luminance. 
% 'whiteLuminance' (optional): Absolute luminance of white (Lw) in cd/m2.
% 'relativeReferenceWhite' (optional): Relative luminance factor of reference white. 0-100. default: 100
% 'relativeBackgroundLuminance' (optional): Relative luminance factor of the background. 0-100. default: 20
% 'Condition' (optional): Viewing conditions. Can be 'dark', 'dim', or 'average'  (default: 'average')
% 'D' (optional): Degree of adaptation override. If unset, standard CAM16 D calculation is used.
%
% Future work: include hue quadrature.
%
% luke 2019-09-13

function appearance = CAM16(xyz,xyzw,varargin)

in = parseInputs(varargin{:}); 

if size(xyz,2)==3
    xyz=xyz';
end
if size(xyz,2) == size(xyz,1)
    warning('Square matrix assumed to be Nx3 orientation.')
end

if size(xyzw,1)==1
    xyzw=xyzw';
end


M16=[ 0.401288 0.650173 -0.051461;...
    -0.250268 1.204414 0.045854;...
    -0.002079 0.048952 0.953127];

if isempty(in.whiteLuminance)
    lw = xyzw(2);
else
    lw = in.whiteLuminance;
end
yw = in.relativeReferenceWhite; % default = 100
yb = in.relativeBackgroundLuminance; % default = 20
if isempty(in.adaptingLuminance)
    La = lw*yb/yw;
else
    La = in.adaptingLuminance;
end

switch in.Condition
    case 'dark'
        F=0.8;
        c=0.525;
        Nc=0.8;
    case 'dim'
        F=0.9;
        c=0.59;
        Nc=0.9;
    case 'average'
        F=1;
        c=0.69;
        Nc=1;
    otherwise
        error('Invalid condition')
end

%% Step 0: Calculate all Values/parameters which are independent of the input values

% Scale inputs to be luminance factor inputs (White Y = 100)
knorm = 100 / xyzw(2);
xyzw = xyzw * knorm;
xyz = xyz * knorm;

rgbw = M16*xyzw;
if isempty(in.D)
    D = F*(1-(1/3.6)*exp((-La-42)/92));
else
    D = in.D;
end

if D>1
    D = 1;
elseif D<0
    D = 0;
end

Dr = D*(xyzw(2)/rgbw(1))+1-D;
Dg = D*(xyzw(2)/rgbw(2))+1-D;
Db = D*(xyzw(2)/rgbw(3))+1-D;

% Schlomer has these equations for Dr, etc, which seem wrong:
% Dr = D*(xyzw(2)/rgbw(1))-1+D;
% Dg = D*(xyzw(2)/rgbw(2))-1+D;
% Db = D*(xyzw(2)/rgbw(3))-1+D;

k = 1/(5*La+1);
Fl = 0.2*k^4*(5*La)+0.1*(1-k^4)^2*(5*La)^(1/3);

n = yb/yw;
z = 1.48+(n)^0.5;
Nbb = 0.725*(1/n)^0.2;
Ncb=Nbb;

Rwc = Dr*rgbw(1);
Gwc = Dg*rgbw(2);
Bwc = Db*rgbw(3);

Raw = 400*((Fl*Rwc/100)^0.42/((Fl*Rwc/100)^0.42+27.13));
Gaw = 400*((Fl*Gwc/100)^0.42/((Fl*Gwc/100)^0.42+27.13));
Baw = 400*((Fl*Bwc/100)^0.42/((Fl*Bwc/100)^0.42+27.13));

Aw = (2*Raw+Gaw+Baw/20)*Nbb;

%% Step 1 Calculate Cone responses

RGB = M16*xyz;

%% Step 2 Color Adaptation
RGBc = [Dr;Dg;Db].*RGB;

%% Step 3 Cone compression
RGBap = 400.*sign(RGBc).*(.01*Fl.*abs(RGBc)).^.42./((.01*Fl.*abs(RGBc)).^.42+27.13);

%% Step 4 Calculate a, b, h, p2p, u

p2pabu = [2,1,1/20;...
    1,-12/11,1/11;...
    1/9,1/9,-2/9;...
    1,1,21/20]*RGBap;
p2p = p2pabu(1,:);
a = p2pabu(2,:);
b = p2pabu(3,:);
u = p2pabu(4,:);
h = mod(atan2d(b,a),360);
% h(h<0)=h(h<0)+360;

%% Step 5: Eccentricity and hue quadrature (skip)

et = .25*(cos(h*pi/180+2)+3.8);

%% Step 6: Achromatic response
A = p2p.*Nbb;

%% Step 7: Correlate of Lightness
J = 100*(A./Aw).^(c.*z);

%% Step 8: Correlate of Brightness
Q = (4./c).*sqrt(J/100).*(Aw+4).*Fl.^(.25);

%% Calculate chroma, colorfulness, saturation
t = ((50000/13).*Nc.*Ncb.*et.*sqrt(a.^2+b.^2))./(u+.305);
alpha = t.^0.9.*(1.64-0.29.^n).^0.73;
C = alpha.*sqrt(J/100);
M = C.*Fl.^.25;
s = 50.*sqrt((alpha.*c)./(Aw+4));

%% Check for imaginary colors
alpha = imag(J') > 0;
alpha = alpha | imag(Q') > 0;
alpha = alpha | imag(C') > 0;
alpha = alpha | imag(M') > 0;

J(alpha) = nan;
Q(alpha) = nan;
C(alpha) = nan;
M(alpha) = nan;
s(alpha) = nan;
h(alpha) = nan;

%% Format output

appearance.J = real(J');
appearance.Q = real(Q');
appearance.C = real(C');
appearance.M = real(M');
appearance.s = real(s');
appearance.h = real(h');

end


function options = parseInputs(varargin)
    persistent parser; 
    if isempty(parser)
        % Set up parser
        parser = inputParser();
        parser.FunctionName = mfilename;

        parser.addParameter('relativeReferenceWhite', 100, ...
            @(x) x>=0 && isreal(x));
        
        parser.addParameter('relativeBackgroundLuminance', 20, ...
            @(x) x>=0 && isreal(x));
        
        parser.addParameter('Condition', 'average', ...
            @(x) ischar(x));
        
        parser.addParameter('D', [], ...
            @(x) isreal(x) && x>=0 && x<=1);
        
        parser.addParameter('adaptingLuminance',[],@(x) isreal(x) && x>=0 )
        
        parser.addParameter('whiteLuminance',[],@(x) isreal(x) && x>=0 )
    end
    parser.parse(varargin{:})
    options = parser.Results;
end