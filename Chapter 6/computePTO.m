clc
clear
close all

%% Step 1. Load the basic variables
T = readtable('CIED65_SPD.xlsx');

XYZtest = [47.5227963500000;50;54.4528875500000]; % A neutral color for evaluating OM (D65 with 50 cd/m2)
XYZwhite = XYZtest .* 4; % Reference white point assumed (D65 with 200 cd/m2)
wavelength = table2array(T(1:end, 1));   % The wavelength range = (390:1:780)nm
spd_D65 = table2array(T(1:end, 2)); % CIE D65 SPD
stdOb = table2array(T(1:end, 3:5)); % CIE 2015 10-deg standard observer
k = 683;

colors2 = flipud([0 0 0;
    1 0 0;
    0 1 0;
    1 0.647 0;
    0 0.7 0;
    0 1 1;
    0 0 1;
    0.6 0 1]);


%% Step 2. Load the categorical observers
T = readtable('CategoricalObservers.xlsx'); 
number_of_catObs = 10;
catObs = struct([]);

for i = 1:number_of_catObs
    idx = (i - 1) .* 3 + 2;
    catObs(i).CMFs = table2array(T(1:end, idx:(idx+2)));
end

%% Step 3. Load test displays. 
T = readtable('ExampleDisplay_N=7.xlsx');
[~, col] = size(T);
number_of_test_displays = (col - 1) ./ 3;
testDisplays = struct([]);

for i = 1:number_of_test_displays
    idx = (i - 1) .* 3 + 2;
    testDisplays(i).spd = table2array(T(1:end, idx:(idx+2))) .* 1;   
end

%% Step 4. Compute the colorfulness area (CA) of the displays
Y = 1000; % Whtie luminance (cd/m2)
white = 200;

for y = 1:number_of_test_displays
       
    M = k .* stdOb' * testDisplays(y).spd;        
    gain = Y ./ sum(M(2, :));
    testDisplays(y).spd = testDisplays(y).spd .* gain;
    M = k .* stdOb' * testDisplays(y).spd;        
    spd = testDisplays(y).spd;
    
    
%     Yb = 20;
    
%    cam(y).val(1) = CAM16(XYZ(:, 1), XYZwhite, 'whiteLuminance', Y, 'relativeBackgroundLuminance', Yb, 'Condition', 'dim');
%     cam(y).val(1) = CAM16(XYZ(:, 1), XYZwhite, 'Condition', 'dim');
%     cam(y).val(2) = CAM16(XYZ(:, 2), XYZwhite, 'Condition', 'dim');
%     cam(y).val(3) = CAM16(XYZ(:, 3), XYZwhite, 'Condition', 'dim');
    
    [CGV(y), tetrahedrons] = computeVCRC(spd);
end

%sliceLChJCh('LCh', 'D:\Yongmin\MATLAB\Chapter 6\data\', 'D:\Yongmin\MATLAB\Chapter 6\slice\' )

%% Step 5. Compute the OMMn of the displays
gain = (k .* stdOb' * spd_D65) \ XYZtest;
SPD_ref_std = spd_D65 .* gain;

OMMn = zeros(number_of_test_displays, number_of_catObs);
perc90 = zeros(number_of_test_displays, 1);

for y = 1:number_of_test_displays
    for x = 1:number_of_catObs
        M_test_ind = (k .* catObs(x).CMFs') * testDisplays(y).spd;
        XYZ_ref_ind = (k .* catObs(x).CMFs') * SPD_ref_std;
        
        RGB_test_ind = M_test_ind \ XYZ_ref_ind;
        SPD_test_ind = (testDisplays(y).spd' .* RGB_test_ind);
        
        OMMn(y, x) = computeOMMn2(SPD_ref_std, SPD_test_ind, stdOb, XYZtest');
    end    
    perc90(y) = percentilenthob(OMMn(y, :)', 0.90);
end

figure;
for i = 1:number_of_test_displays
    scatter(CGV(i), perc90(i), 30, colors2(i, :), 'filled'); hold on;
end

hold off
%xlim([10000 60000]);
ylim([5 16]);
xlabel('CGV');
ylabel('OMM_N');
legend('Rec.709', 'DCI.P3', 'Rec.2020 75%', 'Rec.2020 80%', 'Rec.2020 85%', ...
    'Rec.2020 90%', 'Rec.2020 95%', 'Rec.2020 100%', 'Location', 'northwest');

% figure;
% for i = 1:number_of_test_displays
%     scatter(CGV(i), perc90(i), 30, 'r', 'filled'); hold on;
% end
% 
% hold off
% %xlim([10000 60000]);
% ylim([5 16]);
% xlabel('CGV');
% ylabel('OMM_N');
% legend('Rec.709', 'DCI.P3', 'Rec.2020 75%', 'Rec.2020 80%', 'Rec.2020 85%', ...
%     'Rec.2020 90%', 'Rec.2020 95%', 'Rec.2020 100%', 'Location', 'northwest');







