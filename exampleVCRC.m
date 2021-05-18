clear;
clc;
close all;

%% Step 1. Load test displays. 
T = readtable('ExampleDisplay_N=7.xlsx');
[~, col] = size(T);
number_of_test_displays = (col - 1) ./ 3;
testDisplays = struct([]);

for i = 1:number_of_test_displays
    idx = (i - 1) .* 3 + 2;
    testDisplays(i).spd = table2array(T(1:end, idx:(idx+2))) .* 1;   
end

% We have 7 test displays.
did = 2; % display ID
spd = testDisplays(did).spd;

[volume, tetrahedrons] = computeVCRC(spd);

figure;
tetramesh(tetrahedrons.DT, tetrahedrons.QmAmBm,'FaceAlpha',0.1);
xlabel('a_M');
ylabel('b_M');
zlabel('Q');

title(['Volume: ' num2str(volume)]);