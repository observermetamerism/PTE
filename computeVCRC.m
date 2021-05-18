function [volume, tetrahedrons] = computeVCRC(spd)
    %% Step 2. Load the basic variables
    T = readtable('CIED65_SPD.xlsx');
    stdOb = table2array(T(1:end, 3:5)); % CIE 2015 10-deg standard observer
    k = 683;
    M = (k .* stdOb' * spd);

    %% Make a primary + secondary
    [~, col] = size(M);

    Me = zeros(3, col .* 2);
    Me(:, 1) = M(:, 1);
    k = 2;
    for i = 1:col-1
        Me(:, k) = M(:, i) + M(:, i + 1);
        k = k + 1;
        Me(:, k) = M(:, i + 1);
        k = k + 1;
    end
    Me(:, end) = M(:, 1) + M(:, end);

    %% Assuming 10-bit
    lvs = (0:64:1023);
    lvs = [lvs, 1023] ./ 1023;
    number_of_lvs = length(lvs);
    [~, number_of_primaries] = size(Me);

    for y = 1:number_of_primaries
        primaries(y).XYZs = zeros(number_of_lvs, 3);

        for x = 1:number_of_lvs
            primaries(y).XYZs(x, :) = (lvs(x).^2.2) .* Me(:, y)';
        end
    end

    primaries(number_of_primaries + 1).XYZs = sum(M');

    QMh = zeros(number_of_primaries .* number_of_lvs + 1, 3);

    XYZtest = [47.5227963500000;50;54.4528875500000]; % A neutral color for evaluating OM (D65 with 50 cd/m2)
    XYZwhite = XYZtest .* 4; % Reference white point assumed (D65 with 200 cd/m2)
    k = 1;
    for y = 1:number_of_primaries
        for x = 1:number_of_lvs
            cam = CAM16(primaries(y).XYZs(x, :)', XYZwhite, 'Condition', 'dim');
            QMh(k, :) = [cam.Q, computeMp(cam), cam.h];

            k = k + 1;
        end
    end

    cam = CAM16(primaries(number_of_primaries + 1).XYZs', XYZwhite, 'Condition', 'dim');
    QMh(k, :) = [cam.Q, computeMp(cam), cam.h];

    QMaMb = computeMaMb(QMh(1:end-1, :));
    QpMaMb = QMaMb(17:17:end, :);
    %
    % figure;
    % scatter3(QMaMb(:, 2), QMaMb(:, 3), QMaMb(:, 1), 30, 'k', 'filled');

    % Step 4-1
    DT1 = delaunay(QMaMb(:, 2), QMaMb(:, 3));
    % figure;
    % triplot(DT1, QMaMb(:, 2), QMaMb(:, 3));

    % Step 4-2
    QwMaMb = computeMaMb(QMh(end, :));
    QtMaMb = [QpMaMb; QwMaMb];
    DT2 = delaunay(QtMaMb(:, 2), QtMaMb(:, 3));
    % figure;
    % triplot(DT2, QtMaMb(:, 2), QtMaMb(:, 3));

    % An inner point
    QiMaMb = mean([QtMaMb; 0 0 0]);

    QaMaMb = [QMaMb; QtMaMb; QiMaMb];
    DT3 = DT2 + length(QMaMb);
    DT_all = [DT1; DT3];
    DT_all = [DT_all, repmat(length(QaMaMb), length(DT_all), 1)];

    [row, ~] = size(DT_all);

    volume = 0;

    for x = 1:row
        P1 = QaMaMb(DT_all(x, 1), :);
        P2 = QaMaMb(DT_all(x, 2), :);
        P3 = QaMaMb(DT_all(x, 3), :);
        P4 = QaMaMb(DT_all(x, 4), :);

        volume = volume + (1/6*abs(dot(cross(P2-P1,P3-P1),P4-P1)));
    end

    QzMaMb = [QaMaMb(:, 2), QaMaMb(:, 3), QaMaMb(:, 1)];
    
    tetrahedrons.QmAmBm = QzMaMb;
    tetrahedrons.DT = DT_all;
    
%     figure;
%     tetramesh(DT_all, QzMaMb,'FaceAlpha',0.1);
%     xlabel('a_M');
%     ylabel('b_M');
%     zlabel('Q');
% 
%     disp(volume);
end

function QMaMb = computeMaMb(QMh)
    [row, ~] = size(QMh);
    
    QMaMb = zeros(row, 3);
    
    for x = 1:row
        h = QMh(x, 3);
        M = QMh(x, 2);
        QMaMb(x, :) = [QMh(x, 1), M .* cos(deg2rad(h)), M .* sin(deg2rad(h))];
    end
%     QMhmax = QMaMb(17:17:row, :);
%     QMhmax = [QMhmax; QMh(1, :); QMh(row, :)];
%     
%     QMaMb(row + 1, :) = mean(QMhmax);
end

function Mp = computeMp(cam)
    M = cam.M;

    p(1) = 0.0015;
    p(2) = 0.7183;
    p(3) = 0;
    Mp = polyval(p, M);
end


