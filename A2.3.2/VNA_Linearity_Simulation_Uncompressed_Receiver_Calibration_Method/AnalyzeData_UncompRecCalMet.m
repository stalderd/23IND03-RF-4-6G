% Analyze Linearity Data - Uncompresse Receiver Calibration Method
% Daniel Stalder - 05.05.2025 - 13.08.2025

clear all
close all

%% Variables
vnaDevice = 'Simulation';
name = 'StepAtt';
method = 'UncompresseReceiverCalMethod';

port1 = true; % true = P1, false = P2

if port1
    aypy_Label = 'a1 p1';
    bxpy_Label = 'b2 p1';
    sxy_Label = 'S21';
else
    aypy_Label = 'a2 p2';
    bxpy_Label = 'b1 p2';
    sxy_Label = 'S12';
end

dataSets = {'LinTestDataSet1' , 'LinTestDataSet2', 'LinTestDataSet3', 'LinTestDataSet4', 'LinTestDataSet5', 'LinTestDataSet6', 'LinTestDataSet7'};
nDataSets = length(dataSets);
nMeas = 100;

stepAtt = 0:10:70;
nStepAtt = length(stepAtt);
stepAtt_uncompressed_Index = 5; % Uncompressed test receiver measurement index (below compression influence and above noise floor influence)

power = -30:1:0;
nPower = length(power);
nZero_points_a = 5;         % Number of zero points for the a-Receiver (lowest power levels points, below compression influence)
nZero_points_b_higher = 5;  % Number of zero points for the b-Receiver for the stepAtt measurements with higher b-Receiver power levels compared to stepAtt_uncompressed_Index (lowest power level points, below compression influence)
nZero_points_b_lower = 5;   % Number of zero points for the b-Receiver for the stepAtt measurements with lower b-Receiver power levels compared to stepAtt_uncompressed_Index (highest power level points, above noise floor influence)

plotLimdB = 0.1; % +/- min/max for the linearity plots

%% Load all data
a1p1 = zeros(nPower, nStepAtt, nDataSets, nMeas); % sqrt(mW)
b2p1 = zeros(nPower, nStepAtt, nDataSets, nMeas); % sqrt(mW)
s21 = zeros(nPower, nStepAtt, nDataSets, nMeas);

for i3 = 1:nDataSets
    load([dataSets{i3} '.mat'])
    for i1 = 1:nPower
        for i2 = 1:nStepAtt
            a1p1(i1, i2, i3, :) = a1(:, i1, i2);
            b2p1(i1, i2, i3, :) = b1(:, i1, i2);
            s21(i1, i2, i3, :) = b1(:, i1, i2)./a1(:, i1, i2);
        end
    end
end

%% Mean data from samples
aypy = LinProp(zeros(nPower, nStepAtt, nDataSets)); % sqrt(mW)
bxpy = LinProp(zeros(nPower, nStepAtt, nDataSets)); % sqrt(mW)
sxy = LinProp(zeros(nPower, nStepAtt, nDataSets)); % sqrt(mW)

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nDataSets
            info = sprintf('Power %d, Step Att %d dB', power(i1), stepAtt(i2));
            tp1 = LinProp([reshape(real(a1p1(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(a1p1(i1, i2, i3, :)), [nMeas, 1]) ...
                           reshape(real(b2p1(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(b2p1(i1, i2, i3, :)), [nMeas, 1]) ...
                           reshape(real(s21(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(s21(i1, i2, i3, :)), [nMeas, 1])], 'samples', ['S21, ' info]);
            aypy(i1, i2, i3) = abs(tp1(1) + 1i.*tp1(2));
            bxpy(i1, i2, i3) = abs(tp1(3) + 1i.*tp1(4));
            sxy(i1, i2, i3) = abs(tp1(5) + 1i.*tp1(6));
        end
    end
end

sxy_abs = abs(sxy);

%% Mean data in dBm
aypydBm = 20.*log10(aypy);
bxpydBm = 20.*log10(bxpy);
sxydB = 20.*log10(abs(sxy));

%% Assume ideal b Receiver at specified step attenuator position
%% Calibrate a Receiver (a Receiver is always in the same level range)
cal_aypy_dBm = squeeze(aypydBm(:,stepAtt_uncompressed_Index,:));

sxy_ref_aypy = LinProp(zeros(1, 1, nDataSets));
for ii = 1:nZero_points_a
    sxy_ref_aypy = sxy_ref_aypy + sxy_abs(ii,stepAtt_uncompressed_Index,:) ./ nZero_points_a; % sxy -> b/a
end
aypy_rel_lin = sxy_ref_aypy ./ sxy_abs(:, stepAtt_uncompressed_Index, :); % a_lin = (b_uncompressed_ref / a_ref) * (a / b_uncompressed)

sxy_cor_aypy_rel_lin = sxy_abs .* aypy_rel_lin; % Correct the data for the nonlinearity of the a receiver: sxy_cor_aypy_rel_lin = b / a * a_Lin

bxpy_rel_lin=LinProp(zeros(nPower, nStepAtt, nDataSets)); % Uncompressed test receiver measurement index
bxpy_rel_lin(:,stepAtt_uncompressed_Index,:) = sxy_cor_aypy_rel_lin(:,stepAtt_uncompressed_Index,:) ./ sxy_cor_aypy_rel_lin(:,stepAtt_uncompressed_Index,:);

for i2 = 1:stepAtt_uncompressed_Index-1 % Measurements with higher b-receiver power levels compared to stepAtt_uncompressed_Index
    sxy_ref_bxpy_higher = 0;
    for ii = 1:nZero_points_b_higher
        sxy_ref_bxpy_higher = sxy_ref_bxpy_higher + sxy_cor_aypy_rel_lin(ii,i2,:) ./ nZero_points_b_higher;
    end
    bxpy_rel_lin(:,i2,:) = sxy_cor_aypy_rel_lin(:,i2,:) ./ sxy_ref_bxpy_higher; % b_lin = (b / a_corrected) * (a_corrected_ref / b_ref)
end

for i2 = (stepAtt_uncompressed_Index+1):nStepAtt  % Measurements with lower b-receiver power levels compared to stepAtt_uncompressed_Index
    sxy_ref_bxpy_lower = 0;
    for ii = nPower-nZero_points_b_lower + 1:nPower
        sxy_ref_bxpy_lower = sxy_ref_bxpy_lower + sxy_cor_aypy_rel_lin(ii,i2,:) ./ nZero_points_b_lower;
    end
    bxpy_rel_lin(:,i2,:) = sxy_cor_aypy_rel_lin(:,i2,:) ./ sxy_ref_bxpy_lower; % b_lin = (b / a_corrected) * (a_corrected_ref / b_ref)
end

%% Plots S-Param and normalized S-Param for the different DataSets
mean_sxy_dB = mean(get_value(sxydB), 1); % Calculate the mean of the data along the 1 dimension (DUT)
sxy_dB_norm = sxydB - mean_sxy_dB; % Normalize the data by subtracting the mean

for i3 = 1:nDataSets
    h0 = figure();
    subplot(2,1,1);
    plotu(stepAtt, squeeze(sxydB(:,:,i3)));
    xlabel('StepAtt / dB');
    ylabel(sprintf('|%s| / dB', sxy_Label));
    grid on;

    subplot(2,1,2);
    plotu(stepAtt, squeeze(sxy_dB_norm(:,:,i3)));
    xlabel('StepAtt / dB');
    ylabel(sprintf('|%s| - mean(|%s|) / dB', sxy_Label, sxy_Label));
    grid on;

    a = axes;
    t1 = title({strrep([vnaDevice ', ' name ' @ ' dataSets{i3}], '_', ' '), ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h0, [vnaDevice '_' name '_S-Param_Raw']);
end
ps2pdf([vnaDevice '_' name '_S-Param_Raw'])

%% Plots per DataSet
aypydBm_min = min(double(aypydBm(:)));
aypydBm_max = max(double(aypydBm(:)));
bxpydBm_min = min(double(bxpydBm(:)));
bxpydBm_max = max(double(bxpydBm(:)));

for i3 = 1:nDataSets
    h1 = figure();
    subplot(2,2,1);
    plotv(aypy(:,:,i3), bxpy(:,:,i3));
    xlabel([aypy_Label ' / (V/sqrt(kOhm))']);
    ylabel([bxpy_Label ' / (V/sqrt(kOhm))']);
    grid on;

    subplot(2,2,2);
    plotv(aypydBm(:,:,i3), bxpydBm(:,:,i3));
    xlabel([aypy_Label ' / dBm']);
    ylabel([bxpy_Label ' / dBm']);
    xlim([aypydBm_min aypydBm_max]);
    ylim([bxpydBm_min bxpydBm_max]);
    grid on;

    subplot(2,2,3);
    plotu(cal_aypy_dBm(:,i3), 20.*log10(squeeze(aypy_rel_lin(:,i3))));
    xlabel([aypy_Label ' / dBm']);
    ylabel(['Rel. Non Linearity ' aypy_Label ' / dB']);
    xlim([aypydBm_min aypydBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(2,2,4);
    plotu(bxpydBm(:,:,i3), 20.*log10(bxpy_rel_lin(:,:,i3)));
    xlabel([bxpy_Label ' / dBm']);
    ylabel(['Rel. Non Linearity ' bxpy_Label ' / dB']);
    xlim([bxpydBm_min bxpydBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    a = axes;
    t1 = title({strrep([vnaDevice ', ' name ' @ ' dataSets{i3} ', ' method], '_', ' '), ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h1, [vnaDevice '_' name '_' method]);
end
ps2pdf([vnaDevice '_' name '_' method])

%% Plot Subfunctions
function plotv(x, y) 
plot(double(x), double(y));
end

function plotu(x, y)
k = 2;
x2 = LinProp(x);
y2 = LinProp(y);
x_value = get_value(x2);
y_value = get_value(y2);
y_std = get_stdunc(y2);
y_pos = y_value + k.*y_std;
y_neg = y_value - k.*y_std;

plot(x_value, y_value);
hold on
set(gca,'ColorOrderIndex',1)
plot(x_value, y_pos, 'LineStyle', '-.');
set(gca,'ColorOrderIndex',1)
plot(x_value, y_neg, 'LineStyle', '-.');
hold off
end

function figure2ps(handle, filename)
set(handle,'PaperUnits','Centimeters');
set(handle,'PaperPosition',[0.63 0.63 19.73 28.43]);
set(handle,'PaperSize',[21.00 29.70]);
pause(0.5);
print(handle, '-dpsc', '-r300', '-append', [filename '_temp']);
end

function ps2pdf(filename)
system(['ps2pdf -dEPSCrop -dPDFX ' filename '_temp.ps ' filename '.pdf'],'-echo');
system(['del ' filename '_temp.ps'],'-echo');
end