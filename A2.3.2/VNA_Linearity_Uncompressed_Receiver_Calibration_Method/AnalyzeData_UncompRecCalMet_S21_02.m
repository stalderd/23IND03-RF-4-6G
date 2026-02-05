% Analyze Linearity Data - Uncompresse Receiver Calibration Method
% Daniel Stalder - 05.05.2025

clear all
close all
addpath('D:\svn_work\217.01 S-Parameter\Software\Metas.Vna\source\Matlab');
LoadVNATools();

%% Variables
vnaDevice = 'Agilent_PNA_N5227A_#1';
directory = 'Measurements_01';
name = 'StepAtt60dB(f-f)_SN123456_01';

method = 'UncompresseReceiverCalMethod';
outputfileNameSupplement = 'P1_01';

port1 = true; % true: P1 driving (S21), false = P2 driving (S12)
if port1
    aypy_Label = 'a1 p1';
    bxpy_Label = 'b2 p1';
    sxy_Label = 'S21';
else
    aypy_Label = 'a2 p2';
    bxpy_Label = 'b1 p2';
    sxy_Label = 'S12';
end

freq_GHz = [0.050 1 10 18 30 40 50]; % VNA frequencies to be analised
nFreq = length(freq_GHz);

nMeas = 100; % Number of sweep points per frequency point (CW sweep mode)
nMeasStartPosition = 1; % if source is unstable at the begining the start position can be adjusted to a higher index, where the source is stable
nMeasUsed = nMeas - nMeasStartPosition + 1;

stepAtt = 0:10:60; % Step attenuator positions to be analized (e.g. 0 dB to 60 dB in 10 dB steps)
nStepAtt = length(stepAtt);
stepAtt_uncompressed_Index = [3 3 3 3 3 3 3]; % Uncompressed test receiver measurement index per frequency (below compression influence and above noise floor influence), e.g. 20 dB position

power = -30:1:0; % VNA source power levels to be analized (e.g. -30 dBm to 0 dBm in 1 dB steps)
nPower = length(power);
nZero_points_a = 15;         % Number of zero points for the a-Receiver (lowest power levels points, below compression influence)
nZero_points_b_higher = 10;  % Number of zero points for the b-Receiver for the stepAtt measurements with higher b-Receiver power levels compared to stepAtt_uncompressed_Index (lowest power level points, below compression influence)
nZero_points_b_lower = 10;   % Number of zero points for the b-Receiver for the stepAtt measurements with lower b-Receiver power levels compared to stepAtt_uncompressed_Index (highest power level points, above noise floor influence)

plotLimdB = 0.01; % +/- min/max for the linearity plots

summaryFreqMinExt_GHz = .05;
summaryFreqMinNom_GHz = 1;
summaryFreqMaxNom_GHz = 40;
summaryFreqMaxExt_GHz = 50;

%% Load all data
aypy_all = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
sxy_all = zeros(nPower, nStepAtt, nFreq, nMeas);

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nFreq
            d = LoadVNADataAsStruct([pwd '\' directory '\' name '\' sprintf('%dHz\\%ddBm\\%s_%ddB.vdatb', freq_GHz(i3).*1e9, power(i1), name, stepAtt(i2))]);
            % 'S1,1', 'S1,2', 'S2,1', 'S2,2', 'a1_p1', 'a2_p2'
            if port1
                aypy_all(i1, i2, i3, :) = d.Data(:, 5);
                sxy_all(i1, i2, i3, :) = d.Data(:, 3);
            else
                aypy_all(i1, i2, i3, :) = d.Data(:, 6);
                sxy_all(i1, i2, i3, :) = d.Data(:, 2);
            end
        end
    end
end

bxpy_all = aypy_all .* sxy_all;

%% Mean data from samples
aypy = LinProp(zeros(nPower, nStepAtt, nFreq)); % V/sqrt(kOhm)
bxpy = LinProp(zeros(nPower, nStepAtt, nFreq)); % V/sqrt(kOhm)
sxy = LinProp(zeros(nPower, nStepAtt, nFreq));

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nFreq
            info = sprintf('Power %d, Step Att %d dB', power(i1), stepAtt(i2));
            tpy = LinProp([reshape(real(aypy_all(i1, i2, i3, :)), [nMeasUsed, 1]) reshape(imag(aypy_all(i1, i2, i3, :)), [nMeasUsed, 1]) ...
                           reshape(real(bxpy_all(i1, i2, i3, :)), [nMeasUsed, 1]) reshape(imag(bxpy_all(i1, i2, i3, :)), [nMeasUsed, 1]) ...
                           reshape(real(sxy_all(i1, i2, i3, :)), [nMeasUsed, 1]) reshape(imag(sxy_all(i1, i2, i3, :)), [nMeasUsed, 1])], 'samples', [sxy_Label ', ' info]);
            aypy(i1, i2, i3) = abs(tpy(1) + 1i.*tpy(2));
            bxpy(i1, i2, i3) = abs(tpy(3) + 1i.*tpy(4));
            sxy(i1, i2, i3) = abs(tpy(5) + 1i.*tpy(6));
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
cal_aypy_dBm = LinProp(zeros(nPower, nFreq));
aypy_rel_lin = LinProp(zeros(nPower, nFreq));
sxy_cor_aypy_rel_lin = LinProp(zeros(nPower, nStepAtt, nFreq));
bxpy_rel_lin = LinProp(zeros(nPower, nStepAtt, nFreq));
for i3 = 1:nFreq
    cal_aypy_dBm(:,i3) = squeeze(aypydBm(:,stepAtt_uncompressed_Index(i3),i3));
    sxy_ref_aypy = 0;
    for ii = 1:nZero_points_a
        sxy_ref_aypy = sxy_ref_aypy + sxy_abs(ii,stepAtt_uncompressed_Index(i3),i3) ./ nZero_points_a; % sxy -> b/a
    end
    aypy_rel_lin(:,i3) = sxy_ref_aypy ./ sxy_abs(:, stepAtt_uncompressed_Index(i3), i3); % a_lin = (b_uncompressed_ref / a_ref) * (a / b_uncompressed)
    
    sxy_cor_aypy_rel_lin(:,:,i3) = sxy_abs(:,:,i3) .* aypy_rel_lin(:,i3); % Correct the data for the nonlinearity of the a receiver: sxy_cor_aypy_rel_lin = b / a * a_Lin
    
    for i2 = 1:stepAtt_uncompressed_Index(i3) % Measurements with b-receiver power levels ≥ compared to stepAtt_uncompressed_Index
        sxy_ref_bxpy_higher = 0;
        for ii = 1:nZero_points_b_higher
            sxy_ref_bxpy_higher = sxy_ref_bxpy_higher + sxy_cor_aypy_rel_lin(ii,i2,i3) ./ nZero_points_b_higher;
        end
        bxpy_rel_lin(:,i2,i3) = sxy_cor_aypy_rel_lin(:,i2,i3) ./ sxy_ref_bxpy_higher; % b_lin = (b / a_corrected) * (a_corrected_ref / b_ref)
    end
    
    for i2 = (stepAtt_uncompressed_Index(i3)+1):nStepAtt  % Measurements with b-receiver power levels < compared to stepAtt_uncompressed_Index
        sxy_ref_bxpy_lower = 0;
        for ii = nPower-nZero_points_b_lower + 1:nPower
            sxy_ref_bxpy_lower = sxy_ref_bxpy_lower + sxy_cor_aypy_rel_lin(ii,i2,i3) ./ nZero_points_b_lower;
        end
        bxpy_rel_lin(:,i2,i3) = sxy_cor_aypy_rel_lin(:,i2,i3) ./ sxy_ref_bxpy_lower; % b_lin = (b / a_corrected) * (a_corrected_ref / b_ref)
    end
end


%% Plots S-Param and normalized S-Param for the different Step Attenuator States
mean_sxy_dB = mean(get_value(sxydB), 1); % Calculate the mean of the data along the 1 dimension (DUT)
sxy_dB_norm = sxydB - mean_sxy_dB; % Normalize the data by subtracting the mean

for i2 = 1:nStepAtt
    h0 = figure();
    subplot(2,1,1);
    plotu(freq_GHz, squeeze(sxydB(:,i2,:)));
    xlabel('Frequency / GHz');
    ylabel(sprintf('|%s| / dB', sxy_Label));
    grid on;

    subplot(2,1,2);
    plotu(freq_GHz, squeeze(sxy_dB_norm(:,i2,:)));
    xlabel('Frequency / GHz');
    ylabel(sprintf('|%s| - mean(|%s|) / dB', sxy_Label, sxy_Label));
    grid on;

    a = axes;
    t1 = title({strrep([vnaDevice ', ' name ' @ ' sprintf('%g', stepAtt(i2)) ' dB'], '_', ' '), ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h0, [vnaDevice '_' name '_S-Param_Raw_' outputfileNameSupplement]);
end
ps2pdf([vnaDevice '_' name '_S-Param_Raw_' outputfileNameSupplement])

%% Plots per Frequency
aypydBm_min = min(double(aypydBm(:)));
aypydBm_max = max(double(aypydBm(:)));
bxpydBm_min = min(double(bxpydBm(:)));
bxpydBm_max = max(double(bxpydBm(:)));

for i3 = 1:nFreq
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
    t1 = title({strrep([vnaDevice ', ' name ' @ ' sprintf('%g', freq_GHz(i3)) ' GHz'], '_', ' '), method, ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h1, [vnaDevice '_' name '_' method '_' outputfileNameSupplement]);
end
ps2pdf([vnaDevice '_' name '_' method '_' outputfileNameSupplement])

%% Plots Summary
ii_Ext = (find(freq_GHz == summaryFreqMinExt_GHz):find(freq_GHz == summaryFreqMaxExt_GHz))';
ii_Nom = (find(freq_GHz == summaryFreqMinNom_GHz):find(freq_GHz == summaryFreqMaxNom_GHz))';

hsum = figure();
subplot(2,1,1);
for i3 = ii_Ext
    plotv(cal_aypy_dBm(:,i3), 20.*log10(aypy_rel_lin(:,i3)));
    hold on;
end
xlabel([aypy_Label ' / dBm']);
ylabel(['Rel. Non Linearity ' aypy_Label ' / dB']);
xlim([-35 10]);
ylim([-plotLimdB plotLimdB]);
grid on;
hold off;

subplot(2,1,2);
for i2 = 1:nStepAtt
    for i3 = ii_Ext
        plotv(squeeze(bxpydBm(:,i2,i3)), 20.*log10(squeeze(bxpy_rel_lin(:,i2,i3))));
        hold on;
    end
end
xlabel([bxpy_Label ' / dBm']);
ylabel(['Rel. Non Linearity ' bxpy_Label ' / dB']);
xlim([-60 10]);
ylim([-plotLimdB plotLimdB]);
grid on;
hold off;

a = axes;
t1 = title({strrep([vnaDevice ', ' name  ' @ ' sprintf('%g', summaryFreqMinExt_GHz) ' GHz to ' sprintf('%g', summaryFreqMaxExt_GHz) ' GHz'], '_', ' '), method, ''});
a.Visible = 'off';
t1.Visible = 'on';

figure2ps(hsum, [vnaDevice '_' name '_' method 'Summary_' outputfileNameSupplement]);

hsum = figure();
subplot(2,1,1);
for i3 = ii_Nom
    plotv(cal_aypy_dBm(:,i3), 20.*log10(aypy_rel_lin(:,i3)));
    hold on;
end
xlabel([aypy_Label ' / dBm']);
ylabel(['Rel. Non Linearity ' aypy_Label ' / dB']);
xlim([-35 10]);
ylim([-plotLimdB plotLimdB]);
grid on;
hold off;

subplot(2,1,2);
for i2 = 1:nStepAtt
    for i3 = ii_Nom
        plotv(squeeze(bxpydBm(:,i2,i3)), 20.*log10(squeeze(bxpy_rel_lin(:,i2,i3))));
        hold on;
    end
end
xlabel([bxpy_Label ' / dBm']);
ylabel(['Rel. Non Linearity ' bxpy_Label ' / dB']);
xlim([-60 10]);
ylim([-plotLimdB plotLimdB]);
grid on;
hold off;

a = axes;
t1 = title({strrep([vnaDevice ', ' name  ' @ ' sprintf('%g', summaryFreqMinNom_GHz) ' GHz to ' sprintf('%g', summaryFreqMaxNom_GHz) ' GHz'], '_', ' '), method, ''});
a.Visible = 'off';
t1.Visible = 'on';

figure2ps(hsum, [vnaDevice '_' name '_' method 'Summary_' outputfileNameSupplement]);
ps2pdf([vnaDevice '_' name '_' method 'Summary_' outputfileNameSupplement])

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
