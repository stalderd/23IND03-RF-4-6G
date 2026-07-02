% Analyze Linearity Data - Quadratic Model
% Michael Wollensack METAS - 23.12.2021 - 05.05.2025
% Daniel Stalder METAS - 02.07.2026

clear all
close all
addpath('C:\Users\Public\Documents\Metas.Vna.Tools\Matlab');
LoadVNATools();

%% Variables
vnaDevice = 'Keysight_PNA_N5227B_#2_VDI_WR10_2';
directory = 'Measurements_01';

method = 'opt'; % Quadratic Model (Optimaization)
outputfileNameSupplement = 'P1_8dB'; % P1 -> port 1 is driving, change to P2 if P2 is driving

port1 = true; % true: P1 driving (a1_p1 and b2_p1), false = P2 driving (a2_p2 and b1_p2)
if port1
	twoPortDirStartString = 'P1_';
    aypy_Label = 'a1 p1';
    bxpy_Label = 'b2 p1';
    sxy_Label = 'S21';
else
	twoPortDirStartString = 'P2_';
    aypy_Label = 'a2 p2';
    bxpy_Label = 'b1 p2';
    sxy_Label = 'S12';
end

freq_GHz = 68:2:116; % VNA frequencies to be analyzed in GHz
nFreq = length(freq_GHz);

nMeas = 40; % Number of sweep points per frequency point (CW sweep mode)
nMeasStartPosition = 2; % if source is unstable at the begining the start position can be adjusted to a higher index, where the source is stable
nMeasUsed = nMeas - nMeasStartPosition + 1;

% DUT = {'Thru_01','6dB_236275','10dB_0006','20dB_236276','20dB_236276_10dB_0006', 'Thru_02'};
DUT = {'Thru_01','6dB_236275','20dB_236276','10dB_0006_20dB_236276'}; % sort with increasing attenuation
nDUT = length(DUT);

sourceAttenuator_dB = [0, 0.5, 1, 2, 4, 6, 8, 10, 12, 14, 16, 18, 20, 22, 24, 26, 28, 30]; % low to high
nPower = length(sourceAttenuator_dB);

fitSourceAttenuatorMin = 8; % Lowest source attenuator position in dB used for the fit of the quadratic model (e.g. 8 dB), Start with a low value and increase the value if the fit of the data is not good
fitSourceAttMinIndex = find(sourceAttenuator_dB == fitSourceAttenuatorMin);

plotLimdB = 0.01; % +/- min/max for the linearity plots

covarianceWeigthing = true; % Covariance weighting of the residuals (e.g. true)

%% Load all data
aypy_all = zeros(nPower, nDUT, nFreq, nMeasUsed); % sqrt(mW)
sxy_all = zeros(nPower, nDUT, nFreq, nMeasUsed);

for i1 = 1:nPower
    for i2 = 1:nDUT
        for i3 = 1:nFreq
			file = [pwd '\' directory '\' sprintf('%dHz\\%s%gdB\\%s.vdatb', freq_GHz(i3).*1e9, twoPortDirStartString, sourceAttenuator_dB(i1), DUT{i2})];
            d = LoadVNADataAsStruct(file);
            % 'Sy,y', 'Sy,x', 'ay_py'
            aypy_all(i1, i2, i3, :) = d.Data(nMeasStartPosition:end, 3);
            sxy_all(i1, i2, i3, :) = d.Data(nMeasStartPosition:end, 2);
        end
    end
end

bxpy_all = aypy_all .* sxy_all;

%% Mean data from samples
aypy = LinProp(zeros(nPower, nDUT, nFreq)); % V/sqrt(kOhm)
bxpy = LinProp(zeros(nPower, nDUT, nFreq)); % V/sqrt(kOhm)
sxy = LinProp(zeros(nPower, nDUT, nFreq));

for i1 = 1:nPower
    for i2 = 1:nDUT
        for i3 = 1:nFreq
            info = sprintf('Source Attenuator %d, DUT %s', sourceAttenuator_dB(i1), DUT{i2});
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

%% Compute slopes
pySlope = LinProp(zeros(nDUT, nFreq));
for i2 = 1:nDUT
    for i3 = 1:nFreq
        pySlope(i2, i3) = lscov(aypy(fitSourceAttMinIndex:end, i2, i3), bxpy(fitSourceAttMinIndex:end, i2, i3), 1e9.*get_covariance(bxpy(fitSourceAttMinIndex:end, i2, i3)));
    end
end

%% Fit all at once
% Model for each data set
% b + x_b * b^2 = x_dut_i * (a + x_a * a^2)

xOpt_py = LinProp(zeros(nDUT + 2, nFreq));
for i3 = 1:nFreq
    xStart = [pySlope(:, i3); 0; 0];
    a_i3 = aypy(fitSourceAttMinIndex:end, :, i3);
    b_i3 = bxpy(fitSourceAttMinIndex:end, :, i3);
    pFlat = [a_i3(:); b_i3(:)];
    xOpt_py(:, i3) = optimizer(@objective_function, xStart, pFlat, covarianceWeigthing);
end
x_dut_p1 = xOpt_py(1:nDUT, :);
x_a_p1 = xOpt_py(end - 1,:);
x_b_p1 = xOpt_py(end, :);

%% Fit
aypyfitlin = LinProp(zeros(nPower, nDUT, nFreq));
bxpyfitlin = LinProp(zeros(nPower, nDUT, nFreq));
residuals_p1fit = LinProp(zeros(nPower, nDUT, nFreq));
for i2 = 1:nDUT
    for i3 = 1:nFreq
        aypyfitlin(:, i2, i3) = (bxpy(:, i2, i3) + x_b_p1(i3).*bxpy(:, i2, i3).^2)./x_dut_p1(i2, i3);
        bxpyfitlin(:, i2, i3) = x_dut_p1(i2, i3).*(aypy(:, i2, i3) + x_a_p1(i3).*aypy(:, i2, i3).^2);
        residuals_p1fit(:, i2, i3) = x_dut_p1(i2, i3).*(aypy(:, i2, i3) + x_a_p1(i3).*aypy(:, i2, i3).^2) ...
                                   - (bxpy(:, i2, i3) + x_b_p1(i3).*bxpy(:, i2, i3).^2);
    end
end
aypyfitlindBm = 20.*log10(aypyfitlin);
bxpyfitlindBm = 20.*log10(bxpyfitlin);

%% Model
aypy_model_lin_dBm = zeros(nPower, nFreq);
aypy_model_2nd_dBm = LinProp(zeros(nPower, nFreq));
bxpy_model_lin_dBm = zeros(nPower, nFreq);
bxpy_model_2nd_dBm = LinProp(zeros(nPower, nFreq));
for i3 = 1:nFreq
    aypy_model_lin_dBm(:, i3) = linspace(min(min(double(aypydBm(:, :, i3)))), max(max(double(aypydBm(:, :, i3)))), nPower);
    aypy_model_lin_i3 = 10.^(aypy_model_lin_dBm(:, i3)./20); % swrt(mW)
    aypy_model_2nd_i3 = aypy_model_lin_i3 + x_a_p1(i3).*aypy_model_lin_i3.^2;
    aypy_model_2nd_dBm(:, i3) = 20.*log10(aypy_model_2nd_i3);
    bxpy_model_lin_dBm(:, i3) = linspace(min(min(double(bxpydBm(:, :, i3)))), max(max(double(bxpydBm(:, :, i3)))), nPower);
    bxpy_model_lin_i3 = 10.^(bxpy_model_lin_dBm(:, i3)./20); % swrt(mW)
    bxpy_model_2nd_i3 = bxpy_model_lin_i3 + x_b_p1(i3).*bxpy_model_lin_i3.^2;
    bxpy_model_2nd_dBm(:, i3) = 20.*log10(bxpy_model_2nd_i3);
end
aypy_model_diff_dB = aypy_model_lin_dBm - aypy_model_2nd_dBm;
bxpy_model_diff_dB = bxpy_model_lin_dBm - bxpy_model_2nd_dBm;

%% Plots per Frequency
aypydBm_min = min(double(aypydBm(:)));
aypydBm_max = max(double(aypydBm(:)));
bxpydBm_min = min(double(bxpydBm(:)));
bxpydBm_max = max(double(bxpydBm(:)));

for i3 = 1:nFreq
    h1 = figure();
    subplot(3,2,1);
    plotv(aypy(:,:,i3), bxpy(:,:,i3));
    xlabel([aypy_Label ' / sqrt(mW)']);
    ylabel([bxpy_Label ' / sqrt(mW)']);
    grid on;

    subplot(3,2,2);
    plotv(aypydBm(:,:,i3), bxpydBm(:,:,i3));
    xlabel([aypy_Label ' / dBm']);
    ylabel([bxpy_Label ' / dBm']);
    xlim([aypydBm_min aypydBm_max]);
    ylim([bxpydBm_min bxpydBm_max]);
    grid on;

    subplot(3,2,3);
    plotu([aypydBm(:,:,i3) aypy_model_lin_dBm(:, i3)], [(aypydBm(:,:,i3) - aypyfitlindBm(:,:,i3)) aypy_model_diff_dB(:, i3)]);
    xlabel([aypy_Label ' / dBm']);
    ylabel(['Non Linearity ' aypy_Label ' / dB']);
    xlim([aypydBm_min aypydBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,4);
    plotu([bxpydBm(:,:,i3) bxpy_model_lin_dBm(:, i3)], [(bxpydBm(:,:,i3) - bxpyfitlindBm(:,:,i3)) bxpy_model_diff_dB(:, i3)]);
    xlabel([bxpy_Label ' / dBm']);
    ylabel(['Non Linearity ' bxpy_Label ' / dB']);
    xlim([bxpydBm_min bxpydBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,5);
    plotu(aypydBm(:,:,i3), residuals_p1fit(:,:,i3));
    xlabel([aypy_Label ' / dBm']);
    ylabel('Residuals');
    xlim([aypydBm_min aypydBm_max]);
    %ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,6);
    plotu(bxpydBm(:,:,i3), residuals_p1fit(:,:,i3));
    xlabel([bxpy_Label ' / dBm']);
    ylabel('Residuals');
    xlim([bxpydBm_min bxpydBm_max]);
    %ylim([-plotLimdB plotLimdB]);
    grid on;

    a = axes;
    t1 = title({strrep([vnaDevice ', ' twoPortDirStartString ' @ ' sprintf('%g', freq_GHz(i3)) ' GHz'], '_', ' '), method, ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h1, [vnaDevice '_' method '_' outputfileNameSupplement]);
end
ps2pdf([vnaDevice '_' method '_' outputfileNameSupplement])

%% Objective function
% b + x_b * b^2 = x_dut_i * (a + x_a * a^2)

function f = objective_function(x, p)

nx = length(x);
np = length(p);
np05 = np / 2;

nDUT = nx - 2;
nPower = np05 / nDUT;

x_dut = x(1:nDUT);
x_a = x(end - 1);
x_b = x(end);

a = reshape(p(1:np05), [nPower, nDUT]);
b = reshape(p(np05 + 1:end), [nPower, nDUT]);

f = LinProp(zeros(nPower, nDUT));

for i2 = 1:nDUT
    f(:, i2) = x_dut(i2).*(a(:, i2) + x_a.*a(:, i2).^2) - (b(:, i2) + x_b.*b(:, i2).^2);
end

f = f(:);
end

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
