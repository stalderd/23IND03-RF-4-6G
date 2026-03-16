% Analyze Linearity Data
% Michael Wollensack METAS - 23.12.2021 - 05.05.2025

clear all
close all
addpath('C:\Users\Public\Documents\Metas.Vna.Tools\Matlab');
LoadVNATools();

%% Variables
vnaDevice = 'Agilent_PNA_N5227A_#1';
name = 'StepAtt60dB(f-f)_SN123456_01';
outputfileNameSupplement = '_opt_P2_01';
freq = [50e6 1e9 10e9 18e9 30e9 40e9 50e9]; % VNA frequencies to be measured
nFreq = length(freq);

power = -30:1:-5; % VNA source power levels to be measured (e.g. -30 dBm to -5 dBm in 1 dB steps)
nPower = length(power);
stepAtt = 0:10:60; % Step attenuator values (only used for naming, e.g. 0 dB to 60 dB it 10 dB steps)
nStepAtt = length(stepAtt);
nMeas = 100; % Number of sweep points per frequency point (CW sweep mode)
fitPowerMax = -5; % Highest power level in dBm used for the fit of the quadratic model (e.g. -5 dBm), Start with a high value and decrease the value if the fit of the data is not good
plotLimdB = 0.01; % Scale limit for the plots in dB (e.g. 0.01 dB)

covarianceWeigthing = true; % Covariance weighting of the residuals (e.g. true)

%% Load all datasqrt(mW)
a1p1 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
b2p1 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
a2p2 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
b1p2 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
s21 = zeros(nPower, nStepAtt, nFreq, nMeas);
s12 = zeros(nPower, nStepAtt, nFreq, nMeas);

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nFreq
            d = LoadVNADataAsStruct([pwd '\Measurements_01\' name '\' sprintf('%dHz\\%ddBm\\%s_%ddB.vdatb', freq(i3), power(i1), name, stepAtt(i2))]);
            % 'S1,1', 'S1,2', 'S2,1', 'S2,2', 'a1_p1', 'a2_p2', 'a2/b2_p1', 'a1/b1_p2'
            a1p1(i1, i2, i3, :) = d.Data(:, 5);
            b2p1(i1, i2, i3, :) = d.Data(:, 3).*d.Data(:, 5);
            a2p2(i1, i2, i3, :) = d.Data(:, 6);
            b1p2(i1, i2, i3, :) = d.Data(:, 2).*d.Data(:, 6);
            s21(i1, i2, i3, :) = d.Data(:, 3);
            s12(i1, i2, i3, :) = d.Data(:, 2);
        end
    end
end

%% Mean data from samples
a2p2m = LinProp(zeros(nPower, nStepAtt, nFreq)); % sqrt(mW)
b1p2m = LinProp(zeros(nPower, nStepAtt, nFreq)); % sqrt(mW)

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nFreq
            info = sprintf('Power %d, Step Att %d dB', power(i1), stepAtt(i2));
            tp2 = LinProp([reshape(real(a2p2(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(a2p2(i1, i2, i3, :)), [nMeas, 1]) ...
                           reshape(real(b1p2(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(b1p2(i1, i2, i3, :)), [nMeas, 1]) ...
                           reshape(real(s12(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(s12(i1, i2, i3, :)), [nMeas, 1])], 'samples', ['S12, ' info]);
            a2p2m(i1, i2, i3) = get_value(abs(tp2(1) + 1i.*tp2(2)));
            b1p2m(i1, i2, i3) = abs(tp2(5) + 1i.*tp2(6)).*a2p2m(i1, i2, i3);
        end
    end
end

%% Mean data in dBm
a2p2dBm = 20.*log10(a2p2m);
b1p2dBm = 20.*log10(b1p2m);

%% Compute slopes
fitPowerMaxIndex = find(power == fitPowerMax);
p2Slope = LinProp(zeros(nStepAtt, nFreq));
for i2 = 1:nStepAtt
    for i3 = 1:nFreq
        p2Slope(i2, i3) = lscov(a2p2m(1:fitPowerMaxIndex, i2, i3), b1p2m(1:fitPowerMaxIndex, i2, i3), 1e9.*get_covariance(b1p2m(1:fitPowerMaxIndex, i2, i3)));
    end
end

%% Fit all at once
% Model for each data set
% b + x_b * b^2 = x_dut_i * (a + x_a * a^2)

xOpt_p2 = LinProp(zeros(nStepAtt + 2, nFreq));
for i3 = 1:nFreq
    xStart = [p2Slope(:, i3); 0; 0];
    a_i3 = a2p2m(1:fitPowerMaxIndex, :, i3);
    b_i3 = b1p2m(1:fitPowerMaxIndex, :, i3);
    pFlat = [a_i3(:); b_i3(:)];
    xOpt_p2(:, i3) = optimizer(@objective_function, xStart, pFlat, covarianceWeigthing);
end
x_dut_p2 = xOpt_p2(1:nStepAtt, :);
x_a_p2 = xOpt_p2(end - 1,:);
x_b_p2 = xOpt_p2(end, :);


%% Fit
a2p2fitlin = LinProp(zeros(nPower, nStepAtt, nFreq));
b1p2fitlin = LinProp(zeros(nPower, nStepAtt, nFreq));
residuals_p2fit = LinProp(zeros(nPower, nStepAtt, nFreq));
for i2 = 1:nStepAtt
    for i3 = 1:nFreq
        a2p2fitlin(:, i2, i3) = (b1p2m(:, i2, i3) + x_b_p2(i3).*b1p2m(:, i2, i3).^2)./x_dut_p2(i2, i3);
        b1p2fitlin(:, i2, i3) = x_dut_p2(i2, i3).*(a2p2m(:, i2, i3) + x_a_p2(i3).*a2p2m(:, i2, i3).^2);
        residuals_p2fit(:, i2, i3) = x_dut_p2(i2, i3).*(a2p2m(:, i2, i3) + x_a_p2(i3).*a2p2m(:, i2, i3).^2) ...
                                   - (b1p2m(:, i2, i3) + x_b_p2(i3).*b1p2m(:, i2, i3).^2);
    end
end
a2p2fitlindBm = 20.*log10(a2p2fitlin);
b1p2fitlindBm = 20.*log10(b1p2fitlin);

%% Model
a2p2_model_lin_dBm = zeros(nPower, nFreq);
a2p2_model_2nd_dBm = LinProp(zeros(nPower, nFreq));
b1p2_model_lin_dBm = zeros(nPower, nFreq);
b1p2_model_2nd_dBm = LinProp(zeros(nPower, nFreq));
for i3 = 1:nFreq
    a2p2_model_lin_dBm(:, i3) = linspace(min(min(double(a2p2dBm(:, :, i3)))), max(max(double(a2p2dBm(:, :, i3)))), nPower);
    a2p2_model_lin_i3 = 10.^(a2p2_model_lin_dBm(:, i3)./20); % swrt(mW)
    a2p2_model_2nd_i3 = a2p2_model_lin_i3 + x_a_p2(i3).*a2p2_model_lin_i3.^2;
    a2p2_model_2nd_dBm(:, i3) = 20.*log10(a2p2_model_2nd_i3);
    b1p2_model_lin_dBm(:, i3) = linspace(min(min(double(b1p2dBm(:, :, i3)))), max(max(double(b1p2dBm(:, :, i3)))), nPower);
    b1p2_model_lin_i3 = 10.^(b1p2_model_lin_dBm(:, i3)./20); % swrt(mW)
    b1p2_model_2nd_i3 = b1p2_model_lin_i3 + x_b_p2(i3).*b1p2_model_lin_i3.^2;
    b1p2_model_2nd_dBm(:, i3) = 20.*log10(b1p2_model_2nd_i3);
end
a2p2_model_diff_dB = a2p2_model_lin_dBm - a2p2_model_2nd_dBm;
b1p2_model_diff_dB = b1p2_model_lin_dBm - b1p2_model_2nd_dBm;

%% Plots P2
a2p2dBm_min = min(double(a2p2dBm(:)));
a2p2dBm_max = max(double(a2p2dBm(:)));
b1p2dBm_min = min(double(b1p2dBm(:)));
b1p2dBm_max = max(double(b1p2dBm(:)));

for i3 = 1:nFreq
    h1 = figure();
    subplot(3,2,1);
    plotv(a2p2m(:,:,i3), b1p2m(:,:,i3));
    xlabel('a2 p2 / sqrt(mW)');
    ylabel('b1 p2 / sqrt(mW)');
    grid on;

    subplot(3,2,2);
    plotv(a2p2dBm(:,:,i3), b1p2dBm(:,:,i3));
    xlabel('a2 p2 / dBm');
    ylabel('b1 p2 / dBm');
    xlim([a2p2dBm_min a2p2dBm_max]);
    ylim([b1p2dBm_min b1p2dBm_max]);
    grid on;

    subplot(3,2,3);
    plotu([a2p2dBm(:,:,i3) a2p2_model_lin_dBm(:, i3)], [(a2p2dBm(:,:,i3) - a2p2fitlindBm(:,:,i3)) a2p2_model_diff_dB(:, i3)]);
    xlabel('a2 p2 / dBm');
    ylabel('Non Linearity a2 p2 / dB');
    xlim([a2p2dBm_min a2p2dBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,4);
    plotu([b1p2dBm(:,:,i3) b1p2_model_lin_dBm(:, i3)], [(b1p2dBm(:,:,i3) - b1p2fitlindBm(:,:,i3)) b1p2_model_diff_dB(:, i3)]);
    xlabel('b1 p2 / dBm');
    ylabel('Non Linearity b1 p2 / dB');
    xlim([b1p2dBm_min b1p2dBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,5);
    plotu(a2p2dBm(:,:,i3), residuals_p2fit(:,:,i3));
    xlabel('a2 p2 / dBm');
    ylabel('Residuals');
    xlim([a2p2dBm_min a2p2dBm_max]);
    %ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,6);
    plotu(b1p2dBm(:,:,i3), residuals_p2fit(:,:,i3));
    xlabel('b1 p2 / dBm');
    ylabel('Residuals');
    xlim([b1p2dBm_min b1p2dBm_max]);
    %ylim([-plotLimdB plotLimdB]);
    grid on;

    a = axes;
    t1 = title({strrep([vnaDevice '_P2 @ ' sprintf('%.3f', freq(i3)/1e9) ' GHz'], '_', ' '), ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h1, [vnaDevice outputfileNameSupplement]);
end
ps2pdf([vnaDevice outputfileNameSupplement])

%% Objective function
% b + x_b * b^2 = x_dut_i * (a + x_a * a^2)

function f = objective_function(x, p)

nx = length(x);
np = length(p);
np05 = np / 2;

nStepAtt = nx - 2;
nPower = np05 / nStepAtt;

x_dut = x(1:nStepAtt);
x_a = x(end - 1);
x_b = x(end);

a = reshape(p(1:np05), [nPower, nStepAtt]);
b = reshape(p(np05 + 1:end), [nPower, nStepAtt]);

f = LinProp(zeros(nPower, nStepAtt));

for i2 = 1:nStepAtt
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
y_std = get_stdunc(y2); % Bug fixed (was x2)
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

t = get(handle,'PaperPosition');
set(handle,'PaperUnits','Centimeters'); % R2023a
set(handle,'PaperPosition',[0.63 0.63 19.73 28.43]);
set(handle,'PaperSize',[21.00 29.70]);
print(handle, '-dpsc', '-r300', '-append', [filename '_temp']);
end

function ps2pdf(filename)

system(['ps2pdf -dEPSCrop -dPDFX ' filename '_temp.ps ' filename '.pdf'],'-echo');
system(['del ' filename '_temp.ps'],'-echo');

end
