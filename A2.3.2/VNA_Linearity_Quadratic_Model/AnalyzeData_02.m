% Analyze Linearity Data
% Michael Wollensack METAS - 23.12.2021 - 05.05.2025

clear all
close all
addpath('C:\Users\Public\Documents\Metas.Vna.Tools\Matlab');
LoadVNATools();

%% Variables
vnaDevice = 'Agilent_PNA_N5227A_#1';
name = 'StepAtt60dB(f-f)_SN123456_01';
freq = [50e6 1e9 10e9 18e9 30e9 40e9 50e9]; % VNA frequencies to be measured
nFreq = length(freq);

power = -30:1:0; % VNA source power levels to be measured (e.g. -30 dBm to 0 dBm in 1 dB steps)
nPower = length(power);
stepAtt = 0:10:60; % Step attenuator values (only used for naming, e.g. 0 dB to 60 dB it 10 dB steps)
nStepAtt = length(stepAtt);
nMeas = 100; % Number of sweep points per frequency point (CW sweep mode)
fitPowerMax = -15; % Highest power level in dBm used for the fit of the quadratic model (e.g. -15 dBm)
plotLimdB = 0.01; % Scale limit for the plots in dB (e.g. 0.01 dB)

covarianceWeigthing = true; % Covariance weighting of the residuals (e.g. true)

%% Load all datasqrt(mW)
a1p1 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
b2p1 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
a2p2 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
b1p2 = zeros(nPower, nStepAtt, nFreq, nMeas); % sqrt(mW)
s11 = zeros(nPower, nStepAtt, nFreq, nMeas);
s21 = zeros(nPower, nStepAtt, nFreq, nMeas);
s12 = zeros(nPower, nStepAtt, nFreq, nMeas);
s22 = zeros(nPower, nStepAtt, nFreq, nMeas);

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nFreq
            d = LoadVNADataAsStruct([pwd '\Measurements_01\' name '\' sprintf('%dHz\\%ddBm\\%s_%ddB.vdatb', freq(i3), power(i1), name, stepAtt(i2))]);
            % 'S1,1', 'S1,2', 'S2,1', 'S2,2', 'a1_p1', 'a2_p2', 'a2/b2_p1', 'a1/b1_p2'
            a1p1(i1, i2, i3, :) = d.Data(:, 5);
            b2p1(i1, i2, i3, :) = d.Data(:, 3).*d.Data(:, 5);
            a2p2(i1, i2, i3, :) = d.Data(:, 6);
            b1p2(i1, i2, i3, :) = d.Data(:, 2).*d.Data(:, 6);
            s11(i1, i2, i3, :) = d.Data(:, 1);
            s21(i1, i2, i3, :) = d.Data(:, 3);
            s12(i1, i2, i3, :) = d.Data(:, 2);
            s22(i1, i2, i3, :) = d.Data(:, 4);
        end
    end
end

%% Mean data from samples
a1p1m = LinProp(zeros(nPower, nStepAtt, nFreq)); % sqrt(mW)
b2p1m = LinProp(zeros(nPower, nStepAtt, nFreq)); % sqrt(mW)

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nFreq
            info = sprintf('Power %d, Step Att %d dB', power(i1), stepAtt(i2));
            tp1 = LinProp([reshape(real(a1p1(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(a1p1(i1, i2, i3, :)), [nMeas, 1]) ...
                           reshape(real(b2p1(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(b2p1(i1, i2, i3, :)), [nMeas, 1]) ...
                           reshape(real(s21(i1, i2, i3, :)), [nMeas, 1]) reshape(imag(s21(i1, i2, i3, :)), [nMeas, 1])], 'samples', ['S21, ' info]);
            a1p1m(i1, i2, i3) = get_value(abs(tp1(1) + 1i.*tp1(2)));
            b2p1m(i1, i2, i3) = abs(tp1(5) + 1i.*tp1(6)).*a1p1m(i1, i2, i3);
        end
    end
end

%% Mean data in dBm
a1p1dBm = 20.*log10(a1p1m);
b2p1dBm = 20.*log10(b2p1m);

%% Compute slopes
fitPowerMaxIndex = find(power == fitPowerMax);
p1Slope = LinProp(zeros(nStepAtt, nFreq));
for i2 = 1:nStepAtt
    for i3 = 1:nFreq
        p1Slope(i2, i3) = lscov(a1p1m(1:fitPowerMaxIndex, i2, i3), b2p1m(1:fitPowerMaxIndex, i2, i3), 1e9.*get_covariance(b2p1m(1:fitPowerMaxIndex, i2, i3)));
    end
end

%% Fit all at once
% Model for each data set
% b + x_b * b^2 = x_dut_i * (a + x_a * a^2)

xOpt_p1 = LinProp(zeros(nStepAtt + 2, nFreq));
for i3 = 1:nFreq
    xStart = [p1Slope(:, i3); 0; 0];
    a_i3 = a1p1m(:, :, i3);
    b_i3 = b2p1m(:, :, i3);
    pFlat = [a_i3(:); b_i3(:)];
    xOpt_p1(:, i3) = optimizer(@objective_function, xStart, pFlat, covarianceWeigthing);
end
x_dut_p1 = xOpt_p1(1:nStepAtt, :);
x_a_p1 = xOpt_p1(end - 1,:);
x_b_p1 = xOpt_p1(end, :);


%% Fit
a1p1fitlin = LinProp(zeros(nPower, nStepAtt, nFreq));
b2p1fitlin = LinProp(zeros(nPower, nStepAtt, nFreq));
residuals_p1fit = LinProp(zeros(nPower, nStepAtt, nFreq));
for i2 = 1:nStepAtt
    for i3 = 1:nFreq
        a1p1fitlin(:, i2, i3) = (b2p1m(:, i2, i3) + x_b_p1(i3).*b2p1m(:, i2, i3).^2)./x_dut_p1(i2, i3);
        b2p1fitlin(:, i2, i3) = x_dut_p1(i2, i3).*(a1p1m(:, i2, i3) + x_a_p1(i3).*a1p1m(:, i2, i3).^2);
        residuals_p1fit(:, i2, i3) = x_dut_p1(i2, i3).*(a1p1m(:, i2, i3) + x_a_p1(i3).*a1p1m(:, i2, i3).^2) ...
                                   - (b2p1m(:, i2, i3) + x_b_p1(i3).*b2p1m(:, i2, i3).^2);
    end
end
a1p1fitlindBm = 20.*log10(a1p1fitlin);
b2p1fitlindBm = 20.*log10(b2p1fitlin);

%% Model
a1p1_model_lin_dBm = zeros(nPower, nFreq);
a1p1_model_2nd_dBm = LinProp(zeros(nPower, nFreq));
b2p1_model_lin_dBm = zeros(nPower, nFreq);
b2p1_model_2nd_dBm = LinProp(zeros(nPower, nFreq));
for i3 = 1:nFreq
    a1p1_model_lin_dBm(:, i3) = linspace(min(min(double(a1p1dBm(:, :, i3)))), max(max(double(a1p1dBm(:, :, i3)))), nPower);
    a1p1_model_lin_i3 = 10.^(a1p1_model_lin_dBm(:, i3)./20); % swrt(mW)
    a1p1_model_2nd_i3 = a1p1_model_lin_i3 + x_a_p1(i3).*a1p1_model_lin_i3.^2;
    a1p1_model_2nd_dBm(:, i3) = 20.*log10(a1p1_model_2nd_i3);
    b2p1_model_lin_dBm(:, i3) = linspace(min(min(double(b2p1dBm(:, :, i3)))), max(max(double(b2p1dBm(:, :, i3)))), nPower);
    b2p1_model_lin_i3 = 10.^(b2p1_model_lin_dBm(:, i3)./20); % swrt(mW)
    b2p1_model_2nd_i3 = b2p1_model_lin_i3 + x_b_p1(i3).*b2p1_model_lin_i3.^2;
    b2p1_model_2nd_dBm(:, i3) = 20.*log10(b2p1_model_2nd_i3);
end
a1p1_model_diff_dB = a1p1_model_lin_dBm - a1p1_model_2nd_dBm;
b2p1_model_diff_dB = b2p1_model_lin_dBm - b2p1_model_2nd_dBm;

%% Plots P1
a1p1dBm_min = min(double(a1p1dBm(:)));
a1p1dBm_max = max(double(a1p1dBm(:)));
b2p1dBm_min = min(double(b2p1dBm(:)));
b2p1dBm_max = max(double(b2p1dBm(:)));

for i3 = 1:nFreq
    h1 = figure();
    subplot(3,2,1);
    plotv(a1p1m(:,:,i3), b2p1m(:,:,i3));
    xlabel('a1 p1 / sqrt(mW)');
    ylabel('b2 p1 / sqrt(mW)');
    grid on;

    subplot(3,2,2);
    plotv(a1p1dBm(:,:,i3), b2p1dBm(:,:,i3));
    xlabel('a1 p1 / dBm');
    ylabel('b2 p1 / dBm');
    xlim([a1p1dBm_min a1p1dBm_max]);
    ylim([b2p1dBm_min b2p1dBm_max]);
    grid on;

    subplot(3,2,3);
    plotu([a1p1dBm(:,:,i3) a1p1_model_lin_dBm(:, i3)], [(a1p1dBm(:,:,i3) - a1p1fitlindBm(:,:,i3)) a1p1_model_diff_dB(:, i3)]);
    xlabel('a1 p1 / dBm');
    ylabel('Non Linearity a1 p1 / dB');
    xlim([a1p1dBm_min a1p1dBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,4);
    plotu([b2p1dBm(:,:,i3) b2p1_model_lin_dBm(:, i3)], [(b2p1dBm(:,:,i3) - b2p1fitlindBm(:,:,i3)) b2p1_model_diff_dB(:, i3)]);
    xlabel('b2 p1 / dBm');
    ylabel('Non Linearity b2 p1 / dB');
    xlim([b2p1dBm_min b2p1dBm_max]);
    ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,5);
    plotu(a1p1dBm(:,:,i3), residuals_p1fit(:,:,i3));
    xlabel('a1 p1 / dBm');
    ylabel('Residuals');
    xlim([a1p1dBm_min a1p1dBm_max]);
    %ylim([-plotLimdB plotLimdB]);
    grid on;

    subplot(3,2,6);
    plotu(b2p1dBm(:,:,i3), residuals_p1fit(:,:,i3));
    xlabel('b2 p1 / dBm');
    ylabel('Residuals');
    xlim([b2p1dBm_min b2p1dBm_max]);
    %ylim([-plotLimdB plotLimdB]);
    grid on;

    a = axes;
    t1 = title({strrep([vnaDevice '_P1 @ ' sprintf('%.3f', freq(i3)/1e9) ' GHz'], '_', ' '), ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h1, [vnaDevice '_opt_P1']);
end
ps2pdf([vnaDevice '_opt_P1'])

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
