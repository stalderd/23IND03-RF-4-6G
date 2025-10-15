% Analyze Linearity Data
% Michael Wollensack METAS - 23.12.2021 - 12.08.2025

clear all
close all

%% Variables
vnaDevice = 'Simulation';
name = 'StepAtt';

dataSets = {'LinTestDataSet1', 'LinTestDataSet2', 'LinTestDataSet3', 'LinTestDataSet4', 'LinTestDataSet5', 'LinTestDataSet6', 'LinTestDataSet7'};
nDataSets = length(dataSets);

freq = 1e9;
nFreq = length(freq);

power = -30:1:0;
nPower = length(power);
stepAtt = 0:10:70;
nStepAtt = length(stepAtt);
nMeas = 100;
fitPowerMax = -15;
plotLimdB = 0.1;

covarianceWeigthing = false;

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
a1p1m = LinProp(zeros(nPower, nStepAtt, nDataSets)); % sqrt(mW)
b2p1m = LinProp(zeros(nPower, nStepAtt, nDataSets)); % sqrt(mW)

for i1 = 1:nPower
    for i2 = 1:nStepAtt
        for i3 = 1:nDataSets
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
p1Slope = LinProp(zeros(nStepAtt, nDataSets));
for i2 = 1:nStepAtt
    for i3 = 1:nDataSets
        p1Slope(i2, i3) = lscov(a1p1m(1:fitPowerMaxIndex, i2, i3), b2p1m(1:fitPowerMaxIndex, i2, i3), 1e9.*get_covariance(b2p1m(1:fitPowerMaxIndex, i2, i3)));
    end
end

%% Fit all at once
% Model for each data set
% b + x_b * b^2 = x_dut_i * (a + x_a * a^2)

xOpt_p1 = LinProp(zeros(nStepAtt + 2, nDataSets));
for i3 = 1:nDataSets
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
a1p1fitlin = LinProp(zeros(nPower, nStepAtt, nDataSets));
b2p1fitlin = LinProp(zeros(nPower, nStepAtt, nDataSets));
residuals_p1fit = LinProp(zeros(nPower, nStepAtt, nDataSets));
for i2 = 1:nStepAtt
    for i3 = 1:nDataSets
        a1p1fitlin(:, i2, i3) = (b2p1m(:, i2, i3) + x_b_p1(i3).*b2p1m(:, i2, i3).^2)./x_dut_p1(i2, i3);
        b2p1fitlin(:, i2, i3) = x_dut_p1(i2, i3).*(a1p1m(:, i2, i3) + x_a_p1(i3).*a1p1m(:, i2, i3).^2);
        residuals_p1fit(:, i2, i3) = x_dut_p1(i2, i3).*(a1p1m(:, i2, i3) + x_a_p1(i3).*a1p1m(:, i2, i3).^2) ...
                                   - (b2p1m(:, i2, i3) + x_b_p1(i3).*b2p1m(:, i2, i3).^2);
    end
end
a1p1fitlindBm = 20.*log10(a1p1fitlin);
b2p1fitlindBm = 20.*log10(b2p1fitlin);

%% Model
a1p1_model_lin_dBm = zeros(nPower, nDataSets);
a1p1_model_2nd_dBm = LinProp(zeros(nPower, nDataSets));
b2p1_model_lin_dBm = zeros(nPower, nDataSets);
b2p1_model_2nd_dBm = LinProp(zeros(nPower, nDataSets));
for i3 = 1:nDataSets
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

for i3 = 1:nDataSets
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
    t1 = title({strrep([vnaDevice ' @ ' dataSets{i3}], '_', ' '), ''});
    a.Visible = 'off';
    t1.Visible = 'on';

    figure2ps(h1, [vnaDevice '_opt']);
end
ps2pdf([vnaDevice '_opt'])

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