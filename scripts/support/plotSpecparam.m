% fooof_plot() - Plot a FOOOF model.
%
% Usage:
%   >> fooof_plot(fooof_results)
%
% Inputs:
%   fooof_results   = struct of fooof results
%                       Note: must contain FOOOF model, not just results

function plotSpecparam(specparamResults, vis)
    %% Data Checking
    if ~isfield(specparamResults, 'freqs')
       error('Specparam results struct does not contain model output.')
    end

    %% Set Up
    plt_freqs = specparamResults.freqs;

    % Plot settings
    lw = 2.5;

    %% Create the plots
    if vis; figure;
    else; figure('visible', 'off') ;
    end
    hold on

    % Plot the original data
    data = plot(plt_freqs, specparamResults.power_spectrum, 'black');

    % Plot the full model fit
    model = plot(plt_freqs, specparamResults.fooofed_spectrum, 'red');

    % Plot the aperiodic fit
    ap_fit = plot(plt_freqs, specparamResults.ap_fit, 'b--');

    %% Plot Settings
    % Apply general plot settings
    for plt = [data, model, ap_fit]
        set(plt, 'LineWidth', lw);
    end

    % Set alpha value for model
    model.Color(4) = 0.5;

    legend('Original Spectrum', 'Full Model Fit', 'Aperiodic Fit')
    hold off
end