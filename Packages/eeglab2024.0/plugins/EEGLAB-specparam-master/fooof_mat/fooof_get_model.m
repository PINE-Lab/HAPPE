% fooof_get_model() - Return the model fit values from a FOOOF object
%
% Usage:
%   >> model_fit = fooof_get_model(fm)
%
% Inputs:
%   fm              = FOOOF object
%
% Outputs:
%   model_fit       = model results, in  a struct, including:
%       model_fit.freqs
%       model_fit.power_spectrum
%       model_fit.fooofed_spectrum
%       model_fit.ap_fit
%
% Notes
%   This function is mostly an internal function.
%     It can be called directly by the user if you are interacting with FOOOF objects directly.

function model_fit = fooof_get_model(fm)

    model_fit = struct();

    model_fit.freqs = double(py.array.array('d',fm.freqs));
    model_fit.power_spectrum = double(py.array.array('d', fm.power_spectrum));
    model_fit.fooofed_spectrum = double(py.array.array('d', fm.fooofed_spectrum_));
    model_fit.ap_fit = double(py.array.array('d', py.getattr(fm, '_ap_fit')));