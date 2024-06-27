% fooof_version() - Get FOOOF version information.
%
% Usage:
%   >> fooof_version()
%
% Notes
%   This function returns versions of both the installed FOOOF in Python, and of this Matlab wrapper.
%     The version number of the Matlab FOOOF wrapper is defined in this function.

function fooof_version()

    % Get the version of the Python implementation of FOOOf being used
    fooof_py_version = string(py.pkg_resources.get_distribution('fooof').version);

    % Set the version number of the matlab wrapper
    %   Note: this is where version is defined for this wrapper
    fooof_wrapper_version = "1.0.0";

    % Display versions
    disp('Current version of FOOOF:     ' + fooof_py_version)
    disp('Current version of Wrapper:   ' + fooof_wrapper_version)

end
