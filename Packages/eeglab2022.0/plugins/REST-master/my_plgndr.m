function varargout = my_plgndr ( varargin )

% PLGNDR associated Legendre function
%
% y = plgndr(n,k,x) computes the values of the associated Legendre functions
% of degree N and order K
%
% Implemented as MEX file

% The original implementation was based on "Numerical Recipes in C",
% version 2.0 but has been replaced with an equvalent function from GNU
% Scientific Library,

% Based on FieldTrip 20160222 functions:
% * plgndr by Robert Oostenveld and Thomas Hartmann


% Gets the name of the function and its file name.
funname = mfilename;
funfile = mfilename ( 'fullpath' );
fundir  = fileparts ( funfile );

% The source code must have the same name finished in '.c'.
srcfile = strcat ( funfile, '.c' );

% Tries to compile the code.
try
    % Goes to the source folder.
    curdir  = pwd;
    cd ( fundir );
    
    % Compiles the code.
    mex ( srcfile );
    
    % Returns to the original folder.
    cd ( curdir );
    
catch err
    % Returns to the original folder.
    cd ( curdir );
    
    % Rises an error.
    error ( 'Code for ''%s'' could not be compiled: %s.', funname, err.message );
end

% Executes the newly compilated function.
funhand = str2func ( funname );
[ varargout{ 1: nargout } ] = funhand ( varargin {:} );
