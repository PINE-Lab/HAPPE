function varargout = hlp_scope(assignments, f, varargin)
% Execute a function within a dynamic scope of values assigned to symbols.
% Results... = hlp_scope(Assignments, Function, Arguments...)
%
% This is the only completely reliable way in MATLAB to ensure that symbols that should be assigned
% while a function is running get cleared after the function returns orderly, crashes, segfaults,
% the user slams Ctrl+C, and so on. Symbols can be looked up via hlp_resolve().
%
% In:
%   Assignments : Cell array of name-value pairs or a struct. Values are associated with symbols of
%                 the given names. The names should be valid MATLAB identifiers. These assigments
%                 form a dynamic scope for the execution of the function; scopes can also be
%                 nested, and assignments in inner scopes override those of outer scopes.
%
%   Function  : a function handle to invoke
%
%   Arguments... : arguments to pass to the function
%
% Out:
%   Results... : return value(s) of the function
%
% See also:
%   hlp_resolve
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-05-03

% add a new stack frame with the evaluated assignments & get its unique id
id = make_stackframe(assignments);
% also take care that it gets reclaimed after we're done
reclaimer = onCleanup(@()return_stackframe(id));

% make a function that is tagged by id
func = make_func(id);

% evaluate the function with the id introduced into MATLAB's own stack
[varargout{1:nargout}] = func(f,varargin);



function func = make_func(id)
persistent funccache; % (cached, since the eval() below is a bit slow)
try
    func = funccache.(id);
catch
    func = eval(['@(f,a,frame__' id ')feval(f,a{:})']);
    funccache.(id) = func;
end


function id = make_stackframe(assignments)
% put the assignments into a struct
if iscell(assignments)
    assignments = cell2struct(assignments(2:2:end),assignments(1:2:end),2); end
% get a fresh frame id
global tracking;
try
    id = tracking.stack.frameids.removeLast();
catch
    if ~isfield(tracking,'stack') || ~isfield(tracking.stack,'frameids')        
        % need to create the id repository first
        tracking.stack.frameids = java.util.concurrent.LinkedBlockingDeque();
        for k=50000:-1:1
            tracking.stack.frameids.addLast(sprintf('f%d',k)); end
    else
        if tracking.stack.frameids.size() == 0
            % if this happens then either you have 10.000s of parallel executions of hlp_scope(),
            % or you have a very deep recursion level (the MATLAB default is 500), or your function
            % has crashed 10.000s of times in a way that keeps onCleanup from doing its job, or you have
            % substituted onCleanup by a dummy class or function that doesn't actually work (e.g. on
            % pre-2008a systems).
            error('We ran out of stack frame ids. This should not happen under normal conditions. Please make sure that your onCleanup implementation is not consistently failing to execute.'); 
        end
    end
    id = tracking.stack.frameids.removeLast();
end
% and store the assignments under it
tracking.stack.frames.(id) = assignments;


function return_stackframe(id)
% finally return the frame id again...
global tracking;
tracking.stack.frameids.addLast(id);
