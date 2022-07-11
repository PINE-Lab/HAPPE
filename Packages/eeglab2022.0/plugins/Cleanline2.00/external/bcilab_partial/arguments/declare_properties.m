function declare_properties(varargin)
% Declare properties of a function.
% declare_properties('Property',value,'Property',value, ...)
%
% This construct can be used to declare any key-value properties for a function. It must not follow
% after a call to arg_define() but before it. These properties can then be retrieved using
% arg_report('properties',@myfunction).
%
% BCILAB looks for certain properties in filter (signal processing) and dataset editing functions
% (and possibly other functions in the future), to support improved integration of these filters 
% into the GUI and the evaluation system. See flt_pipeline for more details.
%
% In:
%   Options...: list of name-value pairs.
%
%               Ordering hints applicable for filters (respected by flt_pipeline) are the following:
%                'depends': expresses that some other filter(s) *must* have been called before running 
%                           the filter whose properties are being declared (this is usually due to some 
%                           required meta-data that is supplied by that other filter)
%                'cannot_follow': expresses that the filter being declared cannot follow after some 
%                                 other filter(s) (this is usually due to destructive
%                                 editing performed by those other filters)
%                'cannot_precede': expresses that the filter being declared cannot precede some 
%                                  other filter(s)
%                'follows': expresses that the filter being declared *prefers* to follow some 
%                           other filter(s) (this is usually for numeric or efficiency reasons)
%                'precedes': expresses that the filter being declared *prefers* to precede some 
%                           other filter(s) (this is usually for numeric or efficiency reasons)
%
%               Optional properties respected by the evaluation system (for improved usability):
%                'independent_channels': specify that this filter does not transfer information
%                                        across channels (e.g. channel selection, FIR filter)
%                                        (allows the online system to auto-infer which data channels 
%                                        are actually required by a given BCI model)
%                'independent_trials': specify that this filter does transfer information
%                                      across second-length or larger time scales (on the order of 
%                                      the duration of a trial or larger); this determines whether
%                                      the filter will be called repeatedly for each partition of the 
%                                      data in a cross-validation and other offline analyses
%
%               Further optional properties include:
%                'name': specify the GUI/human-readable name of this filter (respected by flt_pipeline)
%
%
% Examples:
%   function myfunction(...)
%
%   % declare a property called 'name', and assign the string 'MyFunction' to it, and another
%   % property called 'price', with some value attached
%   declare_properties('name','MyFunction','price',1999)
%
%   % in a signal processing function, declare the name of the function as it should appear in the 
%   % GUI panel (default is the name of the file without the flt_ / set_ prefix)
%   declare_properties('name','MyFilter')
%
%   % in a signal processing function, declare that the filter depends on two other (previously applied)
%   % filters (here: flt_ica and set_fit_dipoles)
%   declare_properties('depends',{'flt_ica','set_flt_dipoles'});
%   
%   % in a signal processing function, declare a hint that the filter typically follows some other 
%   % filter in the pipeline (if that other filter was enabled; here: flt_resample)
%   declare_properties('follows','flt_resample');
%
%   % in a signal processing function, declare a hint that the filter typically precedes some other 
%   % filters in the pipeline (if that other filter was enabled; here: flt_iir and flt_fir)
%   declare_properties('precedes',{'flt_iir','flt_fir'});
%
%   % in a signal processing function, declare the constraint that the filter cannot follow some other 
%   % filter
%   declare_properties('cannot_follow','set_makepos');
%
%   % in a signal processing function, declare the constraint that the filter cannot follow some other 
%   % filter
%   declare_properties('cannot_precede','set_makepos');
%
%   % in a signal processing function, declare that the contents of each output channel produced 
%   % by the filter depend only on the contents of the corresponding input channel (i.e. there is not
%   % cross-mixing of channels (default is false)
%   declare_properties('independent_channels',true);
%
%   % in a signal processing function, declare that the contents of each output trial depend not
%   % only on data of the respective input trial, but perhaps on other trials, as well (i.e., there is 
%   % cross-mixing over time at trial granularity); the effect of this is that the respective filter
%   % will be executed once per cross-validation fold
%   declare_properties('independent_trials',false);
%   
% See also:
%   arg_report
%
%                                Christian Kothe, Swartz Center for Computational Neuroscience, UCSD
%                                2010-11-02

if evalin('caller','exist(''arg_nvps'',''var'')')
    warn_once('BCILAB:arg:declare_properties','The function declare_properties must not be be called after arg_define().'); end

try
    throw; %#ok<LTARG> % faster than error()
catch context
    if any(strcmp('arg_report_properties',{context.stack(3:min(6,end)).name}))
        % properties are being requested via arg_report; return them as a struct
        caller = hlp_getcaller();
        
        % collect properties and assign defaults
        properties = hlp_varargin2struct(varargin, ...
            'depends',{}, ...
            'cannot_follow',{}, ...
            'cannot_precede',{}, ...
            'follows',{}, ...
            'precedes',{}, ...
            'independent_channels',[], ...
            'independent_trials',[], ...
            'name', caller);
        
        if isempty(properties.independent_trials) && (strncmp(caller,'flt_',4) || strncmp(caller,'set_',4))
            properties.independent_trials = false;
            % warn about this omission: in practice, this means that this filter and all that come after it
            % in a pipeline will be recomputed for every partition of the data during (nested) cross-validation
            % (perhaps 10x or 100x as often as necessary)
            disp_once('Note: The function "%s" does not declare the property independent_trials; assuming that it is false. This may be a conservative assumption that can have massive performance cost during offline processing. Consider declaring the property explicitly in the declare_properties() clause.',caller);
        end
        
        if isempty(properties.independent_channels) && (strncmp(caller,'flt_',4) || strncmp(caller,'set_',4))
            properties.independent_channels = false;
            % warn about this omission: in practice, this often means that BCILAB has to assume that a given
            % model requires all channels that were present in the source data set -- even if it selects
            % only a subset somewhere in its filter pipeline (e.g., excluding peripheral measures). If such
            % channels are not available in an online stream, the model will fail to work, perhaps
            % unneccessarily.
            disp_once('Note: The function "%s" does not declare the property independent_channels; assuming that it is false. This may be a conservative assumption that can make it unnecessariy difficult to set up online processing. Consider declaring the property explicitly in the declare_properties() clause.',caller);
        end
        
        % report them
        arg_issuereport(properties);
    end
end
