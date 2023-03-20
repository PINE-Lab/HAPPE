function ret=jisubplot(rows,columns,pane,orientation,gap,varargin)
% JISUBPLOT  Set up a figure for multi-pane plotting. Used with NEXTPLOT.
%
%   axis=jisubplot(rows,columns,[pane],[orientation],[gap],['param1','value1',...])
%
%   OVERVIEW: Enhanced replacement for SUBPLOT
%
%       Problem:
%           Multi-axis plotting (using SUBPLOT) is great and essential, but has
%           shortcomings. The main one is that keeping track of the axis indexes
%           can be tedious if you want to do anything more elaborate than plotting
%           in row-order. In addition, the size of all axes are the same, and
%           spacing is not adjustable.
%
%       Solution: JISUBPLOT + NEXTPLOT
%           What you might really want is to set up a grid of subplots once and
%           then just tell the figure when you want to move to the next plot,
%           letting it take care of the details for you.
%
%           That's what the pair of functions JISUBPLOT and NEXTPLOT can do.
%
%           You can move by row, or by column, or arbitrarily. You can create
%           subplots of different sizes. Finally, wouldn't it be nice if
%           when you changed a figure's orientation (e.g. orient tall),
%           it actually changed shape to reflect that?
%
%       Examples:
%
%           figure
%           jisubplot(4,2,1)        % can be used just like SUBPLOT
%           plot(X)
%           nextplot                % advance to next axis
%           plot(something_else)
%           nextplot
%           plot(some_other_thing)
%
%           nextplot('newRow')      % start a new row of axes
%           plot(something_new)
%           nextplot('byCol')       % move down columns
%
%           %% jisubplot / nextplot is especially useful in loops
%           figure
%           jisubplot(4,4,0)
%           for idx = 1:32,
%               nextplot
%               plot(data(idx,:))
%           end
%
%           %% a more complete usage: specify orientation, plot spacing and fontsize
%           figure
%           jisubplot(5,3,0,'tall',[.3 .3],'fontsize',9)
%           nextplot
%
%       See JISUBPLOTDEMO, NEXTPLOT and CURRENTPLOTIS for more advanced usage.
%
%       Description:
%
%       Sets up current figure for multi-pane plotting in a rows x columns grid.
%       To advance among sub-plots, use NEXTPLOT, or specify a value for pane.
%       If pane is unspecified (or=0), doesn't immediately create an axis for
%       plotting, but initializes current figure so that a subsequent call to
%       NEXTPLOT will create the first axis (pane 1). Note, this will also clear
%       the figure. If pane > rows*columns, or NEXTPLOT advances past the current
%       figure, will automatically overflow to a new figure with the same properites.
%
%       If orientation is specified, the figure will resize to have the correct
%       aspect ratio for that type of orientation (e.g. 'tall' or 'landscape').
%       The spacing between plots can be specified (gap), as well as param/value
%       pairs to pass on to each axis.
%
%
%   USAGE:
%
%   axis=jisubplot(rows,columns,[pane],[orientation],[gap],['param1','value1',...])
%
%   INPUTS
%       rows, columns   number of axes in horizontal and vertical directions
%       pane            current sub-plot (if 0 or unspecified, initialze figure)
%                           pane can be specified in two ways:
%                               1) scalar, pick that pane (as in subplot)
%                               2) vector, [toprow leftcol width height] 
%                               row,col of upper left and width & height of axis
%       orientation     'portrait', 'tall', 'landscape' or '' (default is portrait)
%       gap             [x y] space between axes as proportion of axis size
%                           (if unspecified, or [], default to [0.2 0.2].
%                           Note: as is, this is optimized for a high density of
%                           plots, without axis labels. Increase the gap to allow
%                           room for axis labels and titles, etc, use a smaller
%                           axis font size, or use CURRENTPLOTIS to selectively
%                           add labels based on axis position on page.)
%       prop,value      any additional prop/value pairs will be applied to each 
%                           new axis
%
%   OUTPUT
%       axis            handle to current axis (if non-zero value for pane was given)
%
%
%   NOTES
%       
%       if pane is > rows*columns: will automatically overflow to a new figure
%           created with the same properties.
%
%
%   See also NEXTPLOT, CURRENTPLOTIS, JISUBPLOTDEMO.
%
%   John Iversen iversen@nsi.edu
%

%  JRI 4/99    John Iversen, iversen@nsi.edu
%  JRI 4/21/00 added array to store handles to axes. if hold is 'on'
%      activate existing axis, if not, delete and create new one
%  12/16/00 allow separate specification of x & y gap
%  10/11/01 Port to NSI
%  12/11/01 More flexible specification of pane
%       Figure retains its state in appdata, but is historically compatible with
%           older matlab versions lacking appdata (uses userdata in this case).

% Free for all uses, but please retain the following:
%   Original Author: John Iversen
%   john_iversen@post.harvard.edu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set defaults
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

defaultgap = [.25 .25];
usingdefaultgap = 0;
if nargin == 0,
    rows=1;
    columns=1;
    figure
end
if nargin < 3, pane = 0; end
if nargin < 4 || isempty(orientation),
    orientation='portrait';
end
if nargin < 5 || isempty(gap),
    gap=defaultgap;
    usingdefaultgap=1;
end
if nargin >  5, %optional param/value pairs specified
    params = varargin;
else
    params = [];
end

if gap==-1,
    error('use [] for unspecified gap from now on (not -1)')
end

%note, the following refers to system-wide global, G, specific to my personal
%   environment. If this global contains fields G.plot.layout.xplotregion,
%   G.plot.layout.yplotregion and G.plot.layout.superimpose, will use these values;
%   if it's not defined, will set defaults. This can be removed if it conflicts
%   with your personal style, or another use of global G

global G

if isempty(G) | ~isfield(G.plot,'layout') %take care of case when global not defined
    xplotregion = [.05 .95];    %the region on the page in which subplots will fit
    yplotregion = [.05 .95];    %   normalized [0,1] units
else %grab vales from global
    yplotregion = G.plot.layout.yplotregion;
    xplotregion = G.plot.layout.xplotregion;
end
xwidth=xplotregion(2)-xplotregion(1);
ywidth=yplotregion(2)-yplotregion(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%get current state, or initialize if it's a new figure
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

currentfig = gcf; %(creates a figure if there is none)
try
    UD = getappdata(currentfig,'JRI_jisubplotData');
    haveAppdata = 1;
catch
    UD = get(currentfig,'userdata');
    haveAppdata = 0;
end
if isempty(UD), %no userdata means it's a fresh figure, so initialize it
    isVisible = get(currentfig,'visible'); %preserve its (in)visibility
    axisHandles=initializeFigure(currentfig, rows, columns, orientation,callername);
    set(currentfig,'visible',isVisible)
else %if have userdata, unpack it, check for consistency
    axisHandles = UD.axisHandles;
    if rows~=UD.rows | columns~=UD.columns, %if we've changed the number of panes
        axisHandles=zeros(rows,columns); %the old handles aren't useful
    end
    orientation = UD.orientation;
    if nargin>5, %if we've specified extra params, override earlier ones
        params=varargin;
    else
        params = UD.params;
    end
    if usingdefaultgap, %gap wasn't passed in explicitly,
        gap=UD.gap; %so use saved value
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%calculate extent of axis, given the pane number
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   pane is 4-element vector specifying row,col of upper left corner
%       and width and height in panes

%expand scalar value of pane
pl = length(pane);
switch pl,
    case 1,
        if pane ~= 0,
            col=mod(pane-1,columns)+1;
            row=floor((pane-1)/columns)+1;
            pane = [row col 1 1];
        else
            pane = [0 0 1 1];
        end
    case 4,
        %do nothing
    otherwise
        error('Incorrect specification of pane: must be scalar or 4-element')
end

%Calculate axis size
ppp=rows*columns;
xsize=1/columns;
ysize=1/rows;
dx=xsize*gap(1)/2;
dy=ysize*gap(2)/2;

%handle cases where we move to a new figure (if pane moves past current figure)
%   for multi-pane axes, this is computed based on upper left pane number
if pane(1) > 0,
    ulpane = (pane(1)-1) * columns + pane(2);
    fig=max(floor((ulpane-1)/ppp),0)+currentfig;
    %quick fix for advancing off the page in column mode 
    if pane(2)>columns & pane(1)==1,
        fig = currentfig+1;
        ulpane=ppp+1;
    end

    if (fig ~= currentfig), %if we're on a new figure, create and initialize
        ulpane=ulpane-(ppp*(fig-currentfig));
        pane(1)=mod(ulpane-1,columns)+1;
        pane(2)=floor((ulpane-1)/columns)+1;
        %grab suptitle & timestampPrefix from previous figure
        timestampTag = getappdata(currentfig,'JRI_timestampTag');
        axisHandles=initializeFigure(fig,rows,columns,orientation,timestampTag);
        %if previous figure was hidden, hide this one, too
        set(fig,'visible',get(currentfig,'visible'));
        currentfig = fig;
  %      jisuptitle(G.plot.layout.jisuptitlestr)
    end
end

%bounds checking
ulrow = pane(1);
ulcol = pane(2);
w = pane(3);
h = pane(4);
if any([w h]>1),
    if (w > ( 1 + columns - ulcol)) || ...
            (h > (1 + rows - ulrow)),
        error('width or height exceed available space on figure')
    end
end

if (pane(1)>0),	%create an axis
    currentaxis=axisHandles(pane(1),pane(2));

    %for existing axes, use hold state to determine whether to add to axis
    %  or create a new one
    if currentaxis ~= 0,
        try
            superimpose = strcmp(get(currentaxis,'nextplot'),'add') & ...
                strcmp(get(get(currentaxis,'parent'),'nextplot'),'add');
        catch %current axis is invalid (user may have deleted it), so set to zero
            currentaxis = 0;
            superimpose = 0;
        end
    end

    if currentaxis==0 | superimpose==0,	%no current axis or no superimpose:
        %we will make a new axis

        if currentaxis ~=0, %we're not suprimposing, so delete old axis
            try
                delete(currentaxis);
            catch
                currentaxis = 0;
            end
        end

        %create new axis
        row = rows - (pane(1) + (pane(4) - 1));
        col = pane(2) - 1;
        rect=[(col*xsize+dx)*xwidth+xplotregion(1) ... %x
            (row*ysize+dy)*ywidth+yplotregion(1) ... %y
            (xsize*pane(3)-dx*2)*xwidth... %width
            (ysize*pane(4)-dy*2)*ywidth]; %height

        aa=axes('position',rect);
        %axis equal
        axisHandles(pane(1), pane(2))=aa;
        if ~isempty(params)
            set(aa,params{:});
        end

    else %we have a current axis and want to superimpose
        set(currentfig,'CurrentAxes',currentaxis);
        hold on
        aa=currentaxis;
    end
else %pane == 0, don't create an axis
    aa=[];
end

%save figure userdata
UD.rows = rows;
UD.columns = columns;
UD.pane = pane;
UD.gap = gap;
UD.orientation = orientation;
UD.axisHandles = axisHandles;
UD.params = params;

if haveAppdata,
    setappdata(currentfig,'JRI_jisubplotData',UD);
else
    set(currentfig,'userdata',UD);
end

if nargout,
    ret=aa;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function axisHandles=initializeFigure(currentfig, rows, columns,orientation,timestampTag)
%setup new figure's userdata and size
%return matrix of axis handles (initialized to zero)

%optionally initialize text to add to timestamp (only works if have *appdata functions)

axisHandles=zeros(rows,columns);
%set proportions of figure
figure(currentfig); clf reset
try
    rmappdata(currentfig,'JRI_jisubplotData');
end

if nargin < 5
    timestampTag = '';
    try
        rmappdata(currentfig,'JRI_timestampTag')
    end
end

pos=get(currentfig,'position');
pos(1)=16 + 16;     %+ 16*currentfig; % to stagger

%adjust figure size according to screen size, and # screens
saveUnits = get(0,'units');
set(0,'units','pixels');
ss = get(0,'screensize');

% get multiple monitor sizes (updated 3/26/13--different sized monitors)
mp = get(0,'MonitorPositions');
nMonitors = size(mp,1);
if nMonitors > 1,
    mp(2,2) = mp(1,4)-mp(2,4); %correct y origin of second monitor
    figureOrigin = [mp(2,1) + 64 mp(1,4) - 32-80]; %upper left corner of secondary monitor, with room for toolbars
else
    figureOrigin = [mp(3)*1/2 mp(4)-32-80];
end
%want tall oriented figure to be 85% of screen height
max_size = mp(1,4) * 0.85;
%on a big screen, it's too big, so cap it.
max_size = min(max_size,820);


%move it up a bit
%pos(2) = pos(2)-ss(4)/3;

switch orientation
    case 'portrait'

    case 'tall'
        pos(3)= max_size * 0.75;
        pos(4) = max_size;

    case 'landscape'
        pos(3) = max_size;
        pos(4) = max_size * 0.75;
        if numScreens > 1,
         % pos(1) = max(8,pos(1) - 16); %shift slightly left of tall, but not off screen
        end

    case ''
        orientation='portrait';	%default

    otherwise
        error('Invalid orientation')
end

%adjust X to be on left screen (assuming main is right screen)

%pos(1) = pos(1) - ss(3)/2 * (numScreens-1);
pos(1) = figureOrigin(1);

pos(2) = figureOrigin(2) - pos(4); %adjust to lower left of figure

%keep figure on screen: Y > 0, Y+height < screen height
figtop = pos(2)+pos(4)+80; %80 is for toolbars
if figtop > ss(4) + 10,
    pos(2) = 0.95*ss(4) - figtop;
end
if pos(2) < 0,
    pos(2) = ss(4) * 0.05; %keep bottom slightly above screen bottom
end
set(currentfig,'position',pos)
orient(currentfig, orientation);

set(0,'units',saveUnits)

%set figure background to white
set(currentfig,'color',[1 1 1])

if ~isempty(timestampTag),
    setappdata(currentfig,'JRI_timestampTag',timestampTag)
    timestamp
end
