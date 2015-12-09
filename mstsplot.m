% timtopo()   - plot all channels of a data epoch on the same axis
%               and map its scalp map(s) at selected latencies.
% Usage:
%  >> timtopo(data, chan_locs);
%  >> timtopo(data, chan_locs, 'key', 'val', ...);
% Inputs:
%  data       = (channels,frames) single-epoch data matrix
%  chan_locs  = channel location file or EEchanlocs structure.
%               See >> topoplot example for file format.
%
% Optional ordered inputs:
%  'limits'    = [minms maxms minval maxval] data limits for latency (in ms) and y-values
%                 (assumes uV) {default|0 -> use [0 npts-1 data_min data_max];
%                 else [minms maxms] or [minms maxms 0 0] -> use
%                [minms maxms data_min data_max]
%  'plottimes' = [vector] latencies (in ms) at which to plot scalp maps
%                {default|NaN -> latency of maximum variance}
% 'title'      = [string] plot title {default|0 -> none}
% 'plotchans'  = vector of data channel(s) to plot. Note that this does not
%                affect scalp topographies {default|0 -> all}
% 'voffsets'   = vector of (plotting-unit) distances vertical lines should extend
%                above the data (in special cases) {default -> all = standard}
%
% Optional keyword, arg pair inputs (must come after the above):
% 'topokey','val' = optional topoplot() scalp map plotting arguments. See >> help topoplot
%
% Author: Scott Makeig, SCCN/INC/UCSD, La Jolla, 1-10-98
%
% See also: envtopo(), topoplot()

% Copyright (C) 1-10-98 Scott Makeig, SCCN/INC/UCSD, scott@sccn.ucsd.edu
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% 5-31-00 added o-time line and possibility of plotting 1 channel -sm & mw
% 11-02-99 added maplimits arg -sm
% 01-22-01 added to help message -sm
% 01-25-02 reformated help & license, added link -ad
% 03-15-02 add all topoplot options -ad

function [OUTEEG, com] = mstsplot(INEEG)

com = '';

OUTEEG = INEEG;

MAX_TOPOS = 24;
data = INEEG.data;
srate = INEEG.srate;
chan_locs = INEEG.chanlocs;
A = INEEG.A;
idx = INEEG.idx;
times = INEEG.times;
if nargin < 1 %should this be 2?
    help timtopo;
    return
end

[chans,frames] = size(data);
icadefs;

plotframes=find(diff([0; idx']));
encoding = [idx(plotframes)', diff([plotframes; size(idx,2)+1])]; % run length encoding

xmin = times(1);
xmax = times(end);
ymax=max(max(data));
ymin=min(min(data));
%
%%%%%%%%%%%%%%% Compute plot times/frames %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
ntopos = size(A,2);
if ntopos > MAX_TOPOS
    fprintf('timtopo(): too many plottimes - only first %d will be shown!\n',MAX_TOPOS);
    plottimes = plottimes(1:MAX_TOPOS);
    ntopos = MAX_TOPOS;
end

%
%%%%%%%%%%%%%%%%  Compute title and axes font sizes %%%%%%%%%%%%%%%
%
%pos = get(gca,'Position');
pos = [0.1300 0.1100 0.7750 0.8150];
axis('off')
cla % clear the current axes
if pos(4)>0.70
    titlefont= 16;
    axfont = 16;
elseif pos(4)>0.40
    titlefont= 14;
    axfont = 14;
elseif pos(4)>0.30
    titlefont= 12;
    axfont = 12;
elseif pos(4)>0.22
    titlefont= 10;
    axfont = 10;
else
    titlefont= 8;
    axfont = 8;
end

%
%%%%%%%%%%%%%%%% Compute topoplot head width and separation %%%%%%%%%%%%%%%
%
head_sep = 0.2;
topowidth = pos(3)/((6*ntopos-1)/5); % width of each topoplot
if topowidth> 0.25*pos(4) % dont make too large (more than 1/4 of axes width)!
    topowidth = 0.25*pos(4);
end

halfn = floor(ntopos/2);
if rem(ntopos,2) == 1  % odd number of topos
    topoleft = pos(3)/2 - (ntopos/2+halfn*head_sep)*topowidth;
else % even number of topos
    topoleft = pos(3)/2 - ((halfn)+(halfn-1)*head_sep)*topowidth;
end
topoleft = topoleft - 0.01; % adjust left a bit for colorbar

if max(plotframes) > frames |  min(plotframes) < 1
    fprintf('Requested map frame %d is outside data range (1-%d)\n',max(plotframes),frames);
    return
end


%
%%%%%%%%%%%%%%%%%%%%%%% Plot the data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%
%%%%%%%%%%%% site the plot at bottom of the figure %%%%%%%%%%%%%%%%%%
%
figure('Name','Microstates');
%figure()
axdata = axes('Units','Normalized','Position',[pos(1) pos(2) pos(3) 0.6*pos(4)],'FontSize',axfont);
set(axdata,'Color',BACKCOLOR);
limits = get(axdata,'Ylim');
set(axdata,'GridLineStyle',':')
set(axdata,'Xgrid','off')
set(axdata,'Ygrid','on')
axes(axdata)
axcolor = get(gcf,'Color');
set(axdata,'Color',BACKCOLOR);
hold on
cmap = colormap(colorcube);
%image([xmin xmax], [ymin ymax], idx, 'AlphaData',1);
for i=1:(length(plotframes)-1)
    pl=plot(times(plotframes(i):plotframes(i+1)),data(:,plotframes(i):plotframes(i+1))','color',cmap(encoding(i,1),:));    % plot the data
end
pl=plot(times(plotframes(i+1):end),data(:,plotframes(i+1):end)','color',cmap(encoding(i+1,1),:));    % plot the data
xl= xlabel('Latency (ms)');
set(xl,'FontSize',axfont);
yl=ylabel('Potential (\muV)');
set(yl,'FontSize',axfont,'FontAngle','normal');
axis([xmin xmax ymin ymax]);

freezeColors
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw vertical lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
width  = xmax-xmin;
height = ymax-ymin;
lwidth = 1;  % increment line thickness
plotchans = 1:chans;
plottimes = plotframes/srate*1000;
% for t=1:length(plotframes) % draw vertical lines through the data at topoplot frames
%  if length(plotchans)>1
%   l1 = plot([plottimes(t) plottimes(t)],...
%        [min(data(plotchans,plotframes(t))) ...
%        max(data(plotchans,plotframes(t)))],'w'); % white underline behind
%   set(l1,'linewidth',2);
%   l1 = plot([plottimes(t) plottimes(t)],...
%        [min(data(plotchans,plotframes(t))) ...
%        max(data(plotchans,plotframes(t)))],'b'); % blue line
%   set(l1,'linewidth',lwidth);
%  end
% end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Draw oblique lines %%%%%%%%%%%%%%%%%%%%%%%%%%
%
axall = axes('Position',pos,...
    'Visible','Off','FontSize',axfont);   % whole-gca invisible axes
axes(axall)
set(axall,'Color',BACKCOLOR);
axis([0 1 0 1])
axes(axall)
axis([0 1 0 1]);
set(gca,'Visible','off'); % make whole-figure axes invisible


% for t=1:ntopos % draw oblique lines through to the topoplots
%     wherethemstsare = plotframes(encoding(:,1)==t)/srate*1000;
%     for segment=1:length(wherethemstsare)
%         
%         axtp = axes('Units','Normalized','Position',...
%             [pos(1)+topoleft+(t-1)*(1+head_sep)*topowidth ...
%             pos(2)+0.66*pos(4) ...
%             topowidth ...
%             topowidth*(1+head_sep)]); % this will be the topoplot axes
%         % topowidth]); % this will be the topoplot axes
%         axis([-1 1 -1 1]);
%         from = changeunits([wherethemstsare(segment) ,ymax],axdata,axall); % data axes
%         to   = changeunits([0,-0.74],axtp,axall);                % topoplot axes
%         delete(axtp);
%         axes(axall);                                             % whole figure axes
%         l1 = plot([from(1) to(1)],[from(2) to(2)],'b');
%         set(l1,'linewidth',lwidth);
%         
%         hold on
%         set(axall,'Visible','off');
%         axis([0 1 0 1]);
%     end
% end
%
%%%%%%%%%%%%%%%%%%%%%%%%% Plot the topoplots %%%%%%%%%%%%%%%%%%%%%%%%%%
%
topoaxes = zeros(1,ntopos);
for t=1:ntopos
    % [pos(3)*topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
    axtp = axes('Units','Normalized','Position',...
        [topoleft+pos(1)+(t-1)*(1+head_sep)*topowidth ...
        pos(2)+0.66*pos(4) ...
        topowidth topowidth*(1+head_sep)]);
    axes(axtp)                             % topoplot axes
    topoaxes(t) = axtp; % save axes handles
    cla
    
    topoplot(A(:,t),chan_locs, 'electrodes', 'off', 'style', 'map','headrad','rim','hcolor',cmap(t,:),'hlinewidth',5);
    
    
    
    % Else make a 3-D headplot
    %
    % headplot(data(:,plotframes(t)),'chan.spline');
    
    %timetext = [num2str(plottimes(t),'%4.0f')];
    % timetext = [num2str(plottimes(t),'%4.0f') ' ms']; % add ' ms'
    text(0.00,0.80,['Microstate ', num2str(t)],'FontSize',axfont-3,'HorizontalAlignment','Center'); % ,'fontweight','bold');
end

%
% Turn on axcopy()
%

% clicking on ERP pop_up topoplot
% % -------------------------------
% disp('Click on ERP waveform to show scalp map at specific latency');
% 
% dat.times = times;
% dat.erp   = data;
% dat.chanlocs = chan_locs;
% dat.options  = {'electrodes', 'off'};
% dat.srate    = srate;
% dat.axes     = axtp;
% dat.line     = l1;
% 
% cb_code = [ 'tmppos = get(gca, ''currentpoint'');' ...
%     'dattmp = get(gcf, ''userdata'');' ...
%     'set(dattmp.line, ''visible'', ''off'');' ...
%     'axes(dattmp.axes); cla;' ...
%     'latpoint = round((tmppos(1)-dattmp.times(1))/1000*dattmp.srate);' ...
%     'topoplot(dattmp.erp(:,latpoint), dattmp.chanlocs, dattmp.options{:});' ...
%     'title(sprintf(''%.0f ms'', tmppos(1)));' ...
%     'clear latpoint dattmp;' ...
%     ];
% axcopy;
% 
% set(gcf, 'userdata', dat);
% set(axdata, 'ButtonDownFcn', cb_code); %windowbuttondownfcn', cb_code);
% set(pl, 'ButtonDownFcn', cb_code);

%axcopy(gcf, cb_code);
