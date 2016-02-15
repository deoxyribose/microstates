% eegplugin_microstates()
function eegplugin_microstates( fig, try_strings, catch_strings);
if nargin < 3
    error('eegplugin_microstates requires 3 arguments');
end;
% create menus
plotmenu = findobj(fig, 'tag', 'plot');
toolsmenu = findobj(fig, 'tag', 'tools');
submenu = uimenu( toolsmenu, 'label', 'Microstates');
subsubmenu = uimenu( submenu, 'label', 'Plot');

% menu callback
getmicrostates = [try_strings.no_check '[EEG LASTCOM]=pop_getmicrostates(EEG);' catch_strings.new_and_hist];
plotmicrostates2 = [try_strings.no_check '[EEG LASTCOM]=plotmicrostates2(EEG);' catch_strings.new_and_hist];
%plotspectrogram = [try_strings.no_check 'plotspectrograms(EEG, 1);' catch_strings.new_and_hist];
plotmicrostates = [try_strings.no_check 'plotmicrostates(EEG, 1, 0);' catch_strings.new_and_hist];
plotsegmentation = [try_strings.no_check 'plotsegmentation(EEG);' catch_strings.new_and_hist];
microstate_statistics = [try_strings.no_check 'microstate_statistics(EEG,0);' catch_strings.new_and_hist];
plotgtm = [try_strings.no_check 'global_topographic_measures(EEG, 1);' catch_strings.new_and_hist];
pointcloud = [try_strings.no_check 'show_clusters(EEG);' catch_strings.new_and_hist];


% create menu
uimenu( submenu, 'label', 'Estimate microstates', 'callback', getmicrostates);
uimenu( submenu, 'label', 'Plot microstates', 'callback', plotmicrostates2);
uimenu( submenu, 'label', 'Microstate statistics', 'callback', microstate_statistics);
uimenu( subsubmenu, 'label', 'Microstate maps', 'callback', plotmicrostates);
uimenu( subsubmenu, 'label', 'Microstate segmentation', 'callback', plotsegmentation);
uimenu( subsubmenu, 'label', 'Global topographic measures', 'callback', plotgtm);
uimenu( subsubmenu, 'label', 'Point cloud', 'callback', pointcloud);
%uimenu( subsubmenu, 'label', 'Channel spectrograms', 'callback', plotspectrogram);
end