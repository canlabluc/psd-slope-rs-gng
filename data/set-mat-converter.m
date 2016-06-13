files = dir(fullfile('extClfiltCAR-set/*.set'));

for id = 1:numel(files)
    EEG = pop_loadset(files(id).name, 'extClfiltCAR-set/');
    fields_to_remove = {'filepath', 'subject', 'group', 'condition', 'session', 'comments', 'trials', 'xmin', 'xmax', 'times', 'icaact', 'icawinv', 'icasphere', 'icaweights', 'icachansind', 'urchanlocs', 'event', 'urevent', 'eventdescription', 'epoch', 'epochdescription', 'reject', 'stats', 'specdata', 'specicaact', 'splinefile', 'icasplinefile', 'dipfit', 'history', 'saved', 'etc', 'datfile'};
    EEG = rmfield(EEG, fields_to_remove);
    name = EEG.setname;
    srate = EEG.srate;
    data  = EEG.data;
    save(strcat('extClfiltCAR-mat/', files(id).name(1:9), '.mat'), 'name', 'srate', 'data'); 
end