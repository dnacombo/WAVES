%% here I will simulate a traveling wave in an occipital region
% project to channels, do this for N trials, K subjs
% then attempt

dogif = 0;
%% setup environment
addpath(fullfile(cd,'fieldtrip'))
addpath(fullfile(cd,'misc'))
global ft_default
ft_default.trackcallinfo = 'no';
ft_default.showcallinfo = 'no';
ft_defaults

% load EEGLAB & associated plugins
add_eeglab
% eeglab;
rm_frompath('eeglab.*Fieldtrip-lite')
rm_frompath('SEREEGA')

addpath(genpath(fullfile(cd,'SEREEGA')));
%% load sensors geometry
YO = '/home/maximilien.chaumon/owncloud/Lab/Projects/SOLID_MEEG/data/DS/151014/YO_tsss_ica_blink_cardio_corrected.fif';
cfg = [];
cfg.dataset = YO;
cfg.channel = 'meg';
data = ft_preprocessing(cfg);
data.grad = ft_convert_units(data.grad,'mm');
%% create head model
cfg = [];
cfg.method = 'singlesphere';
headmodel = ft_prepare_headmodel(cfg,data.grad);

%% create source model
cfg = [];
% cfg.grid.xgrid      = 'auto';
% cfg.grid.ygrid      = 'auto';
% cfg.grid.zgrid      = 'auto';
cfg.grid.resolution = 20;
cfg.headmodel = headmodel;
lf = ft_prepare_leadfield(cfg,data);

%% create SEREEGA compatible leadfield from fieldtrip leadfield
leadfield = [];
leadfield.leadfield = permute(cat(3,lf.leadfield{lf.inside}),[1 3 2]);
leadfield.orientation = zeros(size(leadfield.leadfield,2),3);
leadfield.pos = lf.pos(lf.inside,:);
for ichan = 1:numel(data.grad.label)
    leadfield.chanlocs(ichan).labels = data.grad.label{ichan};
    leadfield.chanlocs(ichan).X = data.grad.chanpos(ichan,1);
    leadfield.chanlocs(ichan).Y = data.grad.chanpos(ichan,2);
    leadfield.chanlocs(ichan).Z = data.grad.chanpos(ichan,3);
end


%% template leadfield (EEG)
% leadfield = lf_generate_fromfieldtrip( 'montage', 'BioSemi64', 'resolution',10);
% neighbours = ft_prepare_neighbours(struct('method','template','template','biosemi64_neighb.mat'));


sourcesIdx = 564:568;

if dogif
    % h = figure(299);clf
    plot_source_location(sourcesIdx, leadfield,'mode','3d');
    hold on
    ft_plot_sens(data.grad)
    addpath(cdhome('general/plot/export_fig'))
    nstep = 100;
    az = linspace(90,450,nstep);
    for i = 1:numel(az)
        view(az(i),30);
        drawnow
        export_fig(sprintf('tmp%02d',i),'-nocrop');
    end
    !gifski -o sources.gif tmp*.png
    delete tmp*.png
end
%% define signal of interest

modLatency_ori = 500;

ersp.type       ='ersp';
ersp.frequency  = 10;
ersp.amplitude  = 5e-3;
ersp.modulation = 'burst';
ersp.modLatency = modLatency_ori;        % centre of the burst, in ms
ersp.modWidth   = 200;        % width (half duration) of the burst, in ms
ersp.modTaper   = 0.5;        % taper of the burst
ersp = utl_check_class(ersp);

% plot_signal_fromclass(ersp, epochs);
% axis('tight')
allphases = linspace(0,1/2,numel(sourcesIdx));

%% define epochs parameters

epochs.n        = 100;
epochs.length   = 1000;% ms
epochs.srate    = data.fsample; % Hz
epochs.prestim  = 250; % ms

%% define data class with yo data
dat = [];
dat.data = data.trial{1};
dat.data = reshape(dat.data(:,1:epochs.n*epochs.length/1000*epochs.srate),size(dat.data,1),epochs.length/1000*epochs.srate,epochs.n);
dat.index = {'e', ':',':'};
dat.amplitude = 1;
dat.amplitudeType = 'relative';
dat = utl_check_class(dat,'type','data');
%% now one component per source of the traveling wave
% figure(848);clf
noise = [];
noise.type = 'noise';
noise.color = 'brown';
noise.amplitude = 10;
c = [];
for isource = 1:numel(sourcesIdx)
    c(end+1).source = sourcesIdx(isource);      % obtained from the lead field, as above
    
    ersp.modLatency = modLatency_ori + 10*isource;
    c(end).signal = {ersp noise};
    
    c(end).orientation = [0 0 1];
    %      plot_signal_fromclass(ersp, epochs,'newfig',0);
    
end
c = utl_check_component(c, leadfield);

scalpdata1 = generate_scalpdata(c, leadfield, epochs);% + dat.data;


%%
datasim= utl_create_fieldtripdataset(scalpdata1,epochs,leadfield);

cfg = [];
datasim_tl = ft_timelockanalysis(cfg,datasim);

cfg = [];
cfg.layout = ft_prepare_layout(struct('layout','neuromag306mag.lay'));
figure;
ft_multiplotER(cfg,datasim_tl)

figure;
cfg = [];
cfg.layout = ft_prepare_layout(struct('layout','neuromag306mag.lay'));
cfg.xlim = [.22 .23];
ft_topoplotER(cfg,datasim_tl);

