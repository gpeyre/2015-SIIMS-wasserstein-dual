%%
% test for mesh vizualization

test = 'MNE_01';
repin = [test '/visualization/'];
rep = [repin 'results/'];
if not(exist(rep))
    mkdir(rep);
end
addpath('toolbox_visu/');


%%
% Scan directory. First load all datas
A = dir([repin 'LOO*.mat']);

%load([repin 'prototype_graph_0d5']);
load([repin 'mesh']);
load([repin 'edges']);
XY = vertices(:,1:2)';
E = edges'+1;

estim = @(x,i)x(i);
prctile = @(v,p)estim(sort(v(:)), round(length(v(:))*p) );
clamp = @(v,a,b)max(min(v,b),a);
clamp_prc = @(v,eta)clamp( v, prctile(v,eta),prctile(v,1-eta) );
eta = .04;
clampa = @(v)clamp_prc( v, eta );


%%
% Render a few histograms.

load([repin 'data_condensed']);
reph = [repin 'results-histo/'];
if not(exist(reph))
    mkdir(reph);
end
cl = [-1, +1];
nbr = 10;
options.node_width = 100;
if 1
for k=1:2
    I = find(labels==cl(k));
    for i=1:nbr
        a = data(I(i),:);
        a = a/max(a);
        % a = log(a);
        a = clampa(a);
        h = rescale( a );
        clf;
        plot_mne(XY,h);
%        plot_graph_wire(XY,E, [], h, options);
        saveas(gcf, [reph 'histo-cl' num2str(k) '-' num2str(i) '.eps'], 'epsc');
    end
end
end

%%
% Read results.


% original edge length
edgeLengths0 = sqrt( sum( ( XY(:,E(1,:)) - XY(:,E(2,:)) ).^2 ) );

Metric = [];
Dummy = [];
for imat = 1:length(A)
    filename = A(imat).name;
    name = filename(1:end-4);
    load([repin filename]);
    %
    Dummy(:,end+1) = dummyLengths(:);
    Metric(:,end+1) = edgeLengths(:)./edgeLengths0(:);
end
%Dummy = rescale(Dummy);
%Metric = rescale(Metric);

options.edge_width = 10;
options.node_width = 50;

for imat = 13:length(A)  
    filename = A(imat).name;
    name = filename(1:end-4);  
    
    v = rescale( clamp_prc(Dummy(:,imat), .2) );
    w = rescale( clamp_prc(Metric(:,imat), .2) );
        
    opt.CM = gray(256)/2+1/2;
    opt.edge_width = 6;
    clf; hold on;
    plot_mne(XY,v, opt);
    plot_graph_wire(XY,E, w, [], opt);
    saveas(gcf, [rep name '.eps'], 'epsc');
        
    % only edge
    if 0
    clf;
    plot_graph_wire(XY,E,w, [], options);
    saveas(gcf, [rep name '-edges.eps'], 'epsc');
    end
    % only vertices
    if 0
    clf;
    plot_graph_wire(XY,E,[], v, options);
    saveas(gcf, [rep name '-nodes.eps'], 'epsc');
    end
    % both
    if 0
    clf; hold on;
    plot_mne(XY,w);
    plot_graph_wire(XY,E, w, v, options);
    saveas(gcf, [rep name '-edges-nodes.eps'], 'epsc');
    end
end

%%
% Plot average.

v = rescale( clamp_prc(mean(Dummy,2), .2) );
w = rescale( clamp_prc(mean(Metric,2), .2) );

opt.CM = gray(256)/2+1/2;
opt.edge_width = 6;
clf; hold on;
plot_mne(XY,v, opt);
plot_graph_wire(XY,E, w, [], opt);
saveas(gcf, [rep 'MEANS-edges-nodes.eps'], 'epsc');