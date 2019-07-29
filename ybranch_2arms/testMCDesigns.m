% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');

h = appopen('mode');
appevalscript(h,char('load("ybranch.lms");'));
appevalscript(h,char('addpath("../");'));

config = ybranch_init;
baseconts = config.baseconts;

if (triangleDesign == true) % let's try a different starting point!
    tStart = 216;
    tEnd = 850;
    tLen = tEnd - tStart+1;
    baseconts{1}(tStart:tEnd,1) = linspace(baseconts{1}(tStart,1),baseconts{1}(tEnd,1),tLen);
    baseconts{1}(tStart:tEnd,2) = linspace(baseconts{1}(tStart,2),baseconts{1}(tEnd,2),tLen);
    baseconts{1}(tStart-10:tStart+10,2) = smooth(baseconts{1}(tStart-10:tStart+10,2),7);
    baseconts{1}(tEnd-10:tEnd+10,2) = smooth(baseconts{1}(tEnd-10:tEnd+10,2),7);
    baseconts{4}(tStart:tEnd,1) = linspace(baseconts{4}(tStart,1),baseconts{4}(tEnd,1),tLen);
    baseconts{4}(tStart:tEnd,2) = linspace(baseconts{4}(tStart,2),baseconts{4}(tEnd,2),tLen);
    baseconts{4}(tStart-10:tStart+10,2) = smooth(baseconts{4}(tStart-10:tStart+10,2),7);
    baseconts{4}(tEnd-10:tEnd+10,2) = smooth(baseconts{4}(tEnd-10:tEnd+10,2),7);
    config.baseconts = baseconts;
end

[~, norms_fine] = delts_norms(baseconts{:});

%%%% design configuration %%%%

L_FINE_OUT = size(baseconts{1},1);
L_FINE_STOP = 1000;

L_INTERP = 10; % for every point in the design contour, rough contour has L_INTERP points
L_DESN_RATE = 3; % how many points in the design contour does it take before you're allowed to go from 1 to 0 in slope?
STEP_SIZE = 10e-9;

L_DESN_FULL = floor(L_FINE_OUT / L_INTERP);
L_DESN = floor(L_FINE_STOP / L_INTERP);
splinmat = zeros(L_DESN,L_FINE_OUT);
tmpmat = eye(L_DESN_FULL);
for i = 1:L_DESN
    splinmat(i,:) = interp1(linspace(0,1,L_DESN_FULL),tmpmat(i,:),linspace(0,1,L_FINE_OUT),'spline');
end

splinmats = cell(4,1);
splinmats{1} = splinmat * config.signs(1);
splinmats{2} = zeros(L_DESN,size(baseconts{2},1));
splinmats{3} = zeros(L_DESN,size(baseconts{3},1));
splinmats{4} = splinmat * config.signs(4);

N_DESNS = size(desns,1);

N_ROUGH = 100;
ROUGH_AMPL = 10e-9;
ROUGH_CL = 50e-9;
rough_flucts = cell(config.C_CONTS,1);
flucts = cell(config.C_CONTS,1);
    
for k = 1:config.C_CONTS
    randomflucts = randomRoughness(baseconts{k},ROUGH_CL,N_ROUGH);
    L_FINE = size(baseconts{k},1);
    rough_flucts{k} = [zeros(1,L_FINE) ; ROUGH_AMPL * randomflucts .* config.tapers_fine{k}];
    flucts{k} = zeros(N_DESNS * (1+N_ROUGH),L_FINE);
    for m = 1:N_DESNS
        flucts{k}((m-1)*(1+N_ROUGH)+1:m*(1+N_ROUGH),:) = rough_flucts{k} + desns(m,:) * splinmats{k};
    end
end

hist_ss = run_fwds(h, config, flucts);

save("testmcfinal.mat");