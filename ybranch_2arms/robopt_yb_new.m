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
taper_desn = [sin(linspace(0,pi/2,4)).^2 ones(1,L_DESN-8) cos(linspace(0,pi/2,4)).^2];
lim_desn = STEP_SIZE * taper_desn.';
dlim_desn = (lim_desn(1:end-1) + lim_desn(2:end)) / (2*L_DESN_RATE);

%%%% ALGORITHM CONFIGURATION %%%%

ROUGH_AMPL = 10e-9;
ROUGH_CL = 50e-9;

N_ITER = 30;
N_ROUGH = 3;

N_CANDIDATES = 20;
N_TOPSET = 5; % for picking worst-case scenarios

% In this script, we only modify bot_out and top_out, and we do so
% symmetrically for optimization but separately for roughening.

% Also in this script, we only worry about frequency #5: 1550 nm.

% Also in this script, we only allow design evolution along the normal
% vectors to the original base contour.

C_OUTS = config.C_OUTS;
C_CONTS = config.C_CONTS;

hist_roughs = cell(C_CONTS,1);
hist_desns = zeros(N_ITER,L_DESN);
hist_Ss = zeros(N_ITER,N_ROUGH+1,C_OUTS);
hist_adjSs = cell(C_CONTS,1);
hist_erfs = zeros(N_ITER,N_ROUGH+1);
hist_grads = zeros(N_ITER,N_ROUGH+1,L_DESN);

rconts = cell(C_CONTS,1);
corrmats = cell(C_CONTS,1);
for k = 1:C_CONTS
    hist_roughs{k} = zeros(N_ITER,N_ROUGH+1,size(baseconts{k},1));
    hist_adjSs{k} = zeros(N_ITER,N_ROUGH+1,C_OUTS,size(baseconts{k},1));
end

[~, norms_fine] = delts_norms(baseconts{:});

cur_desn = zeros(L_DESN,1);
cur_fine = cell(C_CONTS,1);

for i = 1:N_ITER
    
    hist_desns(i,:) = cur_desn;
    
    cur_fine{1} = splinmat.' * cur_desn * config.signs(1);
    cur_fine{2} = zeros(size(baseconts{2},1),1);
    cur_fine{3} = zeros(size(baseconts{3},1),1);
    cur_fine{4} = splinmat.' * cur_desn * config.signs(4);
    % GEOM SENSITIVE
    
    if i == 1
        % We start the scenario pool with random sampling
        for k = 1:config.C_CONTS
            [randomflucts, corrmats{k}] = randomRoughness(baseconts{k},ROUGH_CL,N_ROUGH);
            hist_roughs{k}(i,2:end,:) = ROUGH_AMPL * randomflucts .* config.tapers_fine{k};
        end
    end
    
    for j = 1:N_ROUGH+1
        for k = 1:config.C_CONTS
            rconts{k} = baseconts{k} + (squeeze(hist_roughs{k}(i,j,:)) + cur_fine{k}) .* norms_fine{k};
        end

        adjstr = run_adjoint(h,config,rconts,norms_fine);
        hist_Ss(i,j,:) = adjstr.ss;
        Ts = abs(adjstr.ss).^2;
        
        % T = ss*, so dT = 2 re{s' s*}
        dTs = cell(config.C_CONTS,1);
        for k = 1:config.C_CONTS
            hist_adjSs{k}(i,j,:,:) = adjstr.adjss{k};
            dTs{k} = squeeze(2 * real(hist_adjSs{k}(i,j,:,:) .* conj(hist_Ss(i,j,:))));
        end
        
        hist_erfs(i,j) = errfunc(Ts);

        % errfuncgrad tells you how much moving in the normal direction to
        % a contour will affect your error function
        % since we're making a symmetric device, the motions are opposed
        hist_grads(i,j,:) = splinmat * (errfuncgrad(dTs{1}, Ts) - errfuncgrad(dTs{4}, Ts)).'; % GEOM SENSITIVE
        
        if j == 1 % we gotta update the scenario set!!
            candidates = cell(config.C_CONTS,1);
            forecastTs = zeros(C_OUTS,N_CANDIDATES) + Ts;
            ratings = zeros(N_CANDIDATES,1);
            for k = 1:C_CONTS
                candidates{k} = ROUGH_AMPL * randomRoughness(baseconts{k},ROUGH_CL,N_CANDIDATES) .* config.tapers_fine{k};
                for m = 1:N_CANDIDATES
                    forecastTs(:,m) = forecastTs(:,m) + sum(dTs{k} .* candidates{k}(m,:),2);
                end
            end
            for m = 1:N_CANDIDATES
                ratings(m) = errfunc(forecastTs(:,m));
            end
            [~, inds] = maxk(ratings, N_TOPSET);
            ratings = -Inf * ones(N_TOPSET,1);
            for m = 1:N_TOPSET
                for l = 1:N_ROUGH+1
                    score = 0;
                    for k = 1:C_CONTS
                        score = score + candidates{k}(inds(m),:) * squeeze(hist_roughs{k}(i,l,:));
                    end
                    ratings(m) = max(ratings(m), score);
                end
            end
            [~, bestind] = min(ratings);
            for k = 1:C_CONTS
                if i > 1
                    hist_roughs{k}(i,:,:) = hist_roughs{k}(i-1,:,:);
                end
                hist_roughs{k}(i,mod(i,N_ROUGH)+2,:) = candidates{k}(inds(bestind),:);
            end
        end
        
        disp("finished rough round "+num2str(i)+"."+num2str(j));
    end
    
    save("databackup.mat");
    
    SCALE_HACK = 1e5;
    nextstep = robustLinprog_deltas(hist_erfs(i,:).',shiftdim(hist_grads(i,:,:),1)/SCALE_HACK,lim_desn*SCALE_HACK,dlim_desn*SCALE_HACK) / SCALE_HACK;
    cur_desn = cur_desn + nextstep;
end

save("datafinal.mat");

function grad_fine = pessimgrad(dTs, Ts) % function is for a SINGLE CONTOUR and, evidently, a single wavelength
    % we avoid the symmetry problem by ignoring loss and not squaring
    % anything
    imbal_weight = 0.3;
    grad_fine = imbal_weight * (dTs(2,:) - dTs(1,:)); % + (dTs{k}(2,:) + dTs{k}(1,:)) * (Tlo + Thi - 1));
end

function out = errfunc(Ts)
    imbal_weight = 0.3;
    Tlo = Ts(1);
    Thi = Ts(2);
    out = (Tlo + Thi - 1)^2 + imbal_weight * (Tlo - Thi)^2;
end

function grad_fine = errfuncgrad(dTs, Ts) % function is for a SINGLE CONTOUR.
    imbal_weight = 0.3;
    Tlo = Ts(1);
    Thi = Ts(2);
    
   	% error function = (Tlo + Thi - 1)^2 + imbal_weight * (Tlo - Thi)^2
    % so gradient = 2(Tlo + Thi - 1)dT +/- imbal_weight * 2(Tlo - Thi)dT
    
    dTlo_coeff = 2 * (Tlo + Thi - 1 + imbal_weight * (Tlo - Thi));
    dThi_coeff = 2 * (Tlo + Thi - 1 - imbal_weight * (Tlo - Thi));
    
    grad_fine = dTlo_coeff * dTs(1,:) + dThi_coeff * dTs(2,:);
end
