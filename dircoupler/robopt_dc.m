% path(path,'/opt/lumerical/mode/7.12.1661-1/api/matlab');

h = appopen('mode');
appevalscript(h,char('load("dircoupler_narrow.lms");'));
appevalscript(h,char('addpath("../");'));

config = dc_init(200, 0.075e-6, 0.575e-6);
baseconts = config.baseconts;

L_FINE = size(baseconts{1},1);

L_INTERP = 10; % for every point in the design contour, rough contour has L_INTERP points

L_DESN = floor(L_FINE / L_INTERP);
splinmat = zeros(L_DESN,L_FINE);
tmpmat = eye(L_DESN);
for i = 1:L_DESN
    splinmat(i,:) = interp1(linspace(0,1,L_DESN),tmpmat(i,:),linspace(0,1,L_FINE),'spline');
end
taper_desn = [sin(linspace(0,pi/2,3)).^2 ones(1,L_DESN-6) cos(linspace(0,pi/2,3)).^2];
taper_fine = taper_desn * splinmat;

%%%% ALGORITHM CONFIGURATION %%%%

ROUGH_AMPL = 3e-9;
ROUGH_CL = 50e-9;

N_ITER = 30;
N_ROUGH = 3;
STEP_SIZE = 5e-9;
USE_OPPOS = false;

% In this script, we only allow symmetrical designs, but asymmetrical roughenings.

% In this script, we only allow design evolution along the normal
% vectors to the original base contour.

C_OUTS = config.C_OUTS;
C_CONTS = config.C_CONTS;
N_FREQS = config.N_FREQS; % careful...

hist_roughs = cell(C_CONTS,1);
hist_desns = zeros(N_ITER,L_DESN,2); % first outside, then inside
hist_Ss = zeros(N_ITER,N_ROUGH+1,C_OUTS,N_FREQS);
hist_adjSs = cell(C_CONTS,1);
hist_erfs = zeros(N_ITER,N_ROUGH+1);
hist_grads = zeros(N_ITER,N_ROUGH+1,L_DESN,2);

rconts = cell(C_CONTS,1);
corrmats = cell(C_CONTS,1);
for k = 1:C_CONTS
    hist_roughs{k} = zeros(N_ITER,N_ROUGH+1,L_FINE);
    hist_adjSs{k} = zeros(N_ITER,N_ROUGH+1,C_OUTS,L_FINE,N_FREQS);
end

[~, norms_fine] = delts_norms(baseconts{:});

cur_desn = zeros(L_DESN,2);
cur_fine = cell(C_CONTS,1);

for i = 1:N_ITER
    
    hist_desns(i,:,:) = cur_desn;
    % GEOM SENSITIVE
    cur_fine{1} = splinmat.' * cur_desn(:,1) * config.signs(1);
    cur_fine{2} = splinmat.' * cur_desn(:,2) * config.signs(2);
    cur_fine{3} = splinmat.' * cur_desn(:,2) * config.signs(3);
    cur_fine{4} = splinmat.' * cur_desn(:,1) * config.signs(4);
    
    if i == 1
        % We start the scenario pool with random sampling
        for k = 1:config.C_CONTS
            [randomflucts, corrmats{k}] = randomRoughness(baseconts{k},ROUGH_CL,N_ROUGH);
            hist_roughs{k}(i,2:end,:) = ROUGH_AMPL * randomflucts .* taper_fine;
        end
    end
    
    for j = 1:N_ROUGH+1
        for k = 1:config.C_CONTS
            rconts{k} = baseconts{k} + (squeeze(hist_roughs{k}(i,j,:)) + cur_fine{k}) .* norms_fine{k};
        end % squeeze considered sketch, but 1x1xL_FINE should work out ok

        adjstr = run_adjoint(h,config,rconts,norms_fine);
        hist_Ss(i,j,:,:) = adjstr.ss;
        Ts = abs(adjstr.ss).^2;
        
        % T = ss*, so dT = 2 re{s' s*}
        dTs = cell(config.C_CONTS,1);
        for k = 1:config.C_CONTS
            hist_adjSs{k}(i,j,:,:,:) = adjstr.adjss{k};
            dTs{k} = shiftdim(2 * real(permute(hist_adjSs{k}(i,j,:,:,:),[1 2 3 5 4]) .* conj(hist_Ss(i,j,:,:))),2);
        end % 1x1xC_OUTxL_FINExN_FREQS times 1x1xC_OUTxN_FREQS --> C_OUT x N_FREQS x L_FINE
        
        hist_erfs(i,j) = errfunc(Ts);

        % errfuncgrad tells you how much moving in the normal direction to
        % a contour will affect your error function
        % since we're making a symmetric device, the motions are opposed
        hist_grads(i,j,:,1) = splinmat * (errfuncgrad(dTs{1}, Ts) - errfuncgrad(dTs{4}, Ts)).'; % GEOM SENSITIVE
        hist_grads(i,j,:,2) = splinmat * (errfuncgrad(dTs{3}, Ts) - errfuncgrad(dTs{2}, Ts)).'; % GEOM SENSITIVE
        
        if j == 1 % we gotta update the scenario set!!    
            grad_fine = cell(config.C_CONTS,1);
            for k = 1:C_CONTS
                % for pessimizing, we use a different error function to avoid
                % the symmetry problem
                grad_fine{k} = pessimgrad(dTs{k}, Ts);
                if i > 1
                    if USE_OPPOS
                        hist_roughs{k}(i,4:end,:,:) = hist_roughs{k}(i-1,2:end-2,:);
                    else
                        hist_roughs{k}(i,3:end,:,:) = hist_roughs{k}(i-1,2:end-1,:);
                    end
                end
                hist_roughs{k}(i,2,:) = (ROUGH_AMPL * sqrt(L_FINE) * grad_fine{k} * corrmats{k} / sqrt(grad_fine{k} * corrmats{k} * grad_fine{k}.')) .* taper_fine;
                if USE_OPPOS
                    hist_roughs{k}(i,3,:) = -1 * hist_roughs{k}(i,2,:);
                end
            end
        end
        
        disp("finished rough round "+num2str(i)+"."+num2str(j));
    end
    
    save("databackup.mat");
    
    SCALE_HACK = 1e5;
    combograds = [shiftdim(hist_grads(i,:,:,1),1) shiftdim(hist_grads(i,:,:,2),1)];
    nextstep = robustLinprog(hist_erfs(i,:).', combograds/SCALE_HACK,STEP_SIZE*SCALE_HACK) / SCALE_HACK;
    
    cur_desn(:,1) = cur_desn(:,1) + nextstep(1:L_DESN) .* taper_desn.';
    cur_desn(:,2) = cur_desn(:,2) + nextstep(L_DESN+1:end) .* taper_desn.';
end

save("datafinal.mat");

function grad_fine = pessimgrad(dTs, Ts) % function is for a SINGLE CONTOUR.
    grad_fine = errfuncgrad(dTs,Ts); % i don't think there's any centering problems to be worried about here.
end

function out = errfunc(Ts)
    out = sum((Ts(1,:)/0.2 - 1).^2 + (Ts(2,:)/0.8 - 1).^2,2);
end

function grad_fine = errfuncgrad(dTs, Ts) % function is for a SINGLE CONTOUR.
    % Ts(1,:) is lo, Ts(2,:) is hi.
   	% each term a has a target b, so d(a/b - 1)^2 = 2/b(a/b - 1) da
    grad_fine = shiftdim(sum(dTs(1,:,:)/0.2 .* (Ts(1,:)/0.2 - 1) + dTs(2,:,:)/0.8 .* (Ts(2,:,:)/0.8 - 1),2),1);
end
