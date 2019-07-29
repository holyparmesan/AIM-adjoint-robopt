function res = analyze_adjStrips_preddelts(hist_adjss, hist_ss, flucts, N_LINES, N_PER, varargin)
% we assume hist_adjss is a cell(C_CONTS,1) made up of [N_MC+1,C_OUTS,L_FINE,N_FREQS] arrays
% and hist_ss is an [N_MC+1,C_OUTS,N_FREQS] array
% and flucts is a cell(C_CONTS,1) made up of [N_MC+1,L_FINE]
% and N_MC = N_LINES * N_PER
% and signs is [C_CONTS], but it should have been taken care of already

C_CONTS = size(hist_adjss,1);
C_OUTS = size(hist_adjss{1},2);
N_FREQS = size(hist_adjss{1},4);

if numel(varargin) > 0
    signs = varargin{1};
else
    signs = ones(C_CONTS,1);
end

dflucts = cell(C_CONTS,1);
for k = 1:C_CONTS
    dflucts{k} = flucts{k}(3:N_PER:end,:) - flucts{k}(2:N_PER:end,:);
end

ssbase = shiftdim(hist_ss(1,:,:),-1);
preddelts = zeros(N_LINES,N_PER,C_OUTS,N_FREQS);
truedelts = zeros(N_LINES,N_PER,C_OUTS,N_FREQS);
ss_trues = zeros(N_LINES,N_PER,C_OUTS,N_FREQS);
ss_preds = zeros(N_LINES,N_PER,C_OUTS,N_FREQS) + ssbase;
ss_lpreds = zeros(N_LINES,N_PER,C_OUTS,N_FREQS) + ssbase;

for i = 1:N_LINES
    ss_trues(i,:,:,:) = hist_ss(N_PER*(i-1)+2:N_PER*i+1,:,:);
    truedelts(i,1,:,:) = ss_trues(i,2,:,:) - ss_trues(i,1,:,:);
    truedelts(i,2:end-1,:,:) = (ss_trues(i,3:end,:,:) - ss_trues(i,1:end-2,:,:))/2;
    truedelts(i,end,:,:) = ss_trues(i,end,:,:) - ss_trues(i,end-1,:,:);
    for k = 1:C_CONTS
        impact = hist_adjss{k}(N_PER*(i-1)+2:N_PER*i+1,:,:,:) .* shiftdim(dflucts{k}(i,:),-1);
        preddelts(i,:,:,:) = preddelts(i,:,:,:) + signs(k) * permute(sum(impact,3),[1 2 4 3]);
        impact_base = hist_adjss{k}(1,:,:,:) .* permute(flucts{k}(N_PER*(i-1)+2:N_PER*i+1,:),[1 3 2]);
        ss_preds(i,:,:,:) = ss_preds(i,:,:,:) + signs(k) * permute(sum(impact_base,3),[1 2 4 3]);
        ss_lpreds(i,:,:,:) = ss_lpreds(i,:,:,:) .* exp(signs(k) * permute(sum(impact_base,3),[1 2 4 3]) ./ ssbase);
    end
end

res = struct;
res.preddelts = preddelts;
res.truedelts = truedelts;
res.ss_lpreds = ss_lpreds;
res.ss_preds = ss_preds;
res.ss_trues = ss_trues;
res.ssbase = ssbase;

q_remainder = ss_trues - ss_preds;
res.quad_coeffs_ind = sum(q_remainder .* linspace(-1,1,N_PER).^2,2) ./ sum(linspace(-1,1,N_PER).^4,2);
res.quad_coeffs_glb = mean(res.quad_coeffs_ind,1);
res.ss_qpreds = ss_lpreds + res.quad_coeffs_glb .* linspace(-1,1,N_PER).^2;

q_dividend = log(ss_trues ./ ss_lpreds);
qlcoffs = sum(q_dividend .* linspace(-1,1,N_PER).^2,2) ./ sum(linspace(-1,1,N_PER).^4,2);
res.ss_qlpreds = ss_lpreds .* exp(linspace(-1,1,N_PER).^2 .* mean(qlcoffs,1));

res.m2err = zeros(3,N_PER,C_OUTS);
res.m2err(1,:,:) = mean(abs(res.ss_trues - res.ssbase).^2);
res.m2err(2,:,:) = mean(abs(res.ss_trues - res.ss_preds).^2);
res.m2err(3,:,:) = mean(abs(res.ss_trues - res.ss_qpreds).^2);

end

function ret = fliprms(arr)
    hlen2 = size(arr,2) / 2;
    avgmags2 = mean(abs(arr).^2,1);
    ret = sqrt((avgmags2(:,hlen2+1:end,:,:) + avgmags2(:,hlen2:-1:1,:,:)) / 2);
end