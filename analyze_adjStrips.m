function res = analyze_adjStrips(adjss, ssbase, hist_ss, flucts, N_LINES, N_PER)
% we assume adjss is a cell(C_CONTS,1) made up of [C_OUTS,L_FINE,N_FREQS] arrays
% ssbase is [C_OUTS,N_FREQS]
% hist_ss is [N_MC,C_OUTS,N_FREQS]
% and flucts is a cell(C_CONTS,1) made up of [N_MC,L_FINE]
% and N_MC = N_LINES * N_PER

C_CONTS = size(adjss,1);
C_OUTS = size(adjss{1},1);
N_FREQS = size(adjss{1},3);


flucts_base = cell(C_CONTS,1);
for k = 1:C_CONTS
    flucts_base{k} = flucts{k}(N_PER:N_PER:end,:);
end

% slopes & quad coefficients are in logspace
adj_slopes = zeros(N_LINES,1,C_OUTS,N_FREQS);
scales = linspace(-1,1,N_PER);

ssbase2 = shiftdim(ssbase,-2);
ss_trues = zeros(N_LINES,N_PER,C_OUTS,N_FREQS);
ss_preds = zeros(N_LINES,N_PER,C_OUTS,N_FREQS) + ssbase2;

for i = 1:N_LINES
    ss_trues(i,:,:,:) = hist_ss(N_PER*(i-1)+1:N_PER*i,:,:);
    for k = 1:C_CONTS
        impact = sum(adjss{k}(:,:,:) .* flucts_base{k}(i,:),2); % C_OUTS x 1 x N_FREQS
        adj_slopes(i,:,:,:) = adj_slopes(i,:,:,:) + permute(impact,[4 2 1 3]) ./ ssbase2;
    end
    ss_preds(i,:,:,:) = ssbase2 .* exp(scales .* adj_slopes(i,:,:,:));
end

q_dividend = log(ss_trues ./ ss_preds);
fit_quadco = sum(q_dividend .* scales.^2,2) ./ sum(scales.^4,2);
fit_slopes = adj_slopes + sum(q_dividend .* scales,2) ./ sum(scales.^2,2);

magbase = abs(ssbase2).^2;
maj_slopes = 2 * real(adj_slopes) .* magbase;
mag_preds = magbase + maj_slopes .* scales;
q_remainder = abs(ss_trues).^2 - mag_preds;
mit_quadco = sum(q_remainder .* scales.^2,2) ./ sum(scales.^4,2);
mit_slopes = maj_slopes + sum(q_remainder .* scales,2) ./ sum(scales.^2,2);

res = struct;
res.ssbase = ssbase2;
res.ss_trues = ss_trues;

res.adj_slopes = adj_slopes;
res.fit_slopes = fit_slopes;
res.fit_quadco = fit_quadco;
res.ss_preds = ss_preds;

res.ss_predpile = zeros(7,N_LINES,N_PER,C_OUTS,N_FREQS);
res.ss_predpile(1,:,:,:,:) = shiftdim(ss_trues,-1);
res.ss_predpile(2,:,:,:,:) = res.ss_predpile(2,:,:,:,:) + shiftdim(ssbase2,-1);
res.ss_predpile(3,:,:,:,:) = shiftdim(ss_preds,-1);
res.ss_predpile(4,:,:,:,:) = shiftdim(ssbase2 .* exp(scales .* fit_slopes),-1);
res.ss_predpile(5,:,:,:,:) = shiftdim(ssbase2 .* exp(scales.^2 .* fit_quadco),-1);
res.ss_predpile(6,:,:,:,:) = shiftdim(ssbase2 .* exp(scales .* adj_slopes + scales.^2 .* fit_quadco),-1);
res.ss_predpile(7,:,:,:,:) = shiftdim(ssbase2 .* exp(scales .* fit_slopes + scales.^2 .* fit_quadco),-1);

res.maj_slopes = maj_slopes;
res.mit_slopes = mit_slopes;
res.mit_quadco = mit_quadco;
res.mpreds = mag_preds;

res.mag_predpile = zeros(7,N_LINES,N_PER,C_OUTS,N_FREQS);
res.mag_predpile(1,:,:,:,:) = shiftdim(abs(ss_trues).^2,-1);
res.mag_predpile(2,:,:,:,:) = res.mag_predpile(2,:,:,:,:) + shiftdim(magbase,-1);
res.mag_predpile(3,:,:,:,:) = shiftdim(mag_preds,-1);
res.mag_predpile(4,:,:,:,:) = shiftdim(magbase + scales .* mit_slopes,-1);
res.mag_predpile(5,:,:,:,:) = shiftdim(magbase + scales.^2 .* mit_quadco,-1);
res.mag_predpile(6,:,:,:,:) = shiftdim(magbase + scales .* maj_slopes + scales.^2 .* mit_quadco,-1);
res.mag_predpile(7,:,:,:,:) = shiftdim(magbase + scales .* mit_slopes + scales.^2 .* mit_quadco,-1);

end