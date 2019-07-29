function res = analyze_adjStrips_reindex(hists_adj, hist_ssa, adjind, rhist_ss, rflucts, fwdind, N_LINES, N_PER, varargin)

C_CONTS = size(hists_adj,1);
adjstr = cell(C_CONTS,1);
flucts = cell(C_CONTS,1);
for k = 1:C_CONTS
    adjstr{k} = shiftdim(hists_adj{k}(adjind,:,:,:),1); % [C_OUTS,L_FINE,N_FREQS]
    if numel(varargin) > 0 % first arg is flucts to be subtracted, second arg is index to use for that flucts set
        if numel(varargin) < 2
            flucts{k} = rflucts{k}(fwdind,:) - varargin{1}{k};
        else
            flucts{k} = rflucts{k}(fwdind,:) - varargin{1}{k}(varargin{2},:);
        end
    else
        flucts{k} = rflucts{k}(fwdind,:);
    end
end
ssbase = shiftdim(hist_ssa(adjind,:,:),1);
hist_ss = rhist_ss(fwdind,:,:);

res = analyze_adjStrips(adjstr, ssbase, hist_ss, flucts, N_LINES, N_PER);

end