[pdelt, tdelt, pss, tss] = analyze_adjStrips(hist_adjss_int8, hist_ss8, flucts, N_LINES, N_PER);

qss = tss - pss;
quad_ind = sum(qss .* linspace(-1,1,N_PER).^2,2) ./ sum(linspace(-1,1,N_PER).^4,2);
quad_glb = mean(quad_ind,1);
ess = qss - quad_glb .* linspace(-1,1,N_PER).^2;
ess_ind = qss - quad_ind .* linspace(-1,1,N_PER).^2;

ssbase = shiftdim(hist_ss8(1,:,:),-1);
errsn_const = fliprms(tss ./ ssbase - 1);
errsn_linear = fliprms(qss ./ ssbase);
errsn_quad_glb = fliprms(ess ./ ssbase);
errsn_quad_ind = fliprms(ess_ind ./ ssbase);

lins = fliprms(linspace(-1,1,N_PER));

function ret = fliprms(arr)
    hlen2 = size(arr,2) / 2;
    avgmags2 = mean(abs(arr).^2,1);
    ret = sqrt((avgmags2(:,hlen2+1:end,:,:) + avgmags2(:,hlen2:-1:1,:,:)) / 2);
end