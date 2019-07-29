load('ropt_reg_0520.mat')
hist_conts = cell(4,1);
hist_conts{1} = shiftdim(baseconts{1},-1) + (hist_desns * splinmat) .* shiftdim(norms_fine{1},-1);
hist_conts{2} = shiftdim(baseconts{2},-1) + zeros(30,1);
hist_conts{3} = shiftdim(baseconts{3},-1) + zeros(30,1);
hist_conts{4} = shiftdim(baseconts{4},-1) - (hist_desns * splinmat) .* shiftdim(norms_fine{4},-1);
subplot(2,4,1);
poly = shiftdim([hist_conts{1}(1,:,:) flip(hist_conts{2}(1,:,:),2) hist_conts{3}(1,:,:) flip(hist_conts{4}(1,:,:),2)],1);
fill(poly(:,1),poly(:,2),'r');
xlim([-1.5e-6, 3e-6]);
axis off;
title("Iter 1");
subplot(2,4,2);
poly = shiftdim([hist_conts{1}(5,:,:) flip(hist_conts{2}(5,:,:),2) hist_conts{3}(5,:,:) flip(hist_conts{4}(5,:,:),2)],1);
fill(poly(:,1),poly(:,2),'r');
xlim([-1.5e-6, 3e-6]);
axis off;
title("Iter 5");
it = 15;
subplot(2,4,3);
poly = shiftdim([hist_conts{1}(it,:,:) flip(hist_conts{2}(it,:,:),2) hist_conts{3}(it,:,:) flip(hist_conts{4}(it,:,:),2)],1);
fill(poly(:,1),poly(:,2),'r');
xlim([-1.5e-6, 3e-6]);
axis off;
title("Iter 15");
subplot(2,4,4);
it = 30;
poly = shiftdim([hist_conts{1}(it,:,:) flip(hist_conts{2}(it,:,:),2) hist_conts{3}(it,:,:) flip(hist_conts{4}(it,:,:),2)],1);
fill(poly(:,1),poly(:,2),'r');
xlim([-1.5e-6, 3e-6]);
axis off;
title("Iter 30");