for i = 1:N_LINES
for k = 1:4
preddelts(i,:,:,:) = preddelts(i,:,:,:) + signs(k) * shiftdim(squeeze(sum(hist_adjss{k}(N_PER*(i-1)+2:N_PER*i+1,:,:,:) .* shiftdim(flucts_base{k}(i,:),-1) / (N_PER/2.0), 3)),-1);
end
end

sslo = zeros(N_LINES,N_PER,N_FREQS);
sshi = zeros(N_LINES,N_PER,N_FREQS);
for i = 1:N_LINES
sslo(i,:,:) = squeeze(hist_ss(N_PER*(i-1)+2:N_PER*i+1,1,:));
sshi(i,:,:) = squeeze(hist_ss(N_PER*(i-1)+2:N_PER*i+1,2,:));
end

truedelts = zeros(N_LINES,N_PER,2,N_FREQS);
for i = 1:N_LINES
truedelts(:,1,1,:) = sslo(:,2,:) - sslo(:,1,:);
truedelts(:,end,1,:) = sslo(:,end,:) - sslo(:,end-1,:);
truedelts(:,2:end-1,1,:) = (sslo(:,3:end,:) - sslo(:,1:end-2,:))/2;
end
for i = 1:N_LINES
truedelts(:,1,2,:) = sshi(:,2,:) - sshi(:,1,:);
truedelts(:,end,2,:) = sshi(:,end,:) - sshi(:,end-1,:);
truedelts(:,2:end-1,2,:) = (sshi(:,3:end,:) - sshi(:,1:end-2,:))/2;
end

for i = 1:N_LINES
hold on; plot(squeeze(preddelts(i,:,2,1)).','o-'); pause; plot(smooth(squeeze(truedelts(i,:,2,1))).','o-'); pause;
end

for i = 1:N_LINES
hold on; plot(squeeze(preddelts(i,:,2,3) ./ sshi(i,:,3)).','o-'); pause; plot(smooth(squeeze(truedelts(i,:,2,3) ./ sshi(i,:,3))).','o-'); pause;
end