load('ropt_reg_0520.mat')

pred_erfs = hist_erfs(1:end-1,:) + sum(permute(hist_desns(2:end,:) - hist_desns(1:end-1,:),[1 3 2]) .* hist_grads(1:end-1,:,:),3);

plot(hist_erfs,'o-')
hold on;
for i = 1:29
for j = 1:4
plot([i; i+1],[hist_erfs(i,j); pred_erfs(i,j)],'o--','color','k');
end
end