function avg = errfigmaker(predpile, ROUGH_AMPL)
% predpile is [7 x N_LINES x N_PER x C_OUT x N_FREQS].
% The seven datasets must be arranged as below.
% we flip-stack the N_PER dimention into the N_LINES dimension square, 
% then square, average over N_LINES and N_FREQS, divide, plot

N_PER = size(predpile,3);
errors = predpile - predpile(1,:,:,:,:);
combined = permute(cat(2,errors(:,:,N_PER/2:-1:1,:,:),errors(:,:,N_PER/2+1:end,:,:)),[1 3 4 5 2]);
restacked = combined(:,:,:,:); % merge the last two dimensions
avg = mean(abs(restacked).^2,4);
scales = linspace(-1,1,N_PER);
ampls = ROUGH_AMPL * 1e9 * scales(N_PER/2+1:end);

figure; hold on;
plot(ampls,avg(3,:) ./ avg(2,:),'o--','DisplayName','Adjoint linear prediction');
plot(ampls,avg(4,:) ./ avg(2,:),'o-','DisplayName','Best-fit linear prediction');
plot(ampls,avg(5,:) ./ avg(2,:),'*-','DisplayName','Best-fit quadratic prediction');
plot(ampls,avg(6,:) ./ avg(2,:),'^--','DisplayName','Adjoint linear+quadratic prediction');
plot(ampls,avg(7,:) ./ avg(2,:),'^-','DisplayName','Best-fit linear+quadratic prediction');
ax = gca;
ax.FontSize=11;
xlabel('Amplitude (nm)','fontsize',14);
ylabel('Fraction of variance unexplained','fontsize',14);
ylim([0 1]);
yticks([0 0.2 0.4 0.6 0.8 1.0]);

end