function makeGreyFig(powers)

figure; hold on;
plot([0 1],[1 0],'k--','HandleVisibility','off');

ax = gca;

for i = 1:3
    for m = 1:10
        if 2*powers(i,m,1,1,5)+powers(i,m,1,2,5) > 2*powers(i,m,20,1,5)+powers(i,m,20,2,5)
%            powers(i,m,1:end,:,:) = powers(i,m,end:-1:1,:,:);
        end
    end
end

for k = [1 5 9]
    for i = 1:3
        plot(squeeze(powers(i,:,:,1,k)).',squeeze(powers(i,:,:,2,k)).','color',[0.8 0.8 0.8],'LineWidth',0.5,'HandleVisibility','off');
    end
    ax.ColorOrderIndex = 1;
    for i = 1:3
        if k == 1
            plot(squeeze(mean(powers(i,:,:,1,k),2)),squeeze(mean(powers(i,:,:,2,k),2)),'-','LineWidth',2);
        else
            plot(squeeze(mean(powers(i,:,:,1,k),2)),squeeze(mean(powers(i,:,:,2,k),2)),'-','LineWidth',2,'HandleVisibility','off');
        end
    end
end

xlim([0.1 0.3]);
ylim([0.65 0.85]);