clear; clc, close('all')

%% Laterality index
load laterality_CRM
figure(1),
subplot 121
hbar = bar(par_crm);
view([90 -90])
L = length(par_crm);
for i=1:L
    S{i} = num2str(i);
end
% set(gca,'Xtick',idx,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
title('Word-recognition','fontsize',16)
% legend([{'fMRI-tValue'};{'MEG'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);

DataArray = par_crm;
Colors = [0.4922    0.0039    0.9063;0.9922    0.8672    0.0039];
figure(2)
subplot 121
[xPositions, yPositions, ~, ~] = UnivarScatter(par_crm,'Label',{'fMRI','meg'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 16);
xlabel('Modality','FontSize', 16);
set(gca,'color','none');
title('Word-recognition','fontsize',16)
disp(['fmri: ',num2str(mean(par_crm(:,1))), '+-', num2str(std(par_crm(:,1)))]);
disp(['meg: ',num2str(mean(par_crm(:,2))), '+-', num2str(std(par_crm(:,2)))]);
set(gca,'FontName','HelveticaNeueLT Std Lt');
hold on
f = [xPositions, yPositions];
for j=1:length(f)
    line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
end


for k = 1:numel(hbar)
    set(hbar(k),'FaceColor',Colors(k,:))
end

%
load laterality_VGA
figure(1)
subplot 122
hbar = bar(par_vga);
view([90 -90])
L = length(par_vga);
for i=1:length(L)
    S{i} = num2str(i);
end
% set(gca,'Xtick',idx,'XtickLabel',S);
xlim([0,L+1]);
ylim([-1,1]);
set(gca,'color','none');
box off
xlabel('Subj');
ylabel('Laterality');
% set(gcf, 'Position', [800   500   1000  500]);
title('Verb-Generation','fontsize',16)
% legend([{'fmri-tValue'};{'meg-hubs'}]);
% set(gca,'FontName','Arial');
set(gca,'FontName','HelveticaNeueLT Std Lt');
disp(['fmri: ',num2str(mean(par_vga(:,1))), '+-', num2str(std(par_vga(:,1)))]);
disp(['meg: ',num2str(mean(par_vga(:,2))), '+-', num2str(std(par_vga(:,2)))]);
% colormap(Colors);
ylabel('Laterality','FontSize', 16);
xlabel('Subj','FontSize', 16);


for k = 1:numel(hbar)
    set(hbar(k),'FaceColor',Colors(k,:))
end

figure(2)
subplot 122
[xPositions, yPositions, Label, RangeCut] = UnivarScatter(par_vga,'Label',{'fmri','meg'},'MarkerFaceColor',Colors);
ylabel('Laterality','FontSize', 16);
set(gca,'color','none');
set(gca,'FontName','HelveticaNeueLT Std Lt');
title('Verb-Generation','fontsize',16)
xlabel('Modality','FontSize', 16);
f = [xPositions, yPositions];
for j=1:length(f)
    line([f(j,1),f(j,2)],[f(j,3),f(j,4)]);
end
