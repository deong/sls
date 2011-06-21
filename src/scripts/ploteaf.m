%
% ploteaf.m
%

spea=load('spea2ts/spea2ts.eaf');
mots=load('mots/mots.eaf');

figure;
hold on;
plot(spea(:,1),spea(:,2),'b-');
plot(mots(:,1),mots(:,2),'r--');
set(gca,'FontSize',14);
legend('SPEA2+TS','MOTS');
print -depsc2 eaf.eps
