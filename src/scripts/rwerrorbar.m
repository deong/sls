%
% rwerrorbar.m
%

% compute the average random walk statistics
rw=zeros(1000,61,2);
rwbo=zeros(1000,61,2);
f=fopen('landscape/rw.out');
for i=1:1000
    for j=1:61
        line=fgets(f);
        D=sscanf(line,'%f')';
        rw(i,j,:)=[D(1),D(2)];
    end
    % eat the blank line
    line=fgets(f);
end
fclose(f);
rwmn=mean(rw);
rwsd=std(rw);

% compute the average random walk between optima statistics
f=fopen('landscape/rwbo.out');
for i=1:1000
    for j=1:61
        line=fgets(f);
        D=sscanf(line,'%f')';
        rwbo(i,j,:)=[D(1),D(2)];
    end
    line=fgets(f);
end
fclose(f);
rwbomn=mean(rwbo);
rwbosd=std(rwbo);

% draw the errorbars for the random walk statistics
figure;
hold on;
errorbar(rwmn(:,:,1),rwsd(:,:,1),'b-');
errorbar(rwbomn(:,:,1),rwbosd(:,:,1),'r--');
set(gca,'FontSize',14);
legend('Random Walk','Random Walk Between Optima','Location','SouthEast');
xlabel('Length of Random Walk');
ylabel('Swap Distance to Nearest Optimum');
print -depsc2 rwerrswapdist.eps

% draw the errorbars for the random walk statistics
figure;
hold on;
errorbar(rwmn(:,:,2),rwsd(:,:,2),'b-');
errorbar(rwbomn(:,:,2),rwbosd(:,:,2),'r--');
set(gca,'FontSize',14);
legend('Random Walk','Random Walk Between Optima','Location','SouthEast');
xlabel('Length of Random Walk');
ylabel('Fitness Distance to Nearest Optimum');
print -depsc2 rwerrfitdist.eps
