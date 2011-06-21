%
% plotrw.m
%

% compute the average random walk statistics
rw=zeros(61,4);
rwbo=zeros(61,4);
f=fopen('landscape/rw.out');
for i=1:1000
    for j=1:61
        line=fgets(f);
        D=sscanf(line,'%f')';
        rw(j,:)=rw(j,:)+D;
    end
    % eat the blank line
    line=fgets(f);
end
fclose(f);
rwmn=rw/1000;

% compute the average random walk between optima statistics
f=fopen('landscape/rwbo.out');
for i=1:1000
    for j=1:61
        line=fgets(f);
        D=sscanf(line,'%f')';
        rwbo(j,:)=rwbo(j,:)+D;
    end
    line=fgets(f);
end
fclose(f);
rwbomn=rwbo/1000;

% get the first 200 points from the random walk
f=fopen('landscape/rw.out');
d1=[];
for i=1:200
    for j=1:61
        line=fgets(f);
        D=sscanf(line,'%f');
        d1=[d1;D'];
    end
    line=fgets(f);
end
fclose(f);

% get the first 200 points from the random walk between optima
f=fopen('landscape/rwbo.out');
d2=[];
for i=1:200
    for j=1:61
        line=fgets(f);
        D=sscanf(line,'%f');
        d2=[d2;D'];
    end
    line=fgets(f);
end
fclose(f);
    
% draw a fdc scatter plot for the random walks
figure;
hold on;
plot(d1(1:6100,1),d1(1:6100,2),'bx');
plot(d2(1:6100,1),d2(1:6100,2),'ro');
set(gca,'FontSize',14);
legend('Random Walk','Random Walk Between Optima','Location','SouthEast');
curraxis=axis;
axis([0,60,curraxis(3),curraxis(4)]);
xlabel('Swap Distance to Nearest Optima');
ylabel('Fitness Distance to Nearest Optima');
print -depsc2 fdcrw.eps

% draw the plot of genotypic distance along the walk
figure;
hold on;
plot(rw(:,1),'k-');
plot(rwbo(:,1),'k--');
set(gca,'FontSize',14);
legend('Random Walk','Random Walk Between Optima');
xlabel('Steps in Random Walk');
ylabel('Swap Distance to Nearest Pareto Optimal Solution');
print -depsc2 swapdistrw.eps

% draw the plot of fitness distance along the walk
figure;
hold on;
plot(rw(:,2),'k-');
plot(rwbo(:,2),'k--');
set(gca,'FontSize',14);
legend('Random Walk','Random Walk Between Optima');
xlabel('Steps in Random Walk');
ylabel('Fitness Distance to Nearest Pareto Optimal Solution');
print -depsc2 fitdistrw.eps
