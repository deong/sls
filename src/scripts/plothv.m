%
% plothv.m
%
% plot the hypervolume of the three algorithms
%

spea2ts=load('spea2ts/trial1.hv');
spea2ts=spea2ts(:,3);
emoeats=load('emoeats/trial1.hv');
emoeats=emoeats(:,3);
emoeatpts=load('emoeatpts/trial1.hv');
emoeatpts=emoeatpts(:,3);

spea2tseval=load('spea2ts/trial1.hv');
spea2tseval=spea2tseval(:,1);
emoeatseval=load('emoeats/trial1.hv');
emoeatseval=emoeatseval(:,1);
emoeatptseval=load('emoeatpts/trial1.hv');
emoeatptseval=emoeatptseval(:,1);

spea2tssize=size(spea2ts,1);
emoeatssize=size(emoeats,1);
emoeatptssize=size(emoeatpts,1);

for i=2:10
    file=sprintf('spea2ts/trial%d.hv',i);
    tmp=load(file);
    
    % check to see if tmp is too big
    if(size(tmp,1)>spea2tssize)
        tmp=tmp(1:spea2tssize,:);
    end
    if(size(tmp,1)<spea2tssize)
        spea2tssize=size(tmp,1);
        spea2ts=spea2ts(1:spea2tssize,:);
        spea2tseval=spea2tseval(1:spea2tssize,:);
    end
    spea2ts=[spea2ts,tmp(:,3)];
    spea2tseval=[spea2tseval,tmp(:,1)];
end

for i=2:10
    file=sprintf('emoeats/trial%d.hv',i);
    tmp=load(file);
    
    % check to see if tmp is too big
    if(size(tmp,1)>emoeatssize)
        tmp=tmp(1:emoeatssize,:);
    end
    if(size(tmp,1)<emoeatssize)
        emoeatssize=size(tmp,1);
        emoeats=emoeats(1:emoeatssize,:);
        emoeatseval=emoeatseval(1:emoeatssize,:);
    end
    emoeats = [emoeats,tmp(:,3)];
    emoeatseval=[emoeatseval,tmp(:,1)];
end

for i=2:10
    file=sprintf('emoeatpts/trial%d.hv',i);
    tmp=load(file);
    
    % check to see if tmp is too big
    if(size(tmp,1)>emoeatptssize)
        tmp=tmp(1:emoeatptssize,:);
    end
    if(size(tmp,1)<emoeatptssize)
        emoeatptssize=size(tmp,1);
        emoeatpts=emoeatpts(1:emoeatptssize,:);
        emoeatptseval=emoeatptseval(1:emoeatptssize,:);
    end
    emoeatpts = [emoeatpts,tmp(:,3)];
    emoeatptseval=[emoeatptseval,tmp(:,1)];
end

% transpose the matrices for use with mean,std
spea2ts=spea2ts';
emoeats=emoeats';
emoeatpts=emoeatpts';
spea2tseval=spea2tseval';
emoeatseval=emoeatseval';
emoeatptseval=emoeatptseval';

% compute mean values
spea2tsmn=mean(spea2ts);
emoeatsmn=mean(emoeats);
emoeatptsmn=mean(emoeatpts);
spea2tsevalmn=mean(spea2tseval);
emoeatsevalmn=mean(emoeatseval);
emoeatptsevalmn=mean(emoeatptseval);

% and standard deviations
spea2tsdv=std(spea2ts);
emoeatsdv=std(emoeats);
emoeatptsdv=std(emoeatpts);

% construct approximations using fewer data points
% to make the error bars less dense
numpoints=size(spea2ts,2);
bigpoints=size(emoeats,2);
%interval=floor(bigpoints/numpoints);
interval=bigpoints/numpoints;
d1=zeros(1,numpoints);
d2=zeros(1,numpoints);
d3=zeros(1,numpoints);
for i=1:numpoints
    d1(i)=emoeatsevalmn(floor(i*interval));
    d2(i)=emoeatsmn(floor(i*interval));
    d3(i)=emoeatsdv(floor(i*interval));
end
bigpoints=size(emoeatpts,2);
interval=bigpoints/numpoints;
d4=zeros(1,numpoints);
d5=zeros(1,numpoints);
d6=zeros(1,numpoints);
for i=1:numpoints
    d4(i)=emoeatptsevalmn(floor(i*interval));
    d5(i)=emoeatptsmn(floor(i*interval));
    d6(i)=emoeatptsdv(floor(i*interval));
end

% and plot the curves with error bars
hold on;
errorbar(spea2tsevalmn,spea2tsmn,spea2tsdv,'r');
errorbar(d1,d2,d3,'b');
errorbar(d4,d5,d6,'c');

% compute the axis boundaries
minx=min([min(min(spea2tseval)),...
          min(min(emoeatseval)),...
          min(min(emoeatptseval))]);
maxx=max([max(max(spea2tseval)),...
          max(max(emoeatseval)),...
          max(max(emoeatptseval))]);
curraxis=axis;
axis([minx,maxx,curraxis(3),curraxis(4)]);

% draw the legend
legend('SPEA2+TS','\epsilon-MOEA+TS','\epsilon-MOEA+TPTS','Location','SouthEast');

% and the title
title('Comparison of Hybrid Multiobjective Evolutionary Algorithms');

% and the labels
xlabel('Evaluations');
ylabel('Hypervolume');

% and save the plot
print -depsc2 hypervolume.eps
