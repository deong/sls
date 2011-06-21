%
% plotfdc.m
%

data=load('landscape/fdc.out');
plot(data(1:1000,1),data(1:1000,2),'k*');
currax=axis;
axis([0,60,currax(3),currax(4)]);
print -deps fdc.eps

