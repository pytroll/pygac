clear all;
close all;

fn='temp1.h5';

prt=hdf5read(fn,'prt');
iprt=hdf5read(fn,'iprt');
lines=hdf5read(fn,'lines');
counts=hdf5read(fn,'counts');
space=hdf5read(fn,'space');
nlin=hdf5read(fn,'nlin');
ncor=hdf5read(fn,'ncor');
tse=hdf5read(fn,'tse');

hold on;

i=find(iprt==1);
plot(lines(i),prt(i),'r-');

i=find(iprt==2);
plot(lines(i),prt(i),'g-');

i=find(iprt==3);
plot(lines(i),prt(i),'b-');

i=find(iprt==4);
plot(lines(i),prt(i),'k-');


plot(lines(1:length(lines)-1),diff(lines)+230,'c-');


i=find(iprt==0);
plot(lines(i),prt(i)+250,'m-');


figure;

subplot(4,1,1);
imagesc(counts-space); colorbar;

subplot(4,1,2);
imagesc(nlin); colorbar;

subplot(4,1,3);
imagesc(ncor); colorbar;

subplot(4,1,4);
imagesc(tse); colorbar;

