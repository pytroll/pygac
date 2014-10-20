clear all;
close all;

fn='temp1.h5';

prt=hdf5read(fn,'prt');
iprt=hdf5read(fn,'iprt');
tprt=hdf5read(fn,'tprt');

hold on;

i=find(iprt==1);
plot(tprt(i),'r-');

i=find(iprt==2);
plot(tprt(i),'g-');

i=find(iprt==3);
plot(tprt(i),'b-');

i=find(iprt==4);
plot(tprt(i),'k-');
