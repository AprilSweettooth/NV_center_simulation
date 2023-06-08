tic
h = NV_cube(1000, 100); % increase lattice size
timeElapsed = toc;
disp('time for execution');
disp(timeElapsed);
h = h((isfinite(h)));
x = histogram(h,30);
xlabel('total coupling');
ylabel('Count');
title('Histogram of NV center coupling');
% [n,x]=hist(h);
% bar(x,n,0.5)
% x.BinWidth = 0.1;
histfit(h);
pd = fitdist(h.','Normal');
disp(pd);       %  95% confidence intervals