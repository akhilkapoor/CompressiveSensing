addpath('./../algorithms');

vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

A = randn( nbMeasures, vecSize ) / sqrt(nbMeasures);
x = zeros(vecSize, 1);
s = 25;
x(1:s) = rand(s,1);
y = A*x;

n = 1;
disp(['Test ', num2str(n)]); n = n+1;
[xStar1,~,~,~,~] = CoSaMP(y, A, s);
res = norm(x - xStar1);
disp(['Using CoSaMP, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar2,~,~,~,~] = CoSaMP(y, A, s, {'maxnbiter',s*2});
res = norm(x - xStar2);
disp(['Using CoSaMP, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar3,~,~,~,~] = CoSaMP(y, A, s-10, {'maxnbiter', s*2, 'adds', s*2});
res = norm(x - xStar3);
disp(['Using CoSaMP, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar4,~,~,~,~] = CoSaMP(y, A, s+5, {'HSS', true});
res = norm(x - xStar4);
disp(['Using CoSaMP, the result is: ', num2str(res)]);
