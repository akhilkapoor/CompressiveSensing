addpath('./../algorithms');

vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

A = randn( nbMeasures, vecSize );
x = zeros(vecSize, 1);
s = 25;
x(1:s) = rand(s,1);
y = A*x;

n = 1;
disp(['Test ', num2str(n)]); n = n+1;
[xStar1,~,~,~,~,~] = omp_sf(y,A);
res = norm(x - xStar1);
disp(['Using original omp_sf, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar2,~,~,~,~,~] = omp(y,A);
res = norm(x - xStar2);
disp(['Using modified omp in sf mode, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar3,~,~,~,~,~] = omp_fast(y,A);
res = norm(x - xStar3);
disp(['Using original omp_fast, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar4,~,~,~,~,~] = omp(y,A,{});
res = norm(x - xStar4);
disp(['Using modified omp in fast mode, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar5,~,~,~,~,~] = omp(y,A,{'fast',true});
res = norm(x - xStar5);
disp(['Using modified omp in fast mode, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar6,~,~,~,~,~] = omp(y,A,{'maxNbIter'});
res = norm(x - xStar6);
disp(['Using modified htp, the result is: ', num2str(res)]);

disp(['Test ', num2str(n)]); n = n+1;
[xStar7,~,~,~,~,~] = omp(y,A, 'fast', false );
res = norm(x - xStar7);
disp(['Using modified omp in sf mode, the result is: ', num2str(res)]);
