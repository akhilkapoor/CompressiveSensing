vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

A = randn( nbMeasures, vecSize );
x = zeros(vecS, 1);
s = 25;
x(1:s) = rand(s,1);
y = A*x;
[xStar1,~,~,~,~,~] = omp_sf(y,A);
res = norm(x - xStar1);
disp(['Using original omp_sf, the result is: ', num2str(res)]);
[xStar2,~,~,~,~,~] = omp(y,A);
res = norm(x - xStar2);
disp(['Using modified omp in sf mode, the result is: ', num2str(res)]);
[xStar1,~,~,~,~,~] = omp_fast(y,A);
res = norm(x - xStar1);
disp(['Using original omp_fast, the result is: ', num2str(res)]);
[xStar2,~,~,~,~,~] = omp(y,A,'fast',true);
res = norm(x - xStar2);
disp(['Using modified omp in fast mode, the result is: ', num2str(res)]);