vecSize     = 1000; % Size of the vector in the 'original' input space
nbMeasures  = 200; % size of the measurement vector

A = randn( nbMeasures, vecSize )/sqrt(nbMeasures);
x = zeros(vecS, 1);
s = 25;
x(1:s) = rand(s,1);
y = A*x;
[xStar1,~,~,~,~,~] = htp_back(y,A,s);
res = norm(x - xStar1);
disp(['Using original htp, the result is: ', num2str(res)]);
[xStar2,~,~,~,~,~] = htp(y,A,s);
res = norm(x - xStar2);
disp(['Using modified htp, the result is: ', num2str(res)]);