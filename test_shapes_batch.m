lambda_list = [0 .5 1 10 20 100];
lambda_list = [0 10 50 100 500];
lambda_list = [200 100 50  0]; % 
Wfixed = [1/2;1/2];
lambda_list = [0 10 20 40 60 80 100 125 150 175 200 250 300]; % 

for il = 1:length(lambda_list)
    lambda = lambda_list(il);
    disp(['----- lambda=' num2str(lambda) ' -----']);
    test_shapes;
end