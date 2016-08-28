lambda_list = 0:.02:.4;

for il = 1:length(lambda_list)
    lambda = lambda_list(il);
    disp(['----- lambda=' num2str(lambda) ' -----']);
    test_meg;
end