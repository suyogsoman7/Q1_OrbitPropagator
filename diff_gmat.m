% Define the states
gmat_end = [4.478764805588552e+03  -5.259844703521667e+03   1.182808046816474e+03   5.651909595207991e+00   5.014110595976435e+00   1.112431091792479e+00];
rkf78_end = [4539.4607371630600028	-5209.7523015010428935	1195.0685478367861379	5.5938825491440403	5.0759900502297768	1.0972179907412563];


% Calculate the differences element-wise
differences = gmat_end - rkf78_end;

% Calculate the percent error
percent_error = (abs(differences) ./ abs(gmat_end)) * 100;

% Display the results
disp('Difference between gmat_end and rkf78_end:');
disp(differences);

disp('Percent error (with gmat_end as the reference):');
disp(percent_error);
