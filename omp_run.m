% authors: nadav goldin, amichai horvitz
% made for SSP course, huji 2015
function [ hit, miss, error, sparisity_level ] = omp_run( d, delta, n, m, noise_level, verbose)
%omp_run( d, delta, n, m, noise_level, verbose)
%   Detailed explanation goes here
% d - vector length
% delta - probability for maximal recovery needed
% n - Gaussian matrix rows number
% m - sparsity level of the signal (l_0 norm)
% noise_level - adds noise to the signal, 0 - none
% verbose - if true will print intermediate results

if verbose == 0
    fd = fopen('NUL:');
else
    fd = 1;
end

if n < m
    error('number of measurements must be bigger than sparsity level');
end


% calculate sparisity_level under limits considerations
%max_sparsity = floor ( d / ( k.*log(d/delta))) -1;


[ theoretical_n, max_sparsity ] = omp_get_limits( d, delta );
%target_sparse = floor( max_sparse_per .* max_sparsity);
target_sparse=m;
x = zeros (1, d);
x=x';

zeros_achieved=0;
while zeros_achieved ~= target_sparse
    j = randi(d);
    
    if x(j) == 0
        zeros_achieved = zeros_achieved +1;
    end
    x(j) = (-1).^(randi(2) - 1) * rand() * randi(100000);
end
sparisity_level=nnz(x);


fprintf(fd, 'Begin build matrix input\n');
fprintf(fd, '------------------------\n');
fprintf(fd, 'vector length: %d\n', d);
fprintf(fd, 'sparisity_level/max_sparsity: %d/%d\n', sparisity_level, max_sparsity);
fprintf(fd, 'n/theoretical n: %d/%d\n', n, theoretical_n);
tic
matrix=randn(n,d);
b=matrix*x;
if noise_level ~= 0
    sigma = noise_level.*norm(b)/sqrt(n);
    z= sigma*randn(n,1);
    b=b+noise_level.*z;
end
fprintf(fd, 'Elasped time %d\n', toc);
fprintf(fd, 'End build matrix input\n');


fprintf(fd, '\nBegin OMP\n');
tic
[ x_hat , indexSet, values, targetMatrix ] = OMP( matrix, b, sparisity_level );
fprintf(fd, 'Elasped time %d\n', toc);

% calculate hit/miss rate and errors
x_nnz = sort(find(x));
x_hat_nnz = sort(find(x_hat));
z = eq(x_nnz, x_hat_nnz);
hit=sum (z(:) == 1 );
miss=sum (z(:) == 0 );
error=norm(x-x_hat)/norm(x);
fprintf(fd, '\nEnd OMP\n');
fprintf(fd,'-------\n');
fprintf(fd, 'Elasped time %d\n', toc);
fprintf(fd, 'hit: %d\n', hit);
fprintf(fd, 'miss: %d\n', miss);
fprintf(fd, 'error: %d\n', error);

if verbose == 0
    fclose('all');
else
    fclose(1);
end


end

