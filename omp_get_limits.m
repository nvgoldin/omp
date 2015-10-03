% authors: nadav goldin, amichai horvitz
% made for SSP course, huji 2015
function [ n_min, m_max ] = omp_get_limits( d, delta )
%OMP_get_limits(d, delta)
% calculates the minimal N rows needed in the Gaussian matrix
% and the maximal sparsity level in the signal, to assure full recovery
% using OMP algorithm, as presented in
% Signal Recovery From Random Measurements Via Orthogonal Matching Pursuit
% http://users.cms.caltech.edu/~jtropp/papers/TG07-Signal-Recovery.pdf

% d - signal length
% delta - probabilty
    % note that k=16 for small m's..
    k = 4;
    m_max = floor ( d / ( k.*log(d/delta))) -1 ;
    n_min = round ( k.* m_max * log( d / delta ) );


end

