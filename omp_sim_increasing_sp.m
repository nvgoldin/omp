% authors: nadav goldin, amichai horvitz
% made for SSP course, huji 2015
% Signal Recovery simulation using Orthogonal Matching Pursuit
% attempts to recover signal in length d, with sparsity level
% of m={1, 2.., 2^(log(d) -1)}
% make sure parpool is set to make things work quickly

% parameters
% signal length
d=2^10;
% after error_thrsehold times in a row both the support was recovered
% and the error is less than 1e-5 the loop is stopped for the given m
error_thrsehold=20;
% probabilty - not used
delta=0.01;


tic


field_names = {'N_vector','hits', 'errors', 'header' ,'sparsity_level'}; 
empty_cells = repmat(cell(1), 1, numel(field_names));
entries = { field_names{:} ; empty_cells{:}};
scatterz = repmat(struct(entries{:}), nextpow2(d) -1,1 );


parfor sparsity_power= 1:(nextpow2(d) - 1)
    sparsity_level=2.^(sparsity_power-1);
    last_index = 0;
    hits=zeros(1, d);
    errors=zeros(1, d);
    errors(errors==0)=nan;
    hits(hits==0)=nan;
    i = sparsity_level;
    no_errors=0;
    last_error=1;
    while i < d && no_errors < error_thrsehold
        verbose = 0;
        [hit, miss, errors(i), ~ ] = omp_run(d, delta, i, sparsity_level, 0, verbose);
        hits(i)= hit;
        if (abs ( errors(i) ) < 1e-5 && hits(i) == sparsity_level )
            if ( last_error == 0 )
                no_errors = no_errors + 1;
                last_index = i;
            else
                no_errors = 1;
                last_error = 0;
            end
        else
            last_error = 1;
            no_errors = 0;
        end

        i = i + 1;
    end
    

    N_vector=1:1:d;
    scatterz(sparsity_power).N_vector = N_vector(1:last_index);
    scatterz(sparsity_power).errors = errors(1:last_index);
    scatterz(sparsity_power).hits = hits(1:last_index);
    scatterz(sparsity_power).sparsity_level = sparsity_level;

    head_begin = sprintf(', d=%d, m=%d', d, sparsity_level);
    if last_index == 0
        head_end = sprintf('Reconstruction FAILED up to N=%i', i);
    else
        head_end = sprintf('Full reconstruction occured %d times in a row at N=%d', error_thrsehold, last_index);
    end
    header = strcat( head_end, head_begin);
    scatterz(sparsity_power).header = header;

end
if length(scatterz) > 6
    vec = [length(scatterz)-5:1:length(scatterz)];
else
    vec = 1:length(scatterz);
end
% print errors
i = 1;
for j = vec
    figure(1)
    subplot(3,2,i)
    scatter(scatterz(j).N_vector, scatterz(j).errors);
    xlabel('N - number of measurements');    
    ylabel('Error');
    grid on;
    title(scatterz(j).header);
    i = i + 1;
    
end

% print support hits
i = 1;
for j = vec
    figure(2)
    subplot(3,2,i)
    scatter(scatterz(j).N_vector, scatterz(j).hits, '+', 'b');
    hline = refline([0 scatterz(j).sparsity_level]);
    set(hline,'LineStyle',':');
    set(hline, 'Color', 'red');
    xlabel('N - number of measurements');
    ylabel('support hits');
    grid on;
    title(scatterz(j).header);
    i=i+1;

end



toc