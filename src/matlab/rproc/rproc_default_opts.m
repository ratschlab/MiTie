function opts = rproc_default_opts()

opts.maxjobs = 2500;
opts.resubmit = 2;
opts.priority = 17;
opts.mem_req_resubmit = [10000 20000];
opts.time_req_resubmit = [1e6 1e6];
%opts.no_result_file=1;
opts.identifier = 'MIP';

