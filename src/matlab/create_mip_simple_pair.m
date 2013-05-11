function [transcripts, weights, how, solvetime, clusters] = create_mip(observed_cov, len, admat, predef_trans, initial, terminal, C, pair_observation, insert_size, param, cov_scale, strand)

if nargin==1
	mip_paths;

	PAR = observed_cov;
	observed_cov = PAR.coverage;
	len = PAR.len;
	admat = PAR.seg_admat;
	predef_trans = PAR.predef_trans; 
	initial = PAR.initial; 
	terminal = PAR.terminal; 
	C = PAR.C;
	pair_observation = PAR.pair_list; 
	insert_size = PAR.insert_size; 
	param = PAR.param;
	cov_scale = PAR.cov_scale;
	strand = PAR.strand;
	fn_res = PAR.fn_res;
end

max_num_tr = param.number_of_additional_transcripts+length(predef_trans);
disp('calling create_mip') ;
use_pair = 0;
if nargin>=9 && ~isempty(pair_observation) && ~isempty(insert_size)
	use_pair = 1;
	tol = insert_size*param.insert_size_tol;
end
if ~exist('param', 'var')
    param.time_limit=60 ;
end ;

param.len_cutoff = min(param.len_cutoff, sum(len));

use_LP_covloss=0;
sort_novel_transcripts=1; % this reduces the run time significantly


%%% loop over multiple RNA-seq samples
%%%
%%% process each sample separately and concat 
%%% only in the end. 
num_samples = size(observed_cov, 1);
num_clusters = min(param.num_clusters, num_samples);
use_cluster = num_clusters>3;
assert(size(admat, 3)==num_samples)

admat_all = admat;

binary_var_index_all = [];
n_of_equalities = 0;

intron_conf = compute_intron_list(admat(:,:,1));

% compute the number of variables
%
% U_st : usage of segment s in transcript t                                 	(binary)
% I_t  : indicator if transcript t has weight>=0 in any sample                  (binary)
% E_str : expected coverage of segment s in transcript t in sample r                 
% W_tr  : weight of transcript t in sample r
% L_sr  : loss of segment s (deviation from the observed coverage) in sample r
% C_ctr : Spliced read expected coverage (c is a subset of sxs) in sample r
% D_cr  : deviation of expected intron coverage from observed intron cov in sample r
%
% var = [ U_st I_t E_str W_tr L_sr C_ctr D_cr]
c = size(intron_conf, 1);
t = max_num_tr;
s = size(observed_cov, 2);
r = num_samples;
nc = num_clusters;

num_var = 0;
% indices of variables in the parameter vector
U_idx = 1:s*t; 						num_var = num_var+length(U_idx);
I_idx = num_var+1:num_var+t;		num_var = num_var+length(I_idx);
E_idx = num_var+1:num_var+s*t*r;	num_var = num_var+length(E_idx);
W_idx = num_var+1:num_var+t*r;		num_var = num_var+length(W_idx);
L_idx = num_var+1:num_var+s*r*2;	num_var = num_var+length(L_idx);
C_idx = num_var+1:num_var+c*t*r;	num_var = num_var+length(C_idx);
D_idx = num_var+1:num_var+c*r*2;	num_var = num_var+length(D_idx);

if use_cluster
	M_idx = num_var+1:num_var+r*nc; num_var = num_var+length(M_idx); % cluster assignment
	m_idx = num_var+1:num_var+t*nc; num_var = num_var+length(m_idx); % cluster centroids
	K_idx = num_var+1:num_var+nc; num_var = num_var+length(K_idx);   % cluster usage
end
if use_pair
	np = size(pair_observation, 1);
	P_idx = num_var+1:num_var+t*np; num_var = num_var+length(P_idx);
	SP_idx = num_var+1:num_var+np; num_var = num_var+length(SP_idx); % slacks for pair observation
end

% Length of variables 
nU = length(U_idx);	assert(nU==s*t);
nI = length(I_idx);	assert(nI==t);
nE = length(E_idx); assert(nE==s*t*r);
nW = length(W_idx); assert(nW==t*r);
nL = length(L_idx);	assert(nL==s*r*2);
nC = length(C_idx);	assert(nC==t*c*r);
nD = length(D_idx); assert(nD==c*r*2);


%                          A1      A2    A3     A4+A5+A6 A8+A9  A10+A11+A12
expected_num_constrains = 2*s*r + s*t + 2*s*t + 3*r*s*t + 2*t + 3*r*t*c; 


%%%
% compute profile matrix: length bias 
profile_mat = ones(t, s);
if 0 % use_profile
	ll = load('/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/hg19_annotations_merged_splice_graph_expr1_max_trans.profile.mat');
	%ll = load('/fml/ag-raetsch/nobackup/projects/mip/human_sim/data_sim_500_alt25/rquant/rquant/genes_expr_weak_bias1_optimalProf_noSeqBias_2012-04-09/profiles.mat')
	if strand=='-'
		ll.profile_weights = ll.profile_weights(end:-1:1, :);
	end
	for k = 1:length(predef_trans)
		tr_len = sum(len(find(predef_trans{k})));
		prof_idx = find(ll.CFG.transcript_len_ranges(:, 1)<=tr_len&ll.CFG.transcript_len_ranges(:, 2)>=tr_len);

		cum_len = 0;
		for j = 1:s
			cpos = cum_len+predef_trans{k}(j)*len(j)/2; 
			cum_len = cum_len+predef_trans{k}(j)*len(j);
			if cpos<tr_len/2
				idx = find(ll.CFG.limits(1:end-1)<=cpos+1&ll.CFG.limits(2:end)>cpos+1);
			else
				idx = find(ll.CFG.limits(1:end-1)<=(tr_len-cpos+1)&ll.CFG.limits(2:end)>(tr_len-cpos+1));
				idx = size(ll.profile_weights, 1)-idx+1;
			end
			profile_mat(k, j) = ll.profile_weights(idx, prof_idx);
		end
	end
end

fprintf('create mip: max_num_trans:%i, segments:%i, samples:%i introns:%i num_var:%i\n', t, s, r, c, num_var);
% Equality constraints: 
% sum up the expected coverage for each segment and penalize the deviation 
% from the observed coverage using a helper variable (l1, ..., ln)
%
%| - - - - - - - - - - - - sample 1 - - - - - - - - - - - - - - - - - | - - - -
%    segment 1 |   segment 2  | ... |  segment n   |...| loss         |
% t1 t2 ... t4 | t1 t2 ... t4 | ... | t1 t2 ... tk |...| l1 l2 ... ln |
%  1  1 ...  1 |  0  0 ...  0 | ... |  0  0 ...  0 |...| -1  0 ... 0  |

% in case of a quadratic penalty:
% sum_t E_str -L_sr = O_sr
% 
% for a linear penalty:
%
%	sum_t E_str -L_sr <= O_sr
% and 
%	-sum_t E_str -L_sr <= -O_sr

A1 = zeros(s*r, num_var);
b1 = zeros(s*r, 1);
cnt = 1;
for i=1:r
	for j = 1:s
		idx = (i-1)*s*t + (j-1)*t + (1:t);
		A1(cnt, E_idx(idx)) = profile_mat(:, j)'; %sum_t E_jti
		A1(cnt, L_idx((i-1)*s+j)) = -1; % -l_ji 
		A1(cnt, L_idx((i-1)*s+j+s*r)) = 1; % -l_ji
		cnt = cnt+1;
	end
	b1((i-1)*s+1:(i-1)*s+s) = observed_cov(i, :); % O_sr
end
clear i j idx cnt 

n_of_equalities = n_of_equalities+size(A1, 1);
A1 = sparse(A1);
	

% Inequality constraints: 
	
% if no segments are selected, W_t has to be 0
% W_x<=\sum_j U_{jx}
A2 = zeros(r*t , num_var);
cnt = 1;
for i=1:r
	for x = 1:t
		A2(cnt, W_idx((i-1)*t+x)) = 1; % W_xi
		for j=1:s
			A2(cnt, U_idx((j-1)*t+x)) = -1; % -U_jx
		end
		cnt = cnt+1;
	end
end
clear i x j cnt 
b2 = zeros(r*t,1);
A2 = sparse(A2);
	
% 
% this takes the connectivity of the splice graph and makes 
% sure that segment j is followed by any of the segments 
% it is connected to in the splice graph G=(N,E)
%
% (1) U_{jx} <= \sum_{k \in \{k | (j,k)\in E\}} U_{cx} 
% (2) U_{jx} <= \sum_{k \in \{k | (k,j)\in E\}} U_{cx} 
%
% for initial (and terminal) segments make sure that there are no
% downstream (upstream) segments used
%
% terminal:
% (3) \sum_{k=j+1}^s U_kx <= (s-j) - (s-j)*U_jx
% <=> \sum_{k=j+1}^s U_kx + (s-j)*U_jx <= (s-j)
%
% initial:
% (4) \sum_{k=1}^{j-1} U_jx <= (j-1) - (j-1)*U_jx
% <=> \sum_{k=1}^{j-1} U_jx + (j-1)*U_jx <= (j-1)

if s<127
	A3 = zeros(2*t*s, num_var, 'int8');
else
	A3 = zeros(2*t*s, num_var, 'int16');
end
b3 = zeros(2*t*s, 1);
cnt = 0;
for x = 1:t
	for j = 1:s
		% get all connected nodes
		cnodes = find(admat(j,:, 1)>=-1);
		children = cnodes(cnodes>j);
		parents = cnodes(cnodes<j);
		if ~isempty(children) && ~terminal(j)
			cnt = cnt+1;
			A3(cnt, U_idx((j-1)*t+x)) = 1; % U_jx
			for k=children
				A3(cnt, U_idx((k-1)*t+x)) = -1; % -U_kx
			end
		elseif isempty(children)
			% make sure there is no downstream segment used if 
			% U_jx is used
			cnt = cnt+1;
			A3(cnt, U_idx((j-1)*t+x)) = s-j; % (s-j)U_jx
			for k=j+1:s
				A3(cnt, U_idx((k-1)*t+x)) = 1; % U_kx
			end
			b3(cnt) = s-j;
		end
		if ~isempty(parents) && ~initial(j)
			cnt = cnt+1;
			A3(cnt, U_idx((j-1)*t+x)) = 1; % U_jx
			for k=parents
				A3(cnt, U_idx((k-1)*t+x)) = -1; % -U_kx
			end
		elseif isempty(parents)
			% make sure there are no upsteam segments if 
			% U_jx is used
			cnt = cnt+1;
			A3(cnt, U_idx((j-1)*t+x)) = j-1; % U_jx
			for k=1:j-1
				A3(cnt, U_idx((k-1)*t+x)) = 1; % U_kx
			end
			b3(cnt) = j-1;
		end
	end
end
A3(cnt+1:end, :) = [];
b3(cnt+1:end) = [];
%A3 = sparse(A3);
if s<127
	A3 = int8tosparse(A3);
else
	A3 = int16tosparse(A3);
end
clear x j cnt

% E_str
% define helper variables computing the expected coverage for each segment
% The expected coverage is given by 
% E_str = U_st * W_tr
% bound them in [0, 1]

% E_str - U_st <= 0
A4 = zeros(r*s*t , num_var, 'int8');
cnt = 1;
for i=1:r
	for j = 1:s
		for k = 1:t
			pos = (i-1)*s*t+(j-1)*t+k;
			A4(cnt, E_idx(pos)) = 1; % E_st
			A4(cnt, U_idx((j-1)*t+k)) = -1; % U_st
			cnt = cnt+1;
		end
	end
end
clear i j k pos cnt
b4 =  zeros(r*s*t, 1);
%A4 = sparse(A4);
A4 = int8tosparse(A4);
	
% E_str + U_st - W_t <= 1
A5 = zeros(s*t*r , num_var, 'int8');
cnt = 1;
for i=1:r
	for j = 1:s
		for k = 1:t
			pos = (i-1)*s*t+(j-1)*t+k;
			A5(cnt, E_idx(pos)) = 1; % E_st 
			A5(cnt, U_idx((j-1)*t+k)) = 1; % U_st
			A5(cnt, W_idx((i-1)*t+k)) = -1; % -W_ki
			cnt = cnt+1;
		end
	end
end
clear i j k pos cnt
b5 = ones(s*t*r, 1);
%A5 = sparse(A5);
A5 = int8tosparse(A5);
	
% -E_st + U_st + W_t <= 1
A6 = zeros(s*t*r , num_var, 'int8');
cnt = 1;
for i=1:r
	for j = 1:s
		for k = 1:t
			pos = (i-1)*s*t+(j-1)*t+k;
			A6(cnt, E_idx(pos)) = -1; % -E_jki
			A6(cnt, U_idx((j-1)*t+k)) = 1; % U_jk
			A6(cnt, W_idx((i-1)*t+k)) = 1; % W_ki
			cnt = cnt+1;
		end
	end
end
clear i j k pos cnt
b6 = ones(s*t*r, 1);
%A6 = sparse(A6);
A6 = int8tosparse(A6);

%
% sort transcripts by weight if there are no predefined transcripts;
% the transcripts will be sorted according to expression of the first 
% sample
% -w_t+w_t+1 <= 0
if 0%sort_novel_transcripts % TODO use this again and check if it makes a difference in the result
    A8 = zeros(t-length(predef_trans)-1, num_var);
	i = 1; % only for the first sample
	cnt = 1;
    for k=length(predef_trans)+1:t-1
        A8(cnt, W_idx((i-1)*t+k)) = -1; % -W_j
        A8(cnt, W_idx((i-1)*t+k+1)) = 1; % W_j+1
		cnt = cnt+1;
    end
    b8 = zeros(length(sort_idx)-1, 1);
	clear cnt k i
else
    A8 = zeros(1, num_var);
    b8 = 0;
end
A8 = sparse(A8);
	
	
%
% I_t
% indicator for transcripts with weight >0 in any sample
% sum_r W_tr - r*I_t <=0
% -I_1 <= -1 predict at least one transcript
A9 = zeros(t+1, num_var);
for j = 1:t
	for i=1:r
		A9(j, W_idx((i-1)*t+j)) = 1; % sum_i W_ji
	end
	A9(j, I_idx(j)) = -r; % I_j
end
clear j i
b9 = zeros(t+1, 1);
if isempty(predef_trans)
	% at least one transcript:
	A9(end, I_idx(1)) = -1; %I_1
	b9(end) = -1;
end
A9 = sparse(A9);
	
%
% S_ss
% spliced read penalty
%
% C_jkt = W_t*U_jt*prod_i=j+1^k-1{1-U_it}*U_kt
%
% (1) C_jkt<=U_jt
% (2) C_jkt<=U_kt
% (3) C_jkt<=1-U_it for all j<i<k
% (4) C_jkt<=W_t-U_jt+1-U_kt+1+ (sum_i U_it) 
% (5) C_jkt>= W_t+U_jt-1+U_kt-1- (sum_i U_it)
%
% compute expected intron coverage - observed intron coverage
% (6) D_jk = sum_t C_jkt - O_jk 
%
% for introns (j,k) not in splicegraph
% do not allow the usage of this intron:
% (7) U_{jt}+U_{kt} <= \sum_{i=j+1}^{k-1} U_{it} +1 

% (1) C_jkt-U_jt<=0					-> A10
% (2) C_jkt-U_kt<=0 				-> A11
% (3) C_jkt+U_it<=1 for all j<i<k	-> A12
A10 = zeros(r*t*c, num_var, 'int8');
b10 = zeros(r*t*c, 1);
A11 = zeros(r*t*c, num_var, 'int8');
b11 = zeros(r*t*c, 1);
if 0
	A12 = zeros(r*t*c*10, num_var, 'int8');
else
	nzmax = r*t*c*10;
	A12_idx = zeros(nzmax, 3);
end
cnt_rows = 0;
nz = 1;
for i=1:r
	intron_conf = compute_intron_list(admat(:,:,i));
	for x = 1:t
		for xx = 1:size(intron_conf, 1)
			j = intron_conf(xx, 1); % intron from segment j
			k = intron_conf(xx, 2); % intron to segment k
			cnt = (i-1)*t*c+(x-1)*c+xx ;
			A10(cnt, U_idx((j-1)*t+x)) = -1; % U_jx
			A10(cnt, C_idx(cnt)) = 1; % C_jkx
			A11(cnt, U_idx((k-1)*t+x)) = -1; % U_kx
			A11(cnt, C_idx(cnt)) = 1; % C_jkx
			for l=j+1:k-1
				cnt_rows = cnt_rows+1;
				if nz>size(A12_idx, 1)-100
					A12_idx = [A12_idx; zeros(size(A12_idx, 1), 3)];
				end
				A12_idx(nz, 1) = cnt_rows;						
				A12_idx(nz, 2) = U_idx((l-1)*t+x);
				A12_idx(nz, 3) = 1;
				nz = nz+1;
				A12_idx(nz, 1) = cnt_rows;						
				A12_idx(nz, 2) = C_idx(cnt);
				A12_idx(nz, 3) = 1;
				nz = nz+1;

				%if cnt_rows>size(A12, 1)
				%	A12 = [A12; zeros(r*t*c*5, num_var, 'int8')];
				%end
				%A12(cnt_rows, U_idx((l-1)*t+x)) = 1; % U_lx
				%A12(cnt_rows, C_idx(cnt)) = 1; % C_jkx
			end
		end
	end
end 
%A12(cnt_rows+1:end, :) = [];
A12_idx(nz:end, :) = [];
A12 = sparse(A12_idx(:, 1), A12_idx(:, 2), A12_idx(:, 3), cnt_rows, num_var);

b12 = ones(size(A12, 1), 1);
clear i x xx j k cnt l cnt_rows nz A12_idx

%A10 = sparse(A10);
%A11 = sparse(A11);
%A12 = sparse(A12);
A10 = int8tosparse(A10);
A11 = int8tosparse(A11);
%A12 = int8tosparse(A12);


% (4)  C_jkt-W_t+U_jt+U_kt-(sum_i U_it) <= 2
A13 = zeros(r*t*c, num_var, 'int8');
for i=1:r
	intron_conf = compute_intron_list(admat(:,:,i));
	for x = 1:t
		for xx = 1:size(intron_conf, 1)
			j = intron_conf(xx, 1);
			k = intron_conf(xx, 2);
			cnt = (i-1)*t*c+(x-1)*c+xx ;
			A13(cnt, U_idx((j-1)*t+x)) = 1; % U_jt
			A13(cnt, U_idx((k-1)*t+x)) = 1; % U_kt
			A13(cnt, C_idx(cnt)) = 1; % C_jkt
			A13(cnt, W_idx((i-1)*t+x)) = -1; % -W_t
			for l=j+1:k-1
				A13(cnt, U_idx((l-1)*t+x)) = -1; %-U_lt
			end
		end
	end 
end
b13 = 2*ones(r*t*c, 1);
%A13 = sparse(A13);
A13 = int8tosparse(A13);

clear i x xx j k cnt l
	
% (5) -C_jkt+W_t+U_jt+U_kt-(sum_i U_it) <= 2
A14 = zeros(r*t*c, num_var, 'int8');
for i=1:r
	intron_conf = compute_intron_list(admat(:,:,i));
	for x = 1:t
		for xx = 1:size(intron_conf, 1)
			j = intron_conf(xx, 1);
			k = intron_conf(xx, 2);
			cnt = (i-1)*t*c+(x-1)*c+xx ;
			A14(cnt, U_idx((j-1)*t+x)) = 1; % U_jx
			A14(cnt, U_idx((k-1)*t+x)) = 1; % U_kx
			A14(cnt, C_idx(cnt)) = -1; % -C_jkx
			A14(cnt, W_idx((i-1)*t+x)) = 1; % W_x
			for l=j+1:k-1
				A14(cnt, U_idx((l-1)*t+x)) = -1; %-U_lx
			end
		end
	end
end
b14 = 2*ones(r*t*c, 1);
%A14 = sparse(A14);
A14 = int8tosparse(A14);

clear i x xx j k cnt l
	
%
% (6) D_jk = sum_t C_jkt - O_jk 
% (6) sum_t C_jkt - D_jk =  O_jk
A15 = zeros(r*c, num_var);
b15 = zeros(r*c, 1);
for i=1:r
	intron_conf = compute_intron_list(admat(:,:,i));
	for xx = 1:size(intron_conf, 1)
		j = intron_conf(xx, 1);
		k = intron_conf(xx, 2);
		for x = 1:t
			cnt = (i-1)*t*c+(x-1)*c+xx ;
			A15((i-1)*c+xx, C_idx(cnt)) = 1; % C_jkx
		end
		A15((i-1)*c+xx, D_idx((i-1)*c+xx)) = -1; % -D_jk1
		A15((i-1)*c+xx, D_idx((i-1)*c+xx+c*r)) = 1; % D_jk2
	end 
	b15((i-1)*c+1:(i-1)*c+c) = intron_conf(:, 3); % observed coverage of intron 
	%if isempty(A15)
	%	A15 = zeros(1, num_var);
	%	b15 = 0;
	%end
end
clear i x xx j k cnt
	
n_of_equalities = n_of_equalities+size(A15, 1);
A15 = sparse(A15);
	
% for introns (j,k) not in splicegraph
% do not allow the usage of this intron:
% (7) U_{jt}+U_{kt} <= \sum_{i=j+1}^{k-1} U_{it} +1 
% (7) U_{jt}+U_{kt}-\sum_{i=j+1}^{k-1} U_{it} <= 1 

% negative formulation
% all not known introns are forbidden
if 0
	A18 = zeros(t*s*(s-1)/2, num_var);
else
	nzmax = t*s*(s-1)/2*10;
	A18_idx = zeros(nzmax, 3);
end
cnt = 0;
nz = 1;
for x = 1:t
	for j = 1:s
		% find the last segment that is connected to segment j
		% constraint matrix A3 will make sure that one of them is 
		% used if the segment is not terminal
		% => if the segment is not terminal introns larger than to 
		% the last connected segment are excluded any way and 
		% we do not need to do this here
		cnodes = find(admat(j,:, 1)>=-1);
		children = cnodes(cnodes>j);
		if terminal(j)
			maxs = s;
		elseif isempty(children)
			% this case is handled by A3
			% it makes sure that no downstream segment is used
			% thus there is also no intron
			continue;
		else
			maxs = max(children-1);
		end

		for k = j+1:maxs 
			if admat(j, k, 1)>=-1 % -2 is a invalid connection
				continue
			end
			cnt = cnt+1;

			if nz>size(A18_idx, 1)-(k-j)-10
				A18_idx = [A18_idx; zeros(size(A18_idx, 1), 3)];
			end

			%A18(cnt, nE+(j-1)*t+x) = 1; % U_jx
			A18_idx(nz, 1) = cnt;
			A18_idx(nz, 2) = U_idx((j-1)*t+x); 
			A18_idx(nz, 3) = 1; 
			nz = nz+1;
			%A18(cnt, nE+(k-1)*t+x) = 1; % U_kx	
			A18_idx(nz, 1) = cnt;
			A18_idx(nz, 2) = U_idx((k-1)*t+x);
			A18_idx(nz, 3) = 1;
			nz = nz+1;
			for l=j+1:k-1
				%A18(cnt, nE+(i-1)*t+x) = -1; %-U_ix
				A18_idx(nz, 1) = cnt;
				A18_idx(nz, 2) = U_idx((l-1)*t+x);
				A18_idx(nz, 3) = -1;
				nz = nz+1;
			end
		end
	end
	clear x j k l
end
%A18(cnt+1:end,:)=[];
A18_idx(nz:end, :) = [];
A18 = sparse(A18_idx(:, 1), A18_idx(:, 2), A18_idx(:, 3), cnt, num_var);
b18=ones(cnt, 1);
%assert(isequal(full(A18_sp), A18))
if isempty(A18_idx)
	A18 = zeros(1, num_var);
	b18 = 0;
end
A18 = sparse(A18);
	
	
%
% make sure that only transcripts with weight >=0 have segments
%
% \sum_j U_{jx}/S <= I_x \forall 1<x<T
%
% and that each used transcript has at least one segment
% I_x <= \sum_j U_{jx} \forall 1<x<T
% 
if 0 %isempty(predef_trans)
	A19 = zeros(2*t, num_var);
	b19 = zeros(2*t, 1);
	for i=1:r
		for x=length(predef_trans)+1:t
			A19(x, I_idx(x)) = 1; % -W_xi
			%A19(t+x, I_idx(x)) = 1; % I_x
			for j=1:s
				A19(x, U_idx((j-1)*t+x)) = 1/s; % U_jx
				%A19(t+x, U_idx((j-1)*t+x)) = -1; % U_jx
			end
		end
	end
	clear x j 
else
	A19=zeros(1, num_var);
	b19=0;
end
A19 = sparse(A19);


% pair stuff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if use_pair
	% P_{i,j,t} = U_{i,t}*U_{j,t}*I_t
	%
	% P_{i,j,t} \in {0,1}
	% P_{i,j,t} - 1/3*U_{i,t} - 1/3*U_{j,t}-1/3*I_t <=0
	% -P_{i,j,t} + U_{i,t} + U_{j,t} + I_t <= 2
	AP1 = zeros(2*np*t, num_var, 'int8');
	bp1 = [zeros(np*t, 1); 2*ones(np*t, 1)];
	obs = pair_observation;
	for x = 1:t
		for y=1:np
			j = obs(y, 1);
			k = obs(y, 2);
			idx = (x-1)*np+y;
			% P_{i,j,t} - 0.5*U_{i,t} - 0.5*U_{j,t} <=0
			AP1(idx, U_idx((j-1)*t+x)) = -1;	
			AP1(idx, U_idx((k-1)*t+x)) = -1;
			AP1(idx, I_idx(x)) = -1;
			AP1(idx, P_idx(idx)) = 3;
			% -P_{i,j,t} + U_{i,t} + U_{j,t} <= 1
			AP1(idx+t*np, U_idx((j-1)*t+x)) = +1;	
			AP1(idx+t*np, U_idx((k-1)*t+x)) = +1;	
			AP1(idx+t*np, I_idx(x)) = +1;
			AP1(idx+t*np, P_idx(idx)) = -1;
		end	
	end
	AP1 = int8tosparse(AP1);

	% make sure at least one transcripts can explain the pair
	% \sum_t P_{i,j,t} + SP_{i,j} >=1  
	% -\sum_t P_{i,j,t} - SP_{i,j} <= -1
	AP2 = zeros(np, num_var);
	bp2 = -ones(np, 1);
	for y=1:np
		j = obs(y, 1);
		k = obs(y, 2);
		AP2(y, SP_idx(y)) = -1;	
		for x = 1:t
			idx = (x-1)*np+y;
			AP2(y, P_idx(idx)) = -1;	
		end	
	end
	AP2 = sparse(AP2);
end	


% AP9: compute the length of new transcripts and make shure it is larger than len_cutoff
len_cutoff_idx=length(predef_trans)+1:t;

AP9 = zeros(length(len_cutoff_idx), num_var);
bp9 = zeros(length(len_cutoff_idx), 1);

if isfield(param, 'len_cutoff') && param.len_cutoff>1,
    fprintf('requiring new transcripts to have length >= %i, if they are used\n', param.len_cutoff);
    cnt = 0;
    for x = len_cutoff_idx % index in transcripts 
        cnt = cnt+1;
        for j = 1:s
            AP9(cnt, U_idx((j-1)*t+x)) = -len(j); % U_ix
        end	
        AP9(cnt, I_idx(x)) = param.len_cutoff ;
    end		
	clear x j
end
AP9 = sparse(AP9);


if use_cluster
	% sum over cluster assignment rows has to be equal to one
	%
	% (sum_{j=1}^nc M_{i,j}) = 1, forall i in {1,..,r}
	AC1 = zeros(r, num_var);
	bc1 = ones(r, 1);
	for i=1:r
		for j=1:nc
			AC1(i, M_idx((i-1)*nc+j)) = 1;
		end
	end
	AC1 = sparse(AC1);
	n_of_equalities = n_of_equalities+r;

	% tie weights to cluster centers if sample belongs to the cluster
	% w_{i,t} \leq m_{j,t} + 1 - M{i,j}
	% w_{i,t} \geq m_{j,t} - 1 + M{i,j}
	%
	% w_{i,t} - m_{j,t} + M{i,j} \leq 1
	% -w_{i,t} + m_{j,t} + M{i,j} \leq 1

	AC2 = zeros(2*r*nc*t, num_var);
	bc2 = ones(2*r*nc*t, 1);
	cnt = 0;
	for i=1:r
		for j=1:nc
			for x=1:t
				cnt = cnt+1;
				AC2(cnt, W_idx((i-1)*t+x)) 	= 1;
				AC2(cnt, m_idx((j-1)*t+x)) 	= -1;
				AC2(cnt, M_idx((i-1)*nc+j)) = 1;
				cnt = cnt+1;
				AC2(cnt, W_idx((i-1)*t+x)) 	= -1;
				AC2(cnt, m_idx((j-1)*t+x)) 	= 1;
				AC2(cnt, M_idx((i-1)*nc+j)) = 1;
			end
		end
	end
	AC2 = sparse(AC2);

	% penalty for number of clusters
	%
	% sum_{i=1}^r M_{1,j} \leq r*K_j, forall j in 1,...,nc
	AC3 = zeros(nc, num_var);
	bc3 = zeros(nc, 1);
	for j=1:nc
		for i=1:r
			AC3(j, M_idx((i-1)*nc+j)) = 1;
			AC3(j, K_idx(j)) = -r;
		end
	end
	AC3 = sparse(AC3);

	% sort clusters according to size
	AC4 = zeros(nc-1, num_var);
	bc4 = zeros(nc-1, 1);
	for j=1:nc-1
		for i=1:r
			AC4(j, M_idx((i-1)*nc+j)) = -1;
			AC4(j, M_idx((i-1)*nc+j+1)) = 1;
		end
	end
	AC4 = sparse(AC4);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% compute objective
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% regularizer for the number of transcripts (linear)
f = zeros(num_var, 1); 
if length(predef_trans)>0
	f(I_idx(1):I_idx(length(predef_trans))) = C.num_transcripts_predef;
	if t>length(predef_trans)
		f(I_idx(length(predef_trans)+1):I_idx(end)) = C.num_transcripts;
	end
else
	f(I_idx) = C.num_transcripts;
end
if use_cluster
	f(K_idx) = C.num_clusters;
end	

%Q = zeros(num_var);
Q = sparse([],[],[],num_var, num_var, num_var);
if ~use_LP_covloss,
	for i=1:r
    	for j=1:s
			readlen = 100;
			x1 = L_idx((i-1)*s+j);
			x2 = L_idx((i-1)*s+j+s*r);
			val = observed_cov(i, j); % this is the normalized coverage
			%val = val*cov_scale(i); % problem: variance is overestimated since this is the mean of many positions
			val = val*cov_scale(i)*len(j)/readlen; % rough estimate of the number of reads in that region

			if strcmp(param.loss, 'nb') || strcmp(param.loss, 'poisson')
				%[left, right] = get_coefficients(val);
				if strcmp(param.loss, 'poisson')
					[ll, lq, rl, rq] = get_coefficients_lq(val, 1, 0, 0);
				else
					[ll, lq, rl, rq] = get_coefficients_lq(val, param.eta1, param.eta2, param.lambda);
				end

        		Q(x1,x1) = 1/r*rq*C.exon_cov*(len(j)/readlen*cov_scale(i))^2;
        		Q(x2,x2) = 1/r*lq*C.exon_cov*(len(j)/readlen*cov_scale(i))^2;
        		f(x1) = 1/r*rl*C.exon_cov*len(j)/readlen*cov_scale(i);
        		f(x2) = 1/r*ll*C.exon_cov*len(j)/readlen*cov_scale(i);
			elseif strcmp(param.loss, 'l2')
				Q(x1,x1) = len(j)*C.exon_cov*beta;
        		Q(x2,x2) = len(j)*C.exon_cov*beta;
			end
		end
    end
	clear i j
end ;
	
for i=1:r
	intron_conf = compute_intron_list(admat(:,:,i));
	for j=1:c
		x1 = D_idx((i-1)*c+j);
		x2 = D_idx((i-1)*c+j+c*r);
		
		val = intron_conf(j, 3); % observed coverage of intron 
		val = val*cov_scale(i);

		if strcmp(param.loss, 'nb') || strcmp(param.loss, 'poisson')
			if strcmp(param.loss, 'poisson')
				[ll, lq, rl, rq] = get_coefficients_lq(val, 1, 0, 0);
			else
				[ll, lq, rl, rq] = get_coefficients_lq(val, param.eta1, param.eta2, param.lambda);
			end

		
			Q(x1,x1) = 1/r*rq*C.introns*(cov_scale(i))^2;
			Q(x2,x2) = 1/r*lq*C.introns*(cov_scale(i))^2; 
			f(x1) = 1/r*rl*C.introns*cov_scale(i);
			f(x2) = 1/r*ll*C.introns*cov_scale(i); 
		elseif strcmp(param.loss, 'l2')
			Q(x1,x1) = C.introns*cov_scale(i)/sum(cov_scale);
			Q(x2,x2) = C.introns*cov_scale(i)/sum(cov_scale); 
		end
	end
end
clear i j

if use_pair
	% deviation of the expected and observed pair coverage
	f(SP_idx) = C.pairs*pair_observation(:, 3);
end
	
lb = -1e20*ones(num_var, 1);
ub = 1e20*ones(num_var, 1);
lb(E_idx) = 0;
ub(E_idx) = 1;
lb(W_idx) = 0; 
ub(W_idx) = 1;
lb(C_idx) = 0;
ub(C_idx) = 1;
lb(L_idx) = 0;
lb(D_idx) = 0;
if use_pair
%	lb(X_idx) = 0;
	lb(SP_idx) = 0;
end
if use_cluster
	lb(m_idx) = 0;
	ub(m_idx) = 1;
end
	
binary_var_index = [U_idx, I_idx]; 
if use_pair
	binary_var_index = [binary_var_index P_idx]; 
end
if use_cluster
	binary_var_index = [binary_var_index M_idx K_idx];
end
	

% fix a set of transcripts 
assert(length(predef_trans)<=t)
for x = 1:length(predef_trans)
	for j = 1:s
		lb(U_idx((j-1)*t+x)) = predef_trans{x}(j);
		ub(U_idx((j-1)*t+x)) = predef_trans{x}(j);
	end
end


% set all but the first weight to zero
% these constraints will be removed one by one 
% during the solving; the last constraint will be removed first
%A20 = zeros(r*(t-1), num_var);
%b20 = zeros(r*(t-1), 1);
cnt = 0;
A20 = zeros(t-1, num_var);
b20 = zeros(t-1, 1);
for x = t:-1:length(predef_trans)+2
	cnt = cnt+1;
	for j = 1:s
		A20(cnt, U_idx((j-1)*t+x)) = 1;
	end
	%for i=1:r
	%	cnt = cnt+1;
	%	A20(cnt, W_idx((i-1)*t+x)) = 1; % -W_t
	%end 
end
clear cnt


% concat constraints
% A1 (if not linear) and A15 and AP8 are equality constraints 
% thats why they have to be on top of the stack
% A20 hast to be in the very end of the A matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~use_pair && ~use_cluster
    A = [A15; A1; A2; A3; A4; A5; A6; A8; A9; A10; A11; A12; A13; A14; A18; A19; AP9; A20];
    b = [b15; b1; b2; b3; b4; b5; b6; b8; b9; b10; b11; b12; b13; b14; b18; b19; bp9; b20];
elseif ~use_pair && use_cluster
    A = [AC1; A15; A1; A2; A3; A4; A5; A6; A8; A9; A10; A11; A12; A13; A14; A18; A19; AP9; AC2; AC3; AC4; A20];
    b = [bc1; b15; b1; b2; b3; b4; b5; b6; b8; b9; b10; b11; b12; b13; b14; b18; b19; bp9; bc2; bc3; bc4; b20];
elseif use_pair && use_cluster
    A = [AC1; A15; A1; A2; A3; A4; A5; A6; A8; A9; A10; A11; A12; A13; A14; A18; A19; AP1; AP2; AC2; AC3; AC4; A20];
    b = [bc1; b15; b1; b2; b3; b4; b5; b6; b8; b9; b10; b11; b12; b13; b14; b18; b19; bp1; bp2; bc2; bc3; bc4; b20];
elseif use_pair && ~use_cluster
    A = [A15; A1; A2; A3; A4; A5; A6; A8; A9; A10; A11; A12; A13; A14; A18; A19; AP1; AP2; A20];
    b = [b15; b1; b2; b3; b4; b5; b6; b8; b9; b10; b11; b12; b13; b14; b18; b19; bp1; bp2; b20];
end ;

% solve min_x x'Qx+f'x
% s.t. 	Ax<=b
% 	 	lb<=x<=ub
solvetime = cputime;
[result, how] = solve_mip(Q, f, A, b, lb, ub, binary_var_index, n_of_equalities, param.time_limit, r, t, s, U_idx, length(predef_trans));
solvetime = cputime-solvetime;

b(length(b)-length(b20)+1:end) = 1;

%colors= 'rg';
%for i = 1:num_is
%	offset = sum(num_i_col(1:i-1));
%	res_E = result(offset+E_idx)';
%	res_U = result(offset+U_idx)';
%	res_W = result(offset+W_idx)';
%	res_L = result(offset+L_idx)'; 
%	res_I = result(offset+I_idx)'; 
%	plot(1:length(res_E), res_E, colors(i)), hold on
%end

for i = 1:r
	res_U = result(U_idx)';
	res_I = result(I_idx)'; 
	res_E = result(E_idx((i-1)*s*t+1:(i-1)*s*t+s*t))';
	res_W = result(W_idx((i-1)*t+1:(i-1)*t+t))';
	res_L = result(L_idx((i-1)*s+1:(i-1)*s+s))'; 
	res_C = result(C_idx((i-1)*c*t+1:(i-1)*c*t+c*t))';
	res_D = result(D_idx((i-1)*c+1:(i-1)*c+c))';
	
	%if use_pair
	%	res_X  = result(X_idx((i-1)*p*t+1:(i-1)*p*t+p*t))';
	%	res_PL = result(PL_idx)';
	%	res_PU = result(PU_idx)';
	%	res_DP = result(DP_idx((i-1)*p+1:(i-1)*p+p))';
	%end
	%pair_objective=sum(res_DP.^2)*C.pairs

	admat = admat_all(:,:,i);
	intron_conf = compute_intron_list(admat);
	if num_var<5000
		print_result;
	end
end

transcripts = {};
weights = [];
res_U = result(U_idx)';
for j = 1:t
	transcripts{j} = [];
	for i = 1:r
		weights(i, j) = result(W_idx((i-1)*t+j));
	end
	for k = 1:s
		%if res_I(j)<0.5
		%	continue
		%end
		if res_U((k-1)*t+j)>0.5
			transcripts{j} = [transcripts{j} k];
		end
	end
end

clusters = zeros(1, r);
if use_cluster
	cluster_mat = zeros(r, nc);
	for i=1:r
		for j=1:nc
			cluster_mat(i, j) = result(M_idx((i-1)*nc+j));
		end
	end
	idx = find(sum(cluster_mat, 1)==0);
	cluster_mat(:,idx) = [];
	cluster_mat
	for j = 1:r
		clusters(j) = find(cluster_mat(j, :)>0.5);
	end
end

if ~isempty(strfind(how, 'No integer feasible solution exists'))
	%keyboard
	how
end

sA20 = size(A20, 1);
A = A(1:size(A, 1)-sA20, :);
b = b(1:size(A, 1));
if ~(all(A*result-b<=1e-3)),
    how=[how '~all(A*result-b<=1e-3)'] 
	%keyboard
end ;
if ~all(result>=lb-1e-3 & result<=ub+1e-3)
    how=[how '~all(result>=lb-1e-3 & result<=ub+1e-3)'] 
	%keyboard
end 
if ~all(abs(res_U)<1e-1|abs(res_U-1)<1e-1)
    how=[how '~all(abs(res_U)<1e-1|abs(res_U-1)<1e-1)'] 
	%keyboard
end ;
if ~all(res_W>-1e-7)
    how=[how '~all(res_W>-1e-7)'] 
	%keyboard
end ;

if exist('fn_res', 'var') && ~isempty(fn_res)
	save(fn_res, 'transcripts', 'weights', 'how', 'solvetime', 'clusters', 'PAR');
end
return


% bugfix stuff

idx = find(A*result-b>0.1)
find(A(idx, :)~=0)

for j = 1:19
	if ~exist(sprintf('A%i', j), 'var')
		continue; 
	end
	eval(sprintf('res = A%i*result -b%i;', j, j)); 
	if any(res>0.1), 
		j
	end
end

for j = 1:length(initial)
	fprintf('%i ', initial(j))
end
fprintf('\n')
for j = 1:length(initial)
	fprintf('%i ', terminal(j))
end
fprintf('\n')
for k = 1:length(predef_trans)
	for j = 1:length(initial)
		fprintf('%i ', predef_trans{k}(j))
	end
	fprintf('\n')
end
function [ll, lq, rl, rq] = get_coefficients_lq(obs, eta1, eta2, lambda)

	%load('/fml/ag-raetsch/nobackup/projects/mip/human_sim/nb_fit/var0.02.mat', 'xpos', 'left_l', 'left_q', 'right_l', 'right_q')
	%load('/fml/ag-raetsch/nobackup/projects/mip/human_sim/nb_fit/var0.02_fit0.1x-3x.mat', 'xpos', 'left_l', 'left_q', 'right_l', 'right_q') one of the files was accidently overwritten

	fname = create_loss_parameters(eta1, eta2, lambda, '~/tmp');
	load(fname, 'xpos', 'left_l', 'left_q', 'right_l', 'right_q')

	% find bin
	for j=1:length(xpos)
		if obs<=xpos(j)
			break;
		end
	end
	if j == 1
		ll = left_l(1);
		lq = left_q(1);
		rl = right_l(1);
		rq = right_q(1);
		return
	end
	if j==length(xpos)
		ll = left_l(end);
		lq = left_q(end);
		rl = right_l(end);
		rq = right_q(end);
		return
	end

	ll = ((obs-xpos(j-1))*left_l(j)+(xpos(j)-obs)*left_l(j-1))/(xpos(j)-xpos(j-1));
	lq = ((obs-xpos(j-1))*left_q(j)+(xpos(j)-obs)*left_q(j-1))/(xpos(j)-xpos(j-1));
	rl = ((obs-xpos(j-1))*right_l(j)+(xpos(j)-obs)*right_l(j-1))/(xpos(j)-xpos(j-1));
	rq = ((obs-xpos(j-1))*right_q(j)+(xpos(j)-obs)*right_q(j-1))/(xpos(j)-xpos(j-1));
return

function [left, right] = get_coefficients(x)

	%load('~raetsch/var/coeff_nb_leftright_v=0.2.mat', 'w1_x', 'w1_y', 'w2_x', 'w2_y')
	load('~raetsch/var/coeff_nb_leftright_v=0.02.mat', 'w1_x', 'w1_y', 'w2_x', 'w2_y')

	% find bin
	for j=1:length(w1_x)
		if x<w1_x(j)
			break;
		end
	end
	if j == 1
		left = w1_y(1);
		right = w2_y(1);
		return
	end
	if j==length(w1_x)
		left = w1_y(end);
		right = w2_y(end);
		return
	end

	left = ((x-w1_x(j-1))*w1_y(j)+(w1_x(j)-x)*w1_y(j-1))/(w1_x(j)-w1_x(j-1));
	right = ((x-w2_x(j-1))*w2_y(j)+(w2_x(j)-x)*w2_y(j-1))/(w2_x(j)-w2_x(j-1));
return

function intron_conf = compute_intron_list(admat)

	s = size(admat, 1);
	intron_conf = zeros(s*(s-1)/2, 3);
	cnt = 0;
	for j = 1:s
		for k = j+1:s
			if admat(j, k)>=0
				% valid intron not just a neighboring segment
				cnt = cnt+1;
				intron_conf(cnt,:) = [j, k, admat(j, k)];
			end	
		end
	end
	intron_conf(cnt+1:end, : ) = [];
return
