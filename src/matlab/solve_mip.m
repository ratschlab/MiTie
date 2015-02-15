function [x2, how] = solve_mip(Q, c, A, b, lb, ub, binary_idx, n_of_equalities, time_limit, r, t, s, U_idx, num_predef)

%addpath ~/svn/tools/cplex/cplex101
%addpath ~/svn/tools/cplex/cplex91
%addpath /cbio/grlab/share/software/ilog/cplex91
%addpath /cbio/grlab/share/git/tools/cplex
%addpath /cbio/grlab/home/jonas/mpghome/svn/projects/MiTie/src/matlab/cplex_gun %use correct version of cplex_license
%addpath /cbio/grlab/home/jonas/mpghome/svn/projects/MiTie/src/matlab %use correct version of cplex_license
addpath matlab/cplex/

disp('calling solve_mip')

global lpenv ;
if ~issparse(A)
	A = sparse(A);
end

lpenv

% get license 
if isempty(lpenv)
    lpenv=cplex_license(1,1)
end ;

if lpenv==0,
    pause(20) ;
    rproc_rerun() ;
    error('no cplex license') ;
end ;

if nnz(Q)==0,
    p_lp=lp_gen(lpenv, c, A, b, lb, ub, n_of_equalities)
else
    p_lp=qp_gen(lpenv, sparse(Q), c, A, b, lb, ub, n_of_equalities)
end ;

[x1,lambda1,how]=lp_resolve(lpenv, p_lp, 1, 'primal');
%tic;[x1,lambda1,how]=lp_resolve(lpenv, p_lp, 1, 'dual');toc;
how

if strcmp(how, 'CPX_STAT_ABORT_OBJ_LIM')
	keyboard
end

x1'*Q*x1+x1'*c
how=lp_set_param(lpenv,'CPX_PARAM_EPGAP', 0.01) ; % for dream6 0.01
how=lp_set_param(lpenv,'CPX_PARAM_EPINT', 0.001) ;% default 1e-5: branch only if varable deviates more than that from an integer 
%how=lp_set_param(lpenv,'CPX_PARAM_HEURFREQ',200) ; how 
%how=lp_set_param(lpenv,'CPX_PARAM_HEURISTIC', -1) ;

how=lp_set_param(lpenv,'CPX_PARAM_PROBE', 3) ; % !!
how=lp_set_param(lpenv,'CPX_PARAM_MIPTHREADS', 8)  ;

how=lp_set_param(lpenv,'CPX_PARAM_TILIM', time_limit) ; 
%how=lp_set_param(lpenv,'CPX_PARAM_NODESEL', 3) ;

%how=lp_set_param(lpenv,'CPX_PARAM_SUBALG', 2) ;
%how=lp_set_param(lpenv,'CPX_PARAM_SUBMIPNODELIM', 10) ;

how=lp_set_param(lpenv,'CPX_PARAM_SYMMETRY', 3) ;
%how=lp_set_param(lpenv,'CPX_PARAM_COVERS', 3) ;

if nnz(Q)==0,
    lp_chgprobtype(lpenv, p_lp, 1)  % CPXPROB_MILP
else
    lp_chgprobtype(lpenv, p_lp, 7)  % CPXPROB_MIQP
end ;

for i=binary_idx,
  lp_chgctype(lpenv, p_lp, i-1,double('B'));
end ;

x2 = [];
if t<2
	[x2,lambda2,how]=lp_resolve(lpenv, p_lp, 3, 'mip');
	how
elseif 0
	% remove constraints W<0 for more and more transcripts
	cnt=0;
	num_rows = size(A,1)
	while cnt<t-1
		[x,lambda2,how2]=lp_resolve(lpenv, p_lp, 3, 'mip');
		how2
		cnt = cnt+1;

		if ~isempty(strfind(how2, 'Time limit exceeded'))
			relgap = getmiprelgap(lpenv, p_lp);
			if isempty(x2) || relgap<0.1
				x2 = x;
				how = how2;
			end
			how = sprintf('%s;time_limit:%i;gap:%.2f', how, cnt, relgap*100)
			break
		end
		from = num_rows-r;
		to = num_rows-1;
		num_rows = num_rows-r;
		h = lp_delrows(lpenv,p_lp,from,to,10);

		x2 = x;
		how = how2;
	end
else
	cnt = 0;
	num_added = 0;
	num_rows = size(A,1)-1; % zero based
	num_var = size(A, 2);

	while all(A(num_rows+1, :)==0)
		% for technical reasons there are as many empty rows in A 
		% as there are annotated transcripts
		disp_ = 0;
		h = lp_delrows(lpenv,p_lp,num_rows,num_rows,disp_);
		num_rows = num_rows-1;
		cnt = cnt+1;
	end

	while 1 % cnt<t
		[x,lambda2,how2]=lp_resolve(lpenv, p_lp, 3, 'mip');
		how2
		cnt = cnt+1;

		if ~isempty(strfind(how2, 'Time limit exceeded'))
			if exist('getmiprelgap', 'file')==3
				relgap = getmiprelgap(lpenv, p_lp);
			else
				relgap = 0.99;
			end
			if isempty(x2) || relgap<0.1
				x2 = x;
				how = how2;
			end
			how = sprintf('%s;time_limit:%i;gap:%.2f', how, cnt, relgap*100)
			break
		end
		if cnt>=t
			x2 = x;
			how = how2;
			break;
		end
		disp_ = 0;
		h = lp_delrows(lpenv,p_lp,num_rows,num_rows,disp_);
		num_rows = num_rows-1;

		% add constraints to fix the previous solution
		A_new = zeros(s, num_var);
		b_new = zeros(s, 1);
		for j = 1:s
			A_new(j, U_idx((j-1)*t+cnt)) = 1;
			b_new(j) = round(x(U_idx((j-1)*t+cnt))); % fix U_{j,cnt} to the result from the previous run
		end
		neq = s;
		A_new = sparse(A_new');
		disp_ = 10;
		how3 = lp_addrows(lpenv,p_lp,A_new,b_new,neq,disp_)
		num_added = num_added+length(b_new);

		x2 = x;
		how = how2;
	end
	if 0%optimal
		h = lp_delrows(lpenv,p_lp,num_rows+1,num_rows+num_added,disp_);
		[x,lambda2,how2]=lp_resolve(lpenv, p_lp, 3, 'mip');
		if ~isempty(strfind(how2, 'Time limit exceeded'))
			x2 = x;
			how = how2;
		end
	end
end
%x2'*Q*x2+x2'*c

%if ~isequal(how, 'OK')
%    keyboard;
%end ;
