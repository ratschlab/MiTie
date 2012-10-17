
run=0 ;
time=[] ;
time_=[] ;
diffs=[] ;
%while(1)

lpenv=cplex_license(0) ;

num_var=10;ceil(rand(1)*400) ;
num_con=10;ceil(rand(1)*1600) ;

c=-rand(num_var, 1) ;
A=rand(num_con, num_var) ;
b=rand(num_con, 1) ;
H=rand(num_var) ; H=H*H' ;

l=zeros(num_var,1) ;
u=10*ones(num_var,1) ;

neq=0;%ceil(num_con*rand(1)) ;


tic;[x,y,how]=qp_solve(lpenv,sparse(H),c,sparse(A),b,l,u,neq,1); t=toc;
%[x,y,how,p_lp]=lp_solve(lpenv,c,sparse(A),b,l,u,neq,0,'primal'); 

cplex_close(lpenv) ; lpenv=0 ;

clear functions
assert(neq==0) ;
tic;[x_,FVAL,EXITFLAG,OUTPUT,y_]=quadprog(sparse(H), c, sparse(A), b, [], [], l, u); t_=toc;

run=run+1 ;
diffs(run)=norm(x-x_) ;

md=mean(diffs)







