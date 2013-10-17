function lpenv = cplex_license(waitflag, num_cpus)
%  lpenv = cplex_license(waitflag, num_cpus)

global lpenv ;

if nargin<1, waitflag = 0 ; end ;
if nargin<2, num_cpus = 1 ; end ;

lpenv = 0 ;
while (lpenv==0)
  nc = num_cpus ;

  while ~lpenv && nc<4,
    
    lpenv = cplex_init_quit(0) ;
    
    if ~lpenv, 
      nc=nc+1 ; 
    end ;
  end ;

  if lpenv == 0 & ~waitflag,
    break ;
  end ;
  if lpenv==0,
    disp('\nwaiting for cplex license') ;
    pause(60)
  end;

end ;

