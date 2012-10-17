function lpenv = cplex_license(waitflag, num_cpus)
%  lpenv = cplex_license(waitflag, num_cpus)

global lpenv ;

if nargin<1, waitflag = 0 ; end ;
if nargin<2, num_cpus = 1 ; end ;

[engine, environment] = determine_engine() ;
if isequal(environment, 'internal')
  envstr1 = 'ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-academic.ilm';
  %envstr1 = 'ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-820980.ilm' ;
  envstr2 = 'ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-703040.ilm' ;
  envstr3 = 'ILOG_LICENSE_FILE=/fml/ag-raetsch/share/software/ilog/licenses/access-587380.ilm' ;
else
  envstr1 = 'ILOG_LICENSE_FILE=/home/galaxy/svn/tools/cplex/cplex91/access-820980.ilm' ;
  envstr2 = 'ILOG_LICENSE_FILE=/home/galaxy/svn/tools/cplex/cplex91/access-703040.ilm' ;
  envstr3 = 'ILOG_LICENSE_FILE=/home/galaxy/svn/tools/cplex/cplex91/access-587380.ilm' ;
end ;

lpenv = 0 ;
while (lpenv==0)
  nc = num_cpus ;

  while ~lpenv && nc<4,
    
    if nc == 0,
      envstr = getenv('ILOG_LICENSE_FILE') ;
    elseif nc == 1,
      envstr = envstr1 ;
    elseif nc == 2,
      envstr = envstr2 ;
    elseif nc > 2,
      envstr = envstr3 ;
    end ;
    
    fprintf('\ntrying license file %s\n', strrep(envstr, 'ILOG_LICENSE_FILE=', '')) ;
    lpenv = cplex_init_quit(0, envstr) ;
    
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

