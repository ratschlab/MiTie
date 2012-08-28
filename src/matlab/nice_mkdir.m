
function [ ok ] = nice_mkdir( dir );
if ~isdir(dir)
  ok = ( 0 == unix( [ 'mkdir -p ', dir ] ) );
  ok = ok & ( 0 == unix( [ 'chmod g+t ', dir ] ) );
end

