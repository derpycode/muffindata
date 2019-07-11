%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% f2v
% Desciption: convert tracer field to vector for use in matrix
% calculations, based on code from Jake Gebbie.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [vector] = f2v(field,i,j,k)

  Nfield = length(i);
  vector = zeros(Nfield,1);
  for ni = 1:Nfield
    vector(ni) = field(k(ni),j(ni),i(ni));
  end