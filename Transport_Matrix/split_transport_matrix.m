function [A1,B,Ii]=split_transport_matrix(A,Ib)
% function to break up A in dC/dt = A*C so that 
% dC/dt = A1*C + B*Cb with C containing interior points only.
% Ib is a vector of indices of BOUNDARY points.
% from Khatiwala code

if ~iscell(A) % A is assumed 2-D
  nb=length(Ib);
  ni=size(A,1)-nb;

% find indices of 'interior' points
  Ii=find(~ismember([1:(nb+ni)]',Ib));
%   Ii=[1:ni+nb]';
%   for j=1:nb
%       Ii=Ii(Ii~=Ib(j));
%   end
  A1=A(Ii,Ii);
  if nargout>1
    B=A(Ii,Ib);
  end
else
  nb=length(Ib);
  nm=length(A);
  ni=size(A{1},1)-nb;

% find indices of 'interior' points
  Ii=find(~ismember([1:(nb+ni)]',Ib));
%   Ii=[1:ni+nb]';
%   for j=1:nb
%       Ii=Ii(Ii~=Ib(j));
%   end
  for im=1:nm, 
     A1{im}=sparse(length(Ii),length(Ii));
     if nargout>1
       B{im}=sparse(length(Ii),length(Ib)); 
     end
  end
  for im=1:nm,
     A1{im}=A{im}(Ii,Ii);
     if nargout>1
       B{im}=A{im}(Ii,Ib); 
     end
  end     
end

