%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% v2f
% Desciption: convert tracer vector to field for use in matrix
% calculations, based on code from Jake Gebbie.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [field] = v2f(vector,i,j,k)

NZ = max(k);
NY = max(j);
NX = max(i);
NVEC = length(vector);
% JDW 19/10/15: temporary hack to make sure full grid is represented
% n.b. this only accounts for 36x36x8/16 and 18x18x16/8 grids crudely!
if NZ<16 & NZ>8
    NZ=16;
elseif NZ<8
    NZ=8;
end

if NY<36 & NZ>18
    NY=36;
elseif NY<18
    NY=18;
end

if NX<36 & NX>18
    NX=36;
elseif NX<18
    NX=18;
end

field = nan(NZ,NY,NX);
for nv = 1:NVEC
 field(k(nv),j(nv),i(nv)) = vector(nv);
end
