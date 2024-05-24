function [F] = assembleforc(iele,nnpe,conn,fe,F)

%FEM_Assemble     assembles the elemental force and stiffness matrices into
%                 the global counterparts

for  i = 1:1:nnpe
    
    ii = conn(iele,i);
    
    F(ii,1) = F(ii,1) + fe(i,1);
    
end

