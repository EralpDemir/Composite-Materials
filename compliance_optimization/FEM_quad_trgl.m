function [xi, eta, w] = FEM_quad_trgl(m)

%------------------------------------------------
% Abscissas (xi, eta) and weights (w)
% for Gaussian integration over a flat triangle
% in the xi-eta plane
%
% Integration is performed with respect
% to the triangle barycentric coordinates
%
% SYMBOLS:
% -------
%
% m: order of the quadrature
%    choose from 1,3,4,5,7,9,12,13
%    Default value is 7
%------------------------------------------------

if ( (m ~= 1) && (m ~=  3) && (m ~=  4) && (m ~= 6) && (m ~= 7) ...
        && (m ~= 9) && (m ~= 12) && (m ~= 13) )
    
    disp('');
    disp(' Gauss_trgl: Chosen number of points');
    disp('   is not available; Will take m=7');
    m = 7;
    
end

if m == 1
    
    xi(1) = 1.0/3.0; eta(1) = 1.0/3.0; w(1) = 1.0;
    
elseif m == 3
    
    xi(1) = 1.0/6.0; eta(1) = 1.0/6.0; w(1) = 1.0/3.0;
    xi(2) = 2.0/3.0; eta(2) = 1.0/6.0; w(2) = w(1);
    xi(3) = 1.0/6.0; eta(3) = 2.0/3.0; w(3) = w(1);
    
elseif m == 4
    
    xi(1) = 1.0/3.0; eta(1) = 1.0/3.0; w(1) = -27.0/48.0;
    xi(2) = 1.0/5.0; eta(2) = 1.0/5.0; w(2) =  25.0/48.0;
    xi(3) = 3.0/5.0; eta(3) = 1.0/5.0; w(3) =  25.0/48.0;
    xi(4) = 1.0/5.0; eta(4) = 3.0/5.0; w(4) =  25.0/48.0;
    
elseif m == 6
    
    al = 0.816847572980459;
    be = 0.445948490915965;
    ga = 0.108103018168070;
    de = 0.091576213509771;
    o1 = 0.109951743655322;
    o2 = 0.223381589678011;
    
    xi(1) = de; xi(2) = al; xi(3) = de; xi(4) = be;
    xi(5) = ga; xi(6) = be;
    
    eta(1) = de;eta(2) = de;eta(3) = al;eta(4) = be;
    eta(5) = be;eta(6) = ga;
    
    w(1) = o1;w(2) = o1;w(3) = o1;w(4) = o2;
    w(5) = o2;w(6) = o2;
    
elseif m == 7
    
    al = 0.797426958353087;
    be = 0.470142064105115;
    ga = 0.059715871789770;
    de = 0.101286507323456;
    o1 = 0.125939180544827;
    o2 = 0.132394152788506;
    
    xi(1) = de;
    xi(2) = al;
    xi(3) = de;
    xi(4) = be;
    xi(5) = ga;
    xi(6) = be;
    xi(7) = 1.0/3.0;
    
    eta(1) = de;
    eta(2) = de;
    eta(3) = al;
    eta(4) = be;
    eta(5) = be;
    eta(6) = ga;
    eta(7) = 1.0/3.0;
    
    w(1) = o1;
    w(2) = o1;
    w(3) = o1;
    w(4) = o2;
    w(5) = o2;
    w(6) = o2;
    w(7) = 0.225;
    
elseif m == 9
    
    al = 0.124949503233232;
    qa = 0.165409927389841;
    rh = 0.797112651860071;
    de = 0.437525248383384;
    ru = 0.037477420750088;
    o1 = 0.205950504760887;
    o2 = 0.063691414286223;
    
    xi(1) = de;
    xi(2) = al;
    xi(3) = de;
    xi(4) = qa;
    xi(5) = ru;
    xi(6) = rh;
    xi(7) = qa;
    xi(8) = ru;
    xi(9) = rh;
    
    eta(1) = de;
    eta(2) = de;
    eta(3) = al;
    eta(4) = ru;
    eta(5) = qa;
    eta(6) = qa;
    eta(7) = rh;
    eta(8) = rh;
    eta(9) = ru;
    
    w(1) = o1;
    w(2) = o1;
    w(3) = o1;
    w(4) = o2;
    w(5) = o2;
    w(6) = o2;
    w(7) = o2;
    w(8) = o2;
    w(9) = o2;
    
elseif m == 12
    
    al = 0.873821971016996;
    be = 0.249286745170910;
    ga = 0.501426509658179;
    de = 0.063089014491502;
    rh = 0.636502499121399;
    qa = 0.310352451033785;
    ru = 0.053145049844816;
    o1 = 0.050844906370207;
    o2 = 0.116786275726379;
    o3 = 0.082851075618374;
    
    xi(1)  = de;
    xi(2)  = al;
    xi(3)  = de;
    xi(4)  = be;
    xi(5)  = ga;
    xi(6)  = be;
    xi(7)  = qa;
    xi(8)  = ru;
    xi(9)  = rh;
    xi(10) = qa;
    xi(11) = ru;
    xi(12) = rh;
    
    eta(1)  = de;
    eta(2)  = de;
    eta(3)  = al;
    eta(4)  = be;
    eta(5)  = be;
    eta(6)  = ga;
    eta(7)  = ru;
    eta(8)  = qa;
    eta(9)  = qa;
    eta(10) = rh;
    eta(11) = rh;
    eta(12) = ru;
    
    w(1)  = o1;
    w(2)  = o1;
    w(3)  = o1;
    w(4)  = o2;
    w(5)  = o2;
    w(6)  = o2;
    w(7)  = o3;
    w(8)  = o3;
    w(9)  = o3;
    w(10) = o3;
    w(11) = o3;
    w(12) = o3;
    
elseif m == 13
    
    al = 0.479308067841923;
    be = 0.065130102902216;
    ga = 0.869739794195568;
    de = 0.260345966079038;
    rh = 0.638444188569809;
    qa = 0.312865496004875;
    ru = 0.048690315425316;
    o1 = 0.175615257433204;
    o2 = 0.053347235608839;
    o3 = 0.077113760890257;
    o4 =-0.149570044467670;
    
    xi(1)  = de;
    xi(2)  = al;
    xi(3)  = de;
    xi(4)  = be;
    xi(5)  = ga;
    xi(6)  = be;
    xi(7)  = qa;
    xi(8)  = ru;
    xi(9)  = rh;
    xi(10) = qa;
    xi(11) = ru;
    xi(12) = rh;
    xi(13) = 1.0/3.0;
    
    eta(1)  = de;
    eta(2)  = de;
    eta(3)  = al;
    eta(4)  = be;
    eta(5)  = be;
    eta(6)  = ga;
    eta(7)  = ru;
    eta(8)  = qa;
    eta(9)  = qa;
    eta(10) = rh;
    eta(11) = rh;
    eta(12) = ru;
    eta(13) = 1.0/3.0;
    
    eta(1)  = de;
    eta(2)  = de;
    eta(3)  = al;
    eta(4)  = be;
    eta(5)  = be;
    eta(6)  = ga;
    eta(7)  = ru;
    eta(8)  = qa;
    eta(9)  = qa;
    eta(10) = rh;
    eta(11) = rh;
    eta(12) = ru;
    eta(13) = 1.0/3.0;
    
    w(1)  = o1;
    w(2)  = o1;
    w(3)  = o1;
    w(4)  = o2;
    w(5)  = o2;
    w(6)  = o2;
    w(7)  = o3;
    w(8)  = o3;
    w(9)  = o3;
    w(10) = o3;
    w(11) = o3;
    w(12) = o3;
    w(13) = o4;
    
end

return
