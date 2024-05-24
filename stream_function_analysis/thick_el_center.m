function  [thick_el] = thick_el_center(numel, nnpe, conn, thick, N0)
%LSC_angel0 calculate the angles at the element centers

thick_el= zeros(numel,1);


for iele = 1:1:numel

    thick_vals_el = zeros(nnpe,1);

    for inod = 1:1:nnpe

        NodeNo = conn(iele,inod);

        thick_vals_el(inod,1) = thick(NodeNo);

    end

    thick_el(iele) = N0*thick_vals_el;

end

    
