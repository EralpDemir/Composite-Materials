# Composite-Materials
Contains Matlab codes for design optimization of variable stiffness composites using lamination parameters within a finite element method framework.

"MAIN.m" is the main runner file.

Them-files with the prefixes contain the functions related to:
- "FEM_": finite elements
- "MAT_": material stiffness
- "LP_" : lamination parameter
- "LSC_": least-squares and curvature constraints
- "OPT_": optimization

The userinputs are given in the m-file the suffix "_inputs". 
i.e. "FEM_inputs.m" contains the inputs related with FEM calculations.

- "/compliance_optimization": in-plane and out-of-plane optimization for compliance/stiffness

    https://doi.org/10.1016/j.compositesb.2019.02.004

    https://doi.org/10.1016/j.compstruct.2019.111280

    https://doi.org/10.1016/j.addma.2020.101728


- "/frequency_optimization": plate natural frequency optimization

    https://doi.org/10.1016/j.compstruct.2021.115151


- "/buckling_optimization": critical buckling load optimization

    submitted


- "/stream_function_analysis": Stream function fit for continuous tow paths
