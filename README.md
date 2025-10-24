This project supports the manuscript with title 
A Corrosion Correction for Minimal Weight Design of Truss Structures under Compliance Constraint
submitted to Technische Mechanik by Anton and Mykola Tkachuk in 2025. The version is finished on 2025-10-24.

The files are split in 4 groups
  1. Solvers
  2. Preparation
  3. Example definitions
  4. Visualization and checking

I. Solvers

  Aim: compact expression of the optimization problem via CVX

  Two versions for fixed width problem (FW with and w/o self-weight)

    S1:: MinVolumeGivenComplRust.m
    S2:: MinVolumeGivenComplRustSelfWeight.m
    (Used in example 1 and 2)

  Two versions for the fixed aspect ratio case (FAR with and w/o proposed scaling)

    S3:: MinVolumeGivenComplAndFreq_fixed_aspect_ratio.m
    S4:: MinVolumeGivenComplAndFreq_fixed_aspect_ratio_scaling_off.m
    (Used in example 1 and 2)

  A version for fixed width problem with elliptical uncertanty (FWUL)
  
    S5:: MinVolumeGivenComplRustMultLoad.m
    (Used in example 3 only)

II. Preparation

  Aim: returns all the quantities to express the total mass, stiffness matrix, mass matrix in a compact way
  K=Qfull*diag(A)*Qfull'
  M=diag(Jfull*A) 
  total_weight= rho*A'*Lvector 

    P1:: SDTruss2D_AssembleMatrRust.m
    P2:: SDTruss2D_AssembleMatrRustMultLoad.m
    P3:: SDTruss2D_AssembleMatrRust_fixed_aspect_ratio.m

III. Example definitions

Example 1
  E1.1:: volume_min_given_compl_2D_17_truss_Ascaled_SingleLoad.m
  (case 0,1,2 depending on values defined == FW, used in figure 6)
  E1.2:: volume_min_given_compl_2D_17_truss_fixed_aspect_ratio.m
  (case 3,4 == FAR, used in figure 6 and 7)
  E1.3:: volume_min_given_compl_2D_17_truss_Ascaled_SW.m
  (case 5,6 == FW + self-weight, used in figure 9)

Example 2

  E2.1:: volume_min_given_compl_2D_33_truss_Ascaled3.m
  (case 1, used in figure 10 middle and right)
  E2.2:: volume_min_given_compl_2D_33_truss_fixed_aspect_ratio.m
  (case 2,3,4)

Example 3

  E3:: volume_min_given_compl_2D_30x30_truss_Cell60cm_MultLoad.m
  (used in figure 12)
  
  another file is given just to show computational time for smaller problem sizes
  volume_min_given_compl_2D_20x20_truss_GridStep60_MultLoad.m

IV. Visualization and checking

    Ploting results with and without colormap
    
    V1:: Truss_thickness_plot2D.m
    V2:: Truss_thickness_plot2D_colormap.m
    
    V3:: degradation_compliance_over_time_17_truss.m
    (used in Figue 8)
    
    V1:: SDTruss2D_AssembleBuckl.m
    (computes stress in each element for given force vector + checks the local buckling conditions)

