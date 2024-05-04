# Fc-aSTI: An improved asymmetric susceptibility tensor imaging model with frequency offset correction
paper link: https://onlinelibrary.wiley.com/doi/full/10.1002/mrm.29494 
Fc-aSTI is a new STI reconstruction method based on the frequency model: _f_= D* X+f_offset

Fc-aSTI releases the symmetry constraint imposed on susceptibility tensors and interprets the orientation-independent part of non-susceptibility-related frequency shifts.  

This repository contains implementations of Fc-aSTI model and the compared STI reconstruction models (STI, Fc-STI, aSTI).    

## Required dependencies: 
1. [STI_Suite V3.0 toolbox](https://people.eecs.berkeley.edu/~chunlei.liu/software.html)
2. Nifti toolbox
3. [Parallel computing toolbox](https://www.mathworks.com/products/parallel-computing.html)  
Additionally, maybe you need to install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) on your system to generate tissue mask and perform registration in the phase preprocessing steps.  

## Usage  
gen_simdata.m is a demo of generating the simulation data in the paper. 
One can run gen_simdata.m successfully as long as Nifti toolbox is installed and added to the path.      
simdata_demo_run.m is a script to run STI, aSTI, Fc-STI, and Fc-aSTI on the simulation data and evaluate the reconstruction results.  

For data acquired with 3D multi-echo GRE sequences, you should do the following preprocessing steps before running STI reconstruction codes in this repository  
(1) extract the tissue mask from magnitude images using FSL Bet    
(2) unwrap the raw phase data using Laplacian-based phase unwrapping in STI_Suite V3.0 toolbox.  
(3) remove the background phase using VSHARP in STI_Suite V3.0 toolbox. (obtain phase_tissue).  
(4) coregister the magnitude images at different orientations to a reference orientation (supine position) using FSL FLIRT.
(5) Apply the transform matrix to the corresponding tissue phase and calculate magnetic field direction. (obtain H_Matrix)  
(6) Once the above steps are completed, you can run different STI models provided in "utils" folder:     

 %%%%%%%Example%%%%%%%%%%%%%%%%  
 addpath('utils');  
 phase_tissue=phase_tissue/(2pi×gamma×B0×TE)     %Normalized field perturbation  
 Fc_aSTI_iteration(tissue_phase,H_Matrix,mask,offset_mask,iter_num,outpath) 
 
 %This function is to solve Fc_aSTI  according to the four steps of Figure 1 in the paper:  
 % An improved asymmetric susceptibility tensor imaging model with frequency offset correction  
 % (Feng et al)  
 % tissue_phase: H×W×D×NumDir, multi-orientation phase data  
 % H_Matrix: NumDir×3, direction vector of the magnetic field  
 % mask: H×W×D, brain mask  
 % offset_mask: H×W×D, mask that defines the regions where frequency offset will be considered.  
 % iter_num: number of iterations, 0 for no iteration  
 % outpath: outpath directory  
 



