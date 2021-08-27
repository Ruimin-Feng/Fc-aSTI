# Fc-aSTI: An improved asymmetric susceptibility tensor imaging model with frequency offset correction  
Fc-aSTI is a new STI reconstruction method based on the frequency model: _f_= D* X+f_offset

Fc-aSTI releases the symmetry constraint imposed on susceptibility tensors and eliminates the orientation-independent frequency shift to mitigate the bias between the current asymmetric STI model and the measured resonance frequency shift.  

This repository contains implementations of Fc-aSTI model and the compared STI reconstruction models (STI, Fc-STI, aSTI).    

## Required dependencies: 
1. [STI_Suite V3.0 toolbox](https://people.eecs.berkeley.edu/~chunlei.liu/software.html)
2. Nifti toolbox
3. [Parallel computing toolbox](https://www.mathworks.com/products/parallel-computing.html)  
Additionally, maybe you need to install [FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki) on your system to generate tissue mask and perform registration in the phase preprocessing steps.  

## Usage  
simulation_demo.m is a demo of simulation experiments in the paper. It includes the following steps:  
(1) generate the ground truth susceptibility tensor elements, mean magnetic susceptibility (MMS), magnetic susceptibility anisotropy (MSA), principal eigenvector (PEV), and the simulation phase data with 12 different head orientations with a certain level of Gaussian noise.   
(2) run Fc-aSTI and the compared models (STI, Fc-STI, and aSTI) on the simulation data.   
(3) derive MMS, MSA, PEV with the recontructed susceptibility tensor.   
(4) calculate quantitative evaluation metrics (RMSE, SSIM, and AE) of the results obtained by the four models.  
One can run simulation_demo.m successfully as long as Nifti toolbox is installed and added to the path.  

For data acquired with 3D multi-echo GRE sequences, you should do the following preprocessing steps before running STI reconstruction codes in this repository  
(1) extract the tissue mask from magnitude images using FSL Bet    
(2) unwrap the raw phase data using Laplacian-based phase unwrapping in STI_Suite V3.0 toolbox.  
(3) remove the background phase using VSHARP in STI_Suite V3.0 toolbox. (obtain phase_tissue).  
(4) coregister the magnitude images at different orientations to a reference orientation (supine position) using FSL FLIRT.
(5) Apply the transform matrix to the corresponding tissue phase and calculate magnetic field direction. (obtain H_Matrix)  
(6) Once the above steps are completed, you can run run_models.m function:     

 %%%%%%%Example%%%%%%%%%%%%%%%%  
 addpath('utils');  
 phase_tissue=phase_tissue/(2pi×gamma×B0×TE)     %Normalized field perturbation  
 run_models(phase_tissue, H_Matrix, mask, folder);    
 % This function runs STI, Fc-STI, aSTI, and Fc-aSTI models.     
 % Parallel Computing Toolbox is used to speed up inversion process.  
 % Inputs:   
 % phase_tissue: N1×N2×N3×N_direction, phase data  
 % H_Matrix: N_direction × 3, each row represents unit direction vector of the applied magnetic field  
 % mask: N1×N2×N3, tissue mask  
 % folder: output directory  
 



