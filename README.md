# Rotary Excitation of non-sinusoidal pulsed magnetic fields: Toward noninvasive direct detection of cardiac conduction
Albertova<sup>1,2</sup>, M. Gram<sup>1,2</sup>, P. M. Blaimer<sup>3</sup>, W.R. Bauer<sup>1</sup>, P.M. Jakob<sup>2</sup>, P. Nordbeck<sup>1</sup>  
<sup>1</sup> Department of Internal Medicine I, University Hospital Würzburg, Würzburg, Germany  
<sup>2</sup> Experimental Physics 5, University of Würzburg, Würzburg, Germany  
<sup>3</sup> Fraunhofer Institute for Integrated Circuits IIS, Würzburg, Germany  
doi: 10.1002/mrm.30190

# Address for correspondence:
Petra Albertova, University of Würzburg,  
Experimental Physics 5,  
Am Hubland, Würzburg,  
DE D-97074,  
petra.albertova@uni-wuerzburg.de,  
https://orcid.org/0000-0003-3646-7019

# Software requirements
- we used Matlab R2022a
- git clone https://github.com/petraalbertova/pulseq_rex_detection
- git clone https://github.com/pulseq/pulseq (we tested with the version from 10 June 2024)
- optional
	- we used the variable density spiral design functions provided by Brian Hargreaves; if you want to change the spiral readouts provided in this project, pull the original code from: http://mrsrl.stanford.edu/~brian/vdspiral/vdspiral.tar.gz
	- we used Shinnar–Le Roux (SLR) optimized pulses; if you want to change parameters, install the SIGPY package: pip install sigpy
	
# Hardware requirements
Please make sure to check the system specifications on your MRI system and adjust critical parameters (maximum gradient strength and slew rate). Run the timing and gradient checks prior to execution on the scanner! seq.checkTiming() and seq.testReport()  
Our measurements were carried out on the following MRI system:
- Siemens MAGNETOM Skyra, 3.0T
- max. gradient strength = 45mT/m
- max. slew rate = 200T/m/s

# Introduction
This folder contains exemplary Pulseq code for the detection of sinusoidal or non-sinusoidal pulsed magnetic fields via Rotary EXcitation (REX) based MRI.
- ahp_pulse.mat contains the pulseq object for an adiabatic half passage excitation pulse which was use for the spin-lock preparation
- spiral_readout.mat contains the gx and gy gradients which were use in our paper
- qrs.mat contains the signal of a cardiac QRS complex
Start Matlab, include the code from the official pulseq repository and open pulseq_main.m
1) enter your MRI system specifications
2) choose parameters for slice excitation or install the sigpy package for SLR pulses
3) keep our spiral readout or define your own with the Hargreaves vds toolbox. We used balanced gy and gy gradients and spoiling in z direction.
4) choose parameters for spin-locking. adapt the spin-lock time tSL and spin-lock amplitude fSL for your application depending on the magnetic stimuli.
5) define the waveforms of the magnetic stimuli. you can select "sine", "gauss", "sinc" or "qrs" and change the spectral properties e.g. with the oscillation frequency or time bandwidth products or the QRS duration.
6) choose a recovery time Trec
7) set the loop structure of your sequence. in our exemplary code, we added 5 dummy magnetization resets to obtain a steady-state prior to spin-locking. a detection experiment uses different relative timings between spin-lock and stimulus. you can change the number of experiments with the number of interactions Ninter. for an acceleration of the sequence, you can perform all interleaves of the spiral acquisition directly after spin-locking with ramped flip angles.
8) read the test report of the sequence to check gradient and slew rate violations.
9) copy the external.seq file to your scanner and perform the scan without changing the slice position.
10) reconstruct the data by using the kspace trajectories calculated by pulseq. we used the fessler toolbox and espirit.
11) calculate a REX detection map by calculating the pixel-wise standard deviation of the different REX weighted images   
