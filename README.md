# Rotary Excitation of non-sinusoidal pulsed magnetic fields: Toward noninvasive direct detection of cardiac conduction
Albertova<sup>1,2</sup>, M. Gram<sup>1,2</sup>, P. M. Blaimer<sup>3</sup>, W.R. Bauer<sup>1</sup>, P.M. Jakob<sup>2</sup>, P. Nordbeck<sup>1</sup>  
1 Department of Internal Medicine I, University Hospital Würzburg, Würzburg, Germany  
2 Experimental Physics 5, University of Würzburg, Würzburg, Germany  
3 Fraunhofer Institute for Integrated Circuits IIS, Würzburg, Germany  

# Address for correspondence:
Petra Albertova, University of Würzburg,  
Experimental Physics 5,  
Am Hubland, Würzburg,  
DE D-97074,  
petra.albertova@uni-wuerzburg.de,  
https://orcid.org/0000-0003-3646-7019

# Introduction
This folder contains exemplary Pulseq code for the detection of sinusoidal or non-sinusoidal pulsed magnetic fields via Rotary EXcitation (REX) based MRI. The measurements were carried out on the following MRI system:
- Siemens MAGNETOM Skyra, 3.0T
- max. gradient strength = 45mT/m
- max. slew rate = 200T/m/s

# Software requirements
- we used Matlab R2022a
- git clone https://github.com/pulseq/pulseq (our version was pulled on 20.01.2024)
- git clone https://github.com/petraalbertova/pulseq_rex_detection
- optional
	- we used the variable density spiral design functions provided by Brian Hargreaves; if you want to change the spiral readouts provided in this project, pull the original code from: http://mrsrl.stanford.edu/~brian/vdspiral/vdspiral.tar.gz
	- we used Shinnar–Le Roux (SLR) optimized pulses; if you want to change parameters, install the SIGPY package: pip install sigpy
	
# Hardware requirements
Please make sure to check the system specifications on your MRI system and adjust critical parameters (maximum gradient strength and slew rate). Run the timing and gradient checks prior to execution on the scanner! seq.checkTiming() and seq.testReport()
