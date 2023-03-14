
PVEDDIES (Optimization Suite):
----------------------------------------------------
Copyright (C) 2022 Paul Ernst, Bulusu Subrahmanyam,
		   Yves Morel, Alexis Chaigneau,
		   Corinne Trott

Primary Author: Paul Ernst
Special Thanks To: Le Vu, L'Hégaret, Nencioli, Pegliasco for code from other applications
Date Of Writing:    1/12/2023

Contact Information:

Paul Ernst:		pernst@seoe.sc.edu

Contents: 

This folder contains all of the necessary pieces to
fully train and run the PVEDDIES algorithms.

Quick Start Guide:

Step 1: You should start from the file Paul_OPT_RunOptimization.m. The first files to modify from there are Paul_OPT_CalibrateNEMODataStructure (so that it matches your file system's model input NC files) and Paul_OPT_CreateDefaultData (this should run without modifications, but ensure all paths are correct). Once these are properly configured you should have a DefaultData_OPT.mat in the style of the example provided. Make sure Paul_OPT_LoadDefaultSettings loads and applies all paths properly.

Step 2: If you have set up the file system correctly and you have the methods you want hooked into the files (default is as in published paper), you can run Paul_RunOptimization.m with your preferred isopycnals at the surface, preferred reference profile location, etc. This will provide you a set of figures as in the published paper (Fig. 1-5) and the optimized identification algorithm. Fig. 6-9 can be created by modifying Paul_OPT_DisplayMethodComparison.m. If you want to fiddle around with the methods for spatial smoothing, you can tinker with spSmoothing.m in the FUNCTIONS folder. As-is it's a basic moving average, flat kernel size. Depending on the speed of your computer, this might take a while (for a single Mac workstation, this took me 2 weeks to run; recommend running on a cluster or heavily reducing the number of days used from full run to just a couple months). Note that if you want to replace the "ground truth method" you will have to go into Paul_OPT_Key.m and Paul_OPT_EddyExtractionMaster.m and insert your method as an option with all relevant inputs.

Step 3: Once you have your optimized method, move to Paul_OPT_RunEddyTracking. You can use Paul_OPT_TS_RSW to select your final depth isopycnals or use some other source. We then use Paul_OPT_MethodPicker to select the tracking method with the parameters you specify. The EDDYTRACKING folder contains the Pegliasco code for cost functions and such. If you want to use your own tracking algorithm, you can hook up the outputs of our detection algorithm to the EDDYTRACKING folder. You can then modify the Paul_OPT_SSMOVIE_PV.m or Paul_OPT_SSMOVIE_SP.m to create nice movies of your IDs over time.

Step 4: You can use Paul_OPT_RunParticleTracking.m with the functions in the PARTICLE_TRACKING folder to obtain a case study like the one we display in our paper.

Each file has a short description of the purpose/inputs/outputs and the main user-tuned inputs are at the top of the file after the header. You may find you need to tinker around with specific values to get the results you want; the code is generally well commented in sections where a user might need to tune values.

PATH DESCRIPTIONS:
Folder COMPARISON:
	This is where each run's files will go for the optimization suite.
Folder EDDYTRACKING:
	This is where the files for the Eddy Tracking code (Pegliasco 2016 version) live.
Folder EXTRACTION:
	This is where the results of the full tracking live (individual files, eddy maps, trajectories)
Folder FIGURES: 
	This is where the figures from all of the files will live.
Folder FUNCTIONS:
	All minor helper functions and packages (m_map, gsw toolbox) live.
Folder PARTICLE_TRACKING:
	The Particle Tracking subfunctions live here.
Folder SIMILARITIES_FINAL:
	The final compiled .mat files for the optimization suite live here.

Feel free to contact me with any other questions.

Terms:	  

The source code (M) is provided under the
terms of the GNU General Public License.

--- end of the file ---