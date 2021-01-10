#!/bin/bash
root -l <<EOF
gROOT->LoadMacro("ParticleType.cxx+")
gROOT->LoadMacro("ResonanceType.cxx+")
gROOT->LoadMacro("Particle.cxx+")
gROOT->LoadMacro("main.cxx+")
Generate()
gROOT->LoadMacro("Analysis.cxx+")
Analize()
EOF
root -l
#The script was tested using Ubuntu 18.04 (Both WSL and OS).
#The script should be located inside the folder containing the particle program. 
#Before running the script "chmod +x runROOT.sh" should be done, then use it normally by doing ./runROOT.sh   
#Replace accordingly the name of your files and function of your main ROOT macro
#After executing your program the ROOT shell will remain interactive in order for you to use a TBrowser and check your newly created plots
#Report any issue to marco.giacalone2@unibo.it