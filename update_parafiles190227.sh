#!/bin/bash

cp parameters.dat parameters.dat.bak
cat parameters.dat.bak | sed 13d | sed 13i'## new_config ; 0:Old Config 1:New Config 2:New config all accept 3:Old config all accept' | sed 11i'## branch_mode !! 0:normal, 1:make branch\n0\n## branch_root ! make branch from this config. default:0\n0\n## branch_num !! number of branches to make \n0' #> parameters.dat


cp inputfile inputfile.bak 
cat inputfile.bak | sed 3i'## branch_use !! in which branch the simulation is performed \n0' #> inputfile

