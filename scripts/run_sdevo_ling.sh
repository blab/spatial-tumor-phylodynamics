#!/bin/sh


#Ling et al
java -jar SDevo.jar xmls/hcc-wes_unidir_state_rep1.xml
java -jar SDevo.jar -seed 11 xmls/hcc-wes_unidir_state_rep2.xml
java -jar SDevo.jar -seed 210 xmls/hcc-wes_unidir_state_rep3.xml
java -jar SDevo.jar -seed 44 xmls/hcc-wes_unidir_state_strict_clock_rep3.xml
java -jar SDevo.jar -seed 33 xmls/hcc-wes_unidir_state_strict_clock_rep2.xml
java -jar SDevo.jar xmls/hcc-wes_unidir_state_strict_clock_rep1.xml
