#!/bin/bash

# Strict clock replacement string
STRICT_CLOCK_STRING='\t\t\t<branchRateModel id="StrictClock.c:H3N2" spec="beast.evolution.branchratemodel.StrictClockModel" clock.rate="@clockRate.c:state_dependent_clock_sim"/>'

#make new xml file for strick clock
STRICT_CLOCK_XML="${1%.xml}_strict_clock.xml"
cat $1 > $STRICT_CLOCK_XML

#then replace state clocks branch model with strict clock
if [ $1 != "" ]; then
  sed -i '/.*branchRateModel.*/,/.*branchRateModel>/c\'"$STRICT_CLOCK_STRING"'' $STRICT_CLOCK_XML

fi
