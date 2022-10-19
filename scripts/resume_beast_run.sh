#!/bin/bash
module load Java/15.0.1
LOGFILE=$(basename "$1" .xml)
if test -f "eden/logs/${LOGFILE}.log"; then
    java -jar Nab7.jar -resume $1
else
    java -jar Nab7.jar $1
fi

