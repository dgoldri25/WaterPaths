#!/bin/bash

JAVA_ARGS="-Xmx1g -classpath MOEAFramework-1.16-Executable.jar"
NUM_SAMPLES=1000
METHOD=latin

RANGES_FILENAME=implementation_ranges.txt
OUTPUT_FILENAME=ImplementationLHS.txt
CSV_FILENAME=ImplementationLHS.csv
java ${JAVA_ARGS} org.moeaframework.analysis.sensitivity.SampleGenerator -m ${METHOD} -n ${NUM_SAMPLES} -p ${RANGES_FILENAME} -o ${OUTPUT_FILENAME}

sed 's/ /,/g' ${OUTPUT_FILENAME} > ${CSV_FILENAME}
