#!/bin/bash
NSEEDS=10
SEEDS=$(seq 1 ${NSEEDS})

for SEED in ${SEEDS}
do
	awk 'BEGIN {FS=" "}; /^#/ {print $0}; /^[^#/]/ {printf("%s %s %s %s %s\n",$111,$112,$113,$114,$115)}' ./runtime/UrbanHeat_S${SEED}.runtime \
	 	> ./objs/UrbanHeat_S${SEED}.obj

done