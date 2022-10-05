#!/usr/bin/bash

ALGN=$1

for TAXA in *_recon.tre
do
	cat ${TAXA%_recon.tre}.fst >> ${ALGN}
done

