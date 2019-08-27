#!/bin/bash

for ii in A002_Xb945f7_X1551.SCI.J1851+0035 A002_Xb945f7_X1551.CAL.J1924-2914; do
	echo $ii
	caom2-repo delete --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo ALMACA $ii
done

for ii in $( ls actual*.xml); do
	echo $ii
	caom2-repo create --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo $ii
done
