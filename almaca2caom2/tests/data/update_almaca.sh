#!/bin/bash

for ii in actual*.xml; do
  echo $ii
  x="${ii/actual_/}"
  y="${x/\.xml/}"
  caom2-repo delete --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo ALMACA $y
  caom2-repo create --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo $ii
done
