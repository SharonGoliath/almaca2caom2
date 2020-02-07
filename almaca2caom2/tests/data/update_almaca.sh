#!/bin/bash

for ii in data/*.actual.xml; do
  echo $ii
  x="${ii/\.actual\.xml/}"
  y="${x/data\//}"
  caom2-repo delete --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo ALMACA $y
  caom2-repo create --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo $ii
done
# caom2-repo delete --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo ALMACA A001_X88b_X23
# caom2-repo create --cert $HOME/.ssl/cadcproxy.pem --resource-id ivo://cadc.nrc.ca/sc2repo data/A001_X88b_X23.actual.xml 
