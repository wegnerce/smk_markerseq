#!/bin/bash

cut -f1 $1 > OTU_IDs
cut -f4 $1 > raw_tax

cut -d : -f2 raw_tax | cut -d , -f1 | sed 's/\"//g' > domain
cut -d : -f3 raw_tax | cut -d , -f1 | sed 's/\"//g' > phylum
cut -d : -f4 raw_tax | cut -d , -f1 | sed 's/\"//g' > class
cut -d : -f5 raw_tax | cut -d , -f1 | sed 's/\"//g' > order
cut -d : -f6 raw_tax | cut -d , -f1 | sed 's/\"//g' > family
cut -d : -f7 raw_tax | cut -d , -f1 | sed 's/\"//g' > genus
cut -d : -f8 raw_tax | cut -d , -f1 | sed 's/\"//g' > species

paste OTU_IDs domain phylum class order family genus species > tax_temp

echo -e "\tDomain\tPhylum\tClass\tOrder\tFamily\tGenus\tSpecies" > header

cat header tax_temp > $2

rm OTU_IDs raw_tax domain phylum class order family genus species tax_temp header
