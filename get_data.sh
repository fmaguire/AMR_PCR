#!/bin/bash

# get the primer data from EUCAST
mkdir -p data/primers
wget -P data/primers https://www.eurl-ar.eu/CustomerData/Files/Folders/25-resourcer/333_primerliste-til-web-07-11-2013-9.pdf

# tidy up the primer data
# manually use "tabula" to convert to tsv format

# get CARD data
mkdir -p data/CARD_prevalence
wget -P CARD_prevalence https://card.mcmaster.ca/download/6/prevalence-v3.0.1.tar.gz 
cd data/CARD_prevalence
tar xvf prevalence-v3.0.1.tar.gz
cd ../..

mkdir -p data/CARD
wget -P data/CARD https://card.mcmaster.ca/download/0/broadstreet-v2.0.2.tar.gz
cd data/CARD
tar xvf broadstreet-v2.0.2.tar.gz
cd ../..
