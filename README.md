# AMR PCR Primers

We are interested in how well different primers cover the AMR gene diversity
we now have access to via CARD prevalance.

Potentially this could be integrated as a test in CARD.

First things first is to get the CARD data and the primer information
from EUCAST tables using the `get_data.sh` script.

Unfortunately, as the primers are only distributed as the PDF of a table we
need to use optical character recognition to parse out the data (e.g. tabula)
once this is done (provided as `data/primers/tabula-333_primerliste-til-web-07-11-2013-9.tsv`)
we can parse it using the notebook and scripts in `amr_pcr.py`.

We then use this notebook and scripts to generate a bash script for 
running VA Mike's VAWare tool to tally the mismatches per AMR gene.

