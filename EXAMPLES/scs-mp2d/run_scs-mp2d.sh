#!/bin/bash

# From PSI4 MP2/aug-cc-pVTZ 
# @DF-RHF Final Energy:  -230.78162420426395
# Same-Spin Energy          =      -0.2399095404665765 [Eh]
# Opposite-Spin Energy      =      -0.7228154950870085 [Eh]


# Reference value:
# ! Total SCS-MP2D energy (hartrees)    =  -231.58762778 


export MP2D_PARAM_PATH=/home/gberan/Dropbox/mp2d_2024
../../MP2D benzene.xyz --scsmp2d --hf=-230.78162420426395 --mp2ss=-0.2399095404665765 --mp2os=-0.7228154950870085 > benzene.scsmp2d

echo "Computed value :" `grep '! Total SCS-MP2D' benzene.scsmp2d | awk '{print $NF}'`
echo "Reference value: -231.58762778"
