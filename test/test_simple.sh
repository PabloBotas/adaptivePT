#!/bin/bash

../adaptSpotEnergies --patient P15_base_plan_simple \
                     --cbct adapt_data_P15/week6_cbct.mha \
                     --vf adapt_data_P15/week6_xform.mha \
                     --outdir results \
                     --output_vf vf.dat
