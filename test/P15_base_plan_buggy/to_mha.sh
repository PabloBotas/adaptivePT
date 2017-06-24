#!/bin/bash

for f in out_base_plan/D*.dose
do
    echo ${f}
    name=$(basename ${f})
    dir=$(dirname ${f})
    outname=${name/.dose/.mha}

    nx=$(grep "nVoxelsX" out_base_plan/geometry.dat | awk '{print $2}')
    ny=$(grep "nVoxelsY" out_base_plan/geometry.dat | awk '{print $2}')
    nz=$(grep "nVoxelsZ" out_base_plan/geometry.dat | awk '{print $2}')
    /opt/gpmc-tools/default/gpmc2xio ${dir}/ast_${name} ${f} $nx $ny $nz 500000 1 >> /dev/null
    plastimatch convert --input-dose-ast ${dir}/ast_${name} \
                        --output-dose-img out_base_plan_mha/${outname} >> /dev/null
done