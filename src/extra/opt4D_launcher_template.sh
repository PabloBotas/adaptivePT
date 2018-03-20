#!/bin/bash

seed=${RANDOM}

echo "Launching Opt4D optimization"

planfile="${1}"
output="${2}"
cd ${output}
opt4D -o ./ \
      --lbfgs --lbfgs_m 20 \
      --linesearch_steps 20 \
      --dont_project_on_bounds \
      --unit_initial_fluence \
      --random_seed ${seed} \
      --add_step_event "UPDATE_LAGRANGE_MULTIPLIERS,20,20" \
      --add_step_event "UPDATE_PENALTY,20,20" \
      --constraint_penalty_multiplier 2  \
      --max_steps 50 \
      --max_time 1200 \
      --write_dose \
      ${planfile} | tee reoptimization.log

# opt4D -o ./ \
#       --lbfgs --lbfgs_m 20 \
#       --linesearch_steps 20 \
#       --dont_project_on_bounds \
#       --unit_initial_fluence \
#       --random_seed ${seed} \
#       --add_step_event "UPDATE_LAGRANGE_MULTIPLIERS,20,20" \
#       --add_step_event "UPDATE_PENALTY,20,20" \
#       --constraint_penalty_multiplier 2  \
#       --stop_fract 0.005 \
#       --stop_fract_memory 5 \
#       --max_time 120 \
#       --write_dose \
#       ${planfile} | tee reoptimization.log

cd - > /dev/null
