#!/bin/bash

seed=${RANDOM}

echo "Launching Opt4D optimization"

planfile="$1"
output="${2//\//}"
opt4D -o "${output}/" \
      --lbfgs --lbfgs_m 20 \
      --linesearch_steps 20 \
      --dont_project_on_bounds \
      --unit_initial_fluence \
      --initial_fluence_noise 0.05 \
      --random_seed ${seed} \
      --add_step_event "UPDATE_LAGRANGE_MULTIPLIERS,20,20" \
      --add_step_event "UPDATE_PENALTY,20,20" \
      --constraint_penalty_multiplier 2  \
      --max_steps 2000 --min_steps 600  \
      --max_time 1200 \
      --write_dose --write_beam_dose \
      ${planfile}
