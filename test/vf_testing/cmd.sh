#!/bin/bash

plastimatch synth --pattern sphere \
                  --dim "3 3 3" \
                  --spacing "1 1 1" \
                  --origin "1 1 1" \
                  --sphere-center "2 2 2" \
                  --sphere-radius "0.5" \
                  --foreground "50" \
                  --background "25" \
                  --output "synth.mha"

# plastimatch synth --pattern sphere \
#                   --dim "101 101 101" \
#                   --spacing "1 1 1" \
#                   --origin "0 0 0" \
#                   --sphere-center "50.5 50.5 50.5" \
#                   --sphere-radius "25 25 25" \
#                   --foreground "50" \
#                   --background "25" \
#                   --output "synth.mha"
