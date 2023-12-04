# Optimizing BLS12-381 pairings

**Step 1**: Measure the time for the current implementation of the Tate pairing. 
Mainly we need to measure the Miller loop (function miller_loop_tate_pairing) 
and the final exponentiation (function: final_exponentiation). 
It would be good to measure each of the two functions in separate, but also together. 
SO we would have a table like the following, where we can report the timings: 

Pairing type | Miller loop | Final expo | Total pairing |
:------------|:-----------:|:----------:|--------------:|
Tate pairing v1 | | | |
