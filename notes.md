# Optimizing BLS12-381 pairings

**Step 1**: Measure the time for the current implementation of the Tate pairing. 
Mainly we need to measure the Miller loop (function miller_loop_tate_pairing) 
and the final exponentiation (function: final_exponentiation). 
