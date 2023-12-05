# Optimizing BLS12-381 pairings

## **Step 1**: Measure the time for the current implementation of the Tate pairing. 

Mainly we need to measure the Miller loop (function miller_loop_tate_pairing) 
and the final exponentiation (function: final_exponentiation). 
It would be good to measure each of the two functions in separate, but also together. 
SO we would have a table like the following, where we can report the timings: 

Pairing type | Miller loop | Final expo | Total pairing |
:------------|:-----------:|:----------:|--------------:|
Tate pairing v1 | | | |

## **Step 2: Final expo** 

Assume that the output of the Miller function is $f \in \mathbb{F}_{p^{12}}$. 
Then we need to compute $f^{(p^{12} - 1)/r}$, where $(p^{12} - 1)/r$ is an exponent of size 4317-bit.
Raising an element $ f $ to such an expoenent is a very expensive operation. 
Hence we need to apply certain tricks to speed-up this process. We write: 

$$ \dfrac{p^{12} - 1}{r} = \dfrac{(p^6 - 1)(p^2 + 1)(p^4 - p^2 + 1)}{r} $$

Now instead of raising $f$ to the left hand-side, we can raise it to the right hand-side. 
That is, wee need to compute the value: 

$$ f^{\dfrac{(p^6 - 1)(p^2 + 1)(p^4 - p^2 + 1)}{r}} = f^{(p^6 - 1)(p^2 + 1)\dfrac{p^4 - p^2 + 1}{r}} $$

In order to do this, we first raise $f$ to the exponent $(p^6 - 1)(p^2 + 1)$. 
We call this the *easy part* of the final exponentiation. 
Then we raise the result to the exponent $(p^4 - p^2 + 1)/r$. 
We call this the *hard part* of the final exponentiation. 
