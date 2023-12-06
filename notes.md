# Optimizing BLS12-381 pairings

## **Step 1**: Measure the time for the current implementation of the Tate pairing. 

Mainly we need to measure the Miller loop (function miller_loop_tate_pairing) 
and the final exponentiation (function: final_exponentiation). 
It would be good to measure each of the two functions in separate, but also together. 
SO we would have a table like the following, where we can report the timings: 

Pairing type | Miller loop | Final expo | Total pairing |
:------------|:-----------:|:----------:|--------------:|
Tate pairing v1 | 0.114225 | 0.058757 | 0.172982 |

## **Step 2: Final expo** 

Assume that the output of the Miller function is $f \in \mathbb{F}_{p^{12}}$. 
Then we need to compute $f^{(p^{12} - 1)/r}$, where $(p^{12} - 1)/r$ is an exponent of size 4317-bit.
Raising an element $f$ to such an expoenent is a very expensive operation. 
Hence we need to apply certain tricks to speed-up this process. We write: 

$$ \dfrac{p^{12} - 1}{r} = \dfrac{(p^6 - 1)(p^2 + 1)(p^4 - p^2 + 1)}{r} $$

Now instead of raising $f$ to the left hand-side, we can raise it to the right hand-side. 
That is, wee need to compute the value: 

$$ f^{\dfrac{(p^6 - 1)(p^2 + 1)(p^4 - p^2 + 1)}{r}} = f^{(p^6 - 1)(p^2 + 1)\dfrac{p^4 - p^2 + 1}{r}} $$

In order to do this, we first raise $f$ to the exponent $(p^6 - 1)(p^2 + 1)$. 
We call this the *easy part* of the final exponentiation. 
Then we raise the result to the exponent $(p^4 - p^2 + 1)/r$. 
We call this the *hard part* of the final exponentiation. 

No the function final_exponentiation_BLS12 will be split in two parts, the **easy part** and the **hard part** and it will look like: 

```r
final_exponentiation_BLS12(f) {
  f <- final_exp_easy_k12(f)
  f <- final_exp_hard_BLS12(f)
}
```

We examine each of the functions final_exp_easy_k12 and final_exp_hard_BLS12 in separate. 

### Final exponentiation: easy part

On input the Miller function $f$, we need to compute the value: 

$$ f^{(p^6 - 1)(p^2 + 1)} $$

This is done with the function final_exp_easy_k12 as follows: 

```r
final_exp_easy_k12(f) {
  inv_f <- 1/f         // inv_f = f^(-1)
  f     <- f^(p^6)     // f = f^(p^6)
  f     <- f * inv_f   // f = f^(p^6) * f^(-1) = f^(p^6 - 1)
  f2    <- f^(p^2)     // f2 = f^(p^2) = [f^(p^6 - 1)]^(p^2)
  f     <- f2 * f      // f = f2 * f = [f^(p^6 - 1)]^(p^2) * [f^(p^6 - 1)] = f^[(p^6 - 1)(p^2 + 1)]
  return(f)
}
```
### Final exponentiation: hard part

On input the output $f$ of the function final_exp_easy_k12, we need to compute the value: 

$$ f^{\dfrac{p^4 - p^2 + 1}{r}} $$

By constraction and by the definition of the parameters, we know that $r$ divides $p^4 - p^2 + 1$. 
Then the exponent $(p^4 - p^2 + 1)/r$ can be written as follows: 

$$\dfrac{p^4 - p^2 + 1}{r} = \lambda_3 p^3 + \lambda_2 p^2 + \lambda_1p + \lambda_0$$

where the $\lambda_i$ are polynomials in the seed $u \in \mathbb{Z}$ for every $i = 0, 1, 2, 3$. 
Obtaining these $\lambda_i$ is out of the scope for now, we assume they are fixed as follows: 

$$ \lambda_3 = (u - 1)^2, \quad \lambda_2 = u\lambda_3, \quad \lambda_1 = u \lambda_2 - \lambda_3, \quad \lambda_0 = u \lambda_1 + 3 $$

Observe that the seed for the curve BLS12-381 is $u = -(2^{63} + 2^{62} + 2^{60} + 2^{57} + 2^{48} + 2^{16})$. 
Hence it is negative and can be written as $u = - u_0$, where $u_0 = 2^{63} + 2^{62} + 2^{60} + 2^{57} + 2^{48} + 2^{16}$. 
Thus, the above $\lambda_i$ can be written in terms of $u_0$ as follows: 

$$ \lambda_3 = (-u_0 - 1)^2 = (u_0 + 1)^2, \quad \lambda_2 = -u_0\lambda_3, \quad \lambda_1 = -u_0\lambda_2 - \lambda_3 = -(u_0\lambda_2 + \lambda_3), \quad \lambda_0 = - u_0 \lambda_1 + 3 $$

Based on this, we need to compute the exponentiation: 

$$ f^{\lambda_3 p^3 + \lambda_2 p^2 + \lambda_1p + \lambda_0} = f^{\lambda_0 + p (\lambda_1 + p (\lambda_2 + p \lambda_3))} $$

This is done using the following function: 

```r
final_exp_hard_BLS12(f, u0) {
  f3 <- f^(u0 + 1)         // f3 = f^(u0 + 1)
  f3 <- f3^(u0 + 1)        // f3 = f^(u0 + 1)^2 = f^(λ3)
  res <- f3^p              // res = f3^p = f^(λ3*p)
  f2 <- f3^u0              // f2 = f^[u0(u0 + 1)^2] = f^[u0*λ3]
  f2 <- 1/f2               // f2 = f2^(-1) = f^[-u0*λ3] = f^(λ2)           
  res <- res * f2          // res = f^(λ3*p) * f^(λ2) = f^(λ2 + pλ3)
  f1 <- f2^u0              // f1 = f2^u0 = f^(u0*λ2)
  f1 <- f1 * f3            // f1 = f1 * f3 = f^(u0*λ2) * f^(λ3) = f^[u0*λ2 + λ3]
  f1 <- 1/f1               // f1 = f1^(-1) = f^[-(u0*λ2 + λ3)] = f^(λ1)
  res <- res^p             // res = res^p = f^[p(λ2 + pλ3)]
  res <- res * f1          // res = res * f1 = f^[p(λ2 + pλ3)] * f^(λ1) = f^[λ1 + p(λ3*p + λ2)]
  f0 <- f1^u0              // f0 = f1^u0 = f^(u0*λ1)
  f0 <- 1/f0               // f0 = f0^(-1) = f^(-u0*λ1)
  f0 <- f0 * f^2 * f       // f0 = f0 * f^3 = f^(-u0*λ1 + 3) = f^(λ0)
  res <- res^p             // res = res^p = res^[p(λ1 + p(λ3*p + λ2))]
  res <- res * f0          // res = res * f0 = f^[p(λ1 + p(λ3*p + λ2))] * f^(λ0) = f^[λ0 + p(λ1 + p(λ3*p + λ2))]    
  return(res)
}
```
