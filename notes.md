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
*Remark 1*. One expensive operation in the above algorithms is raising an element to the exponent $p$, where the size of $p$ is 381-bit. 
In order to speedup these exponentiations, we use the $\texttt{frobenius()}$ function that exists in the Sage $\texttt{Integer}$ folder. 
In particular we have the following equivalence: 

```r
from sage.all import Integer

f^(p^i) = f.frobenius(i) 
```

for all $i = 1, 2, \ldots, k - 1$, where in this case $k = 12$. 
Replacing the exponentiation to powers of $p$ with the $\texttt{frobenius()}$ function will speedup the easy part and hard part of the final exponentiation. 

*Remark 2*. Another expensive operation in the hard part of the final exponentiation is computing the inverse of an element in $\mathbb{F}_{p^{12}}$. 
Our construction offers one important property to do this. In particular, computing the inverse of an element $f$ is equivalent to raising this element to the power $p^6$. 
In other words, we have: 

```r
1/f = f.frobenius(6) 
```
Applying this trick will also speedup the final exponentiation.

## **Step 3: Addition/doubling and computation of lines** 

Usually in ECC implementations use the projective coordinate system to represent elliptic curve points and for executing the addition, doubling and line computation processes. 
The reason is that with projective coordinates we avoid the inverstions that appear both in addition and doubling. 
Recall that a point $P = (x_P, y_P)$ when written like this, it is in *affine form*, and satisfies the affine Weierstrass elliptic curve equation: 

$$E: y^2 = x^3 + ax + b$$

meaning that sabstituting the coordinates of $P$ in the equation of $E$ we get: $y_P^2 - x_P^3 - ax_P + b = 0$. 

A point on a curve in *projective form* has three coordinates and it is written as $P = (X_P : Y_P : Z_P)$. 
Then the corresponding projective Weierstrass elliptic curve equation is: 

$$E: Y^3Z = X^3 + aXZ^2 + bZ^3$$

meaning that sabstituting the coordinates of $P$ in the equation of $E$ we get: $y_P^2Z_P - x_P^3 - ax_PZ_P^2 + bZ_P^3 = 0$. 

What is relevant in pairing implementations is adding two points, doubling a point and computing the line functions $\ell$ and $v$. The algorithms for doing these computations are given below. 

### Doubling step

Given a point $P = (X_P : Y_P : Z_P)$ on a curve $E$, we want to compute the point $R = [2]P = (X_R : Y_R : Z_R)$. 
The coordinates of this new point are defined as follows: 

$$ X_R = 9X_P^4 - 8X_PY_P, \quad Y_R = 9X_R^3(4Y_R^2 - 3X_R^3) - 8Y_R^4, \quad Z_R = 2Y_RZ_R$$

In the Miller loop, in the doubling step, we need to compute the line $\ell_{R}(Q)$ (the tangent line that passes through $R$ and is evaluated at $Q = (x_Q, y_Q)$) and the line $v_R(Q)$ (the vertical line that passes through $R$ and is evaluated at $Q = (x_Q, y_Q)$).
Note that the second point $Q$ is in affine form, since we don't really do any computation on it. 
Then the lines $\ell_{R}(Q)$ and $v_R(Q)$ are:

$$\ell_R(Q) = 4Y_R^2 (2y_QY_RZ_R^3 - 3x_QX_R^2Z_R^2 + 3X_R^3 - 2Y_R^2), \quad v_R(Q) = 4x_QY_R^2Z_R^2 - 9X_R^4 + 8X_RY_R^2$$

The algorithm for efficiently computing the double of a point $R$ and the two lines $\ell_{R}(Q)$ and $v_R(Q)$ is the following (note that for BLS12-381, $a = 0$):

```r
DBL_step(R, Q, a = 0) {
  xR <- R[0], yR <- R[1], zR <- R[2]         
  xQ <- Q[0], yQ <- Q[1]
  T1 <- zR^2                                                                           
  A <- yR^2                                                                            
  B <- xR*A                            // B = xRyR^2                          
  C <- 3*xR^2                          // C = 3*xR^2 
  X3 = C^2 - 8*B                       // X3 = 9*xR^2 - 8*xR*yR^2  
  Z3 = (yR + zR)^2 - A - T1            // Z3 = 2*yR*zR
  Y3 = C*(4*B - X3) - 8*A^2            // Y3 = 9*xR^3(4*yR^2 - 3*xR^3) - 8*yR^4
  T3 = Z3^2                            // T3 = 4*yR^2*zR^2
  R = [X3, Y3, Z3]                          
  l = T3*(yQ*Z3 - C*xQ) + Y3 + C*X3    // l = 4*yR^2(2*yQ*yRzR^3 - 3*xQ*xR^2*zR^2 + 3*xR^3 - 2*yR^2)
  v = (xQ*T3 - X3)                     // v = 4*xQ*yR^2*zR^2 - 9*xR^4 + 8*xR*yR^2
  return(R, l, v)
}
```

### Addition step

Given a point $R = (X_R : Y_R : Z_R)$ and a point $P = (x_P, y_P)$ on a curve $E$, we want to compute the point $R = R + P = (X_R : Y_R : Z_R)$. 
Note that in the pairing computation, the point $P$ is in affine form, since we don't really do any computation on it. 
The coordinates of this new point are defined as follows: 

$$ X_R = 4(y_PZ_R^3 - Y_R)^2 - 4(x_PZ_R^2 + X_R)(x_PZ_R^2 - X_R)^2 $$ 

$$ Y_R = 8(x_PZ_R^2 + 2X_R)(y_PZ_R^3 - Y_R)(x_PZ_R^2 - X_R)^2 - 8(y_PZ_R^3 - Y_R)^3 - 8Y_R(x_PZ_R^2 - X_R)^3$$

$$ Z_R = 2Z_R(x_PZ_R^2 - X_R) $$

In the Miller loop, in the doubling step, we need to compute the line $\ell_{R,P}(Q)$ (the tangent line that passes through $R$ and $P$ and is evaluated at $Q = (x_Q, y_Q)$) and the line $v_{R + P}(Q)$ (the vertical line that passes through $R + P$ and is evaluated at $Q = (x_Q, y_Q)$).
Note that the second point $Q$ is in affine form, since we don't really do any computation on it. 
Then the lines $\ell_{R,P}(Q)$ and $v_{R + P}(Q)$ are:

$$\ell_{R, P}(Q) = (2y_QZ_R - 2y_PZ_R)(x_PZ_R^2 - X_R) - (2x_Q - 2x_P)(y_PZ_R^3 - Y_R)$$

$$v_{R + P}(Q) = 4x_QZ_R^2(x_PZ_R^2 - X_R)^2 -  4(y_PZ_R^3 - Y_R)^2 - 4(x_PZ_R^2 + X_R)(x_PZ_R^2 - X_R)^2$$

The algorithm for efficiently computing the sum of two points $R = (X_R : Y_R : Z_R)$ and $P = (x_P, y_P)$ and the two lines $\ell_{R,P}(Q)$ and $v_{R + P}(Q)$ is the following:

```r
ADD_step(R, P, Q) {
  xR <- R[0], yR <- R[1], zR <- R[2]
  xP <- P[0], yP <- P[1]     
  xQ <- Q[0], yQ <- Q[1]
  T1 <- zR^2
  B <- xP*T1                           // B = xP*zR^2      
  R2 <- yP^2                      
  D <- ((yP + zR)^2 - R2 - T1)*T1      // D = 2*yP*zR^3
  H <- B - xR                          // H = xP*zR^2 - xR
  I <- H^2                             // I = (xP*zR^2 - xR)^2
  E <- 4*I                             // E = 4*(xP*zR^2 - xR)^2
  J <- H*E                             // J = 4*(xP*zR^2 - xR)^3
  r <- D - 2*yR                        // r = 2*yP*zR^3 - 2*yR
  v <- xR*E                            // v = 4*xR*(xP*zR^2 - xR)^2
  X3 <- r^2 - J - 2*v                  // X3 = (2*yP*zR^3 - 2*yR)^2 - 4*(xP*zR^2 - xR)^3 - 8*xR*(xP*zR^2 - xR)^2
  Y3 <- r*(v - X3) - 2*yR*J            // Y3 = 8*(xP*zR^2 + 2*xR)(yP*zR^3 - yR)(xP*zR^2 - xR)^2 - 8*(yP*zR^3 - yR)^3 - 8*yR*(xP*zR^2 - xR)^3
  Z3 <- ((zR + H)^2 - T1 - I)          // Z3 = 2*zR*(xP*zR^2 - xR)
  T3 <- Z3^2                           // T3 = 4*zR^2*(xP*zR^2 - xR)^2
  R <- [X3, Y3, Z3]
  l <- (yQ - yP)*Z3 - (xQ - xP)*r      // l = (2*yQ*zR - 2*yP*zR)(xP*zR^2 - xR) - (2*xQ - 2*xP)(yP*zR^3 - yR)
  v <- (xQ*T3 - X3)                    // v = 4*xQ*zR^2*(xP*zR^2 - xR)^2 -  4*(yP*zR^3 - yR)^2 - 4*(xP*zR^2 + xR)*(xP*zR^2 - xR)^2
  return(R, l, v)
}
```
## **Step 4: Optimal ate pairing**

The most efficient type of pairings on elliptic curves is the *optimal ate* pairing. 
One can derive the optimal ate pairing from the tate pairing with just a few small modifications. 
We describe these modifications in the case of optimal ate pairing on the BLS12-381 curve. 

### Reverse the order of points

Let $P \in E(\mathbb{F}_p)[r]$ and $Q \in E(\mathbb{F}_q)[r]$, where $q = p^{12}$. 
That is, $P$ and $Q$ are both points of order $r$, but their coordinates are defined over different fields. 
In the case of the Tate pairing, the first input is the point $P$ and the second input is $Q$. 
This means that the doubling and addition operations involve the point $P$ and hence they are performed over $\mathbb F_p$, while the point $Q$ is used only to evaluate the lines $\ell$ and $v$. 
In the case of the optimal ate pairing, the input points are reversed. 
That is, the first input is the point $Q$ and the second input is $P$. 
This means that the doubling and addition operations will now involve the point $Q$ and hence they will be performed over the big field $\mathbb F_q$, while the point $P$ is used only to evaluate the lines $\ell$ and $v$. 

However, in the case of the optimal ate pairing, the point $Q$ needs to be fixed a bit differently. 
In addition to $Q$ having order $r$ on the big curve $E/F_q$, we also need to make sure that raising the coordinates of $Q$ to the power $p$ is equal to the point $[p]Q$. 
In other words, for the point $Q = (x_Q, y_Q)$, let $Q' = (x_Q^p, y_Q^p)$, we must have $Q' = [p]Q$. 
Fixing a point $Q$ such that it satisfies this property can be done with the following script: 

```r
Q = E12.random_element()  // Generate random point Q on the big curve E(Fq)
h12 = n12 / r ** 2        // Find the cofactor h12 of the big curve 
Q = h12 * Q               // This ensures that Q has order r
Q' = Q
xQ' = Q[0].frobenius()    // Raise x-coordinate of Q to the power p
yQ' = Q[1].frobenius()    // Raise y-coordinate of Q to the power p
Q' = E12(xQ',yQ')         // Set Q' as a point on the big curve E(Fq)
Q = Q' - Q                // Set Q = Q' - Q
}
```

Using this script, the new point $Q$ satisfies the property: $(x_Q^p, y_Q^p) = [p] (x_Q, y_Q)$. 
Note that the first three lines are the same as in the Tate pairing implementation. 

### Reduce the length of the Miller loop

The formula for computing the optimal ate pairing on BLS12 curves is: 

$$ e(Q, P) = f_{Q,|u|}(P)^{\dfrac{p^{12} - 1}{r}} $$

The above notation means that the length of the Miller loop of the optimal ate pairing for BLS12 curves depends on the binary representation of the seed $u$ and more precisely, because the seed can be negative, it depends on the binary representation of the absolute value of the seed, i.e. $|u|$.
The optimal ate pairing essentially requires much less doubling and addition steps. 
In particular, it requires $\log_2(|u|)-1$ doubling steps and $\text{hw}(|u|) - 1$ addition steps, as opposed to $\log_2(r)-1$ doubling steps and $\text{hw}(r) - 1$ addition steps in the case of the Tate pairing. 

Assuming that the size of $|u|$ is $\log_2(|u|) = n$ and its binary representation is $|u| = (b_{n}, b_{n - 1}, \ldots, b_1, b_0)$, the algorithm for computing the Miller function of the optimal ate pairing is given as follows. 

```r
optimal_ate_Miller_function_BLS12(Q, P, u)
  u0 <- abs(u)
  R <- Q
  f1 = 1                          // Two components for the Miller function: f1 and f2
  f2 = 1
  for i = b_{n - 1} to b_0 do: 
    R <- [2]R                     // Doubling depends on the point Q over the big field
    f1 <- f1^2 * l
    f2 <- f2^2 * v
    if b_i == 1 then:  
      R <- R + Q                  // Addition depends on the point Q over the big field
      f1 <- f1 * l
      f2 <- f2 * v
  if u < 0 then: f <- f2/f1      
  else: f <- f1/f2
  return(f)
```

Note that when the seed is negative, we need to compute the inverse of the Miller function. 

Alternatively, in order to avoid even this one inversion of the Miller function at the end when the seed is negative, we can do the following. 
When the seed is negative, we can set the point $Q$ as $-Q$ and then **we do not need any inversion of the Miller function in the end**. 
The algorithm for computing the optimal ate pairing in this case should be as follows: 

```r
optimal_ate_Miller_function_BLS12(Q, P, u)
  u0 <- abs(u)
  if u < 0 then: Q <- -Q = (xQ, -yQ)
  R <- Q
  f1 = 1                              // Two components for the Miller function: f1 and f2
  f2 = 1
  for i = b_{n - 1} to b_0 do: 
    R <- [2]R                         // Doubling depends on the point Q over the big field
    f1 <- f1^2 * l
    f2 <- f2^2 * v
    if b_i == 1 then:  
      R <- R + Q                      // Addition depends on the point Q over the big field
      f1 <- f1 * l
      f2 <- f2 * v
  f <- f1/f2
  return(f)
```
