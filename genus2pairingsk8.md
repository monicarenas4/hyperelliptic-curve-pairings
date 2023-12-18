# Genus 2 pairings for $k = 8$

## General notes on genus 2 hyperelliptic curves

Elliptic curves are genus 1 curves and in the context of cryptographic applications they are represented in short Weierstrass form over a prime field $\mathbb F_p$ as:

$$ E/\mathbb F_p: y^2 = f(x) $$

where $f(x) = x^3 + ax + b$. Note that the degree of the polynomial $f$ is $\deg(f) = 3$. 

A genus 2 hyperelliptic curve $C$ over a prime field $\mathbb F_p$ is represented by an equation of the form 

$$C/\mathbb F_p: y^2 = f(x)$$

where in the cases we are interested in, the polynomial $f$ has the form $f(x) = x^5 + ax^4 + bx^3 + cx^2 + dx + e$, hence has degree $\deg(f) = 5$ and coefficients in $\mathbb F_p$. 

In the case of elliptic curves, the set $E(\mathbb F_p)$ defines a group. 
That is, the set of all points on the curve, together with the point at infinity $\mathcal O$, forms a group and the operation in the group is point addition. 
We write this group as: 

$$ E(\mathbb F_p) = \lbrace P = (x, y) \in \mathbb F_p \times \mathbb F_p : y^2 \equiv (x^3 + ax + b) \bmod p \rbrace \cup \lbrace \mathcal O \rbrace $$

For genus 2 curves, we assume the same set; i.e. the set of all points $P = (x, y)$ that satisfy the equation of the genus 2 curve. 
This is written as follows:

$$ C(\mathbb F_p) = \lbrace P = (x, y) \in \mathbb F_p \times \mathbb F_p : y^2 \equiv f(x) \bmod p \rbrace \cup \lbrace \mathcal O \rbrace $$

The problem with genus 2 is that such a set **is not** a group. 
Then if we want to develop DLP-based protocols using genus 2 curves, we need a different object associated to the curve $C$. 
We call this *special* object the *Jacobian of a genus 2 curve* and we denoted as $J(\mathbb F_p)$.   

### The Jacobian 

The elements in the Jacobian **are not points on the curve**, but they are derived from two points on the curve. 
In particular, let $P_1 = (x_1, y_1)$ and $P_2 = (x_2, y_2)$ be two points on the curve $C$. 
We compute the polynomial $u(x)$ as: 

$$ u(x) = (x - x_1)(x - x_2) = x^2 - (x_1 + x_2)x + x_1x_2 = x^2 + u_1x + u_0 $$

where $u_1 = - (x_1 + x_2)$ and $u_0 = x_1x_2$. 
Then we search for a polynomial $v(x)$ such that it satisfies both relations: 

$$ v(x_1) = y_1 \quad \text{and} \quad v(x_2) = y_2 $$

Solving this system yields a polynomial $v(x) = v_1x + v_0$, where 

$$ v_1 = \dfrac{y_2 - y_1}{x_2 - x_1} \quad \text{and} \quad v_0 = \dfrac{(x_2 - x_1)y_1 - (y_2 - y_1)x_1}{x_2 - x_1} $$

Then the element $D = [u(x), v(x)]$ is an element in the Jacobian $J(\mathbb F_p)$. 
This representation $D = [u(x), v(x)]$ is called *Mumford representation*. 

So, when we choose a random point in the Jacobian $J(\mathbb F_p)$, this is actually a pair of polynomials $D = [u(x), v(x)]$ where each polynomial has coefficients in $\mathbb F_p$. 
A few remarks based on the above definition of Jacobian elements: 
- The polynomial $u(x)$ is *monic*. This means that the leading coefficient (in other words, the coefficient of $x^2$) will always be 1.
- The polynomial $u(x)$ divides the polynomial $f(x) - v(x)^2$.
- The degree of $u(x)$ is $\deg(u) \leq 2$ and the degree of $v(x)$ is $\deg(v) < \deg (u)$. In other words, the polynomial $u(x)$ can be quadratic, linear, or constant, while the polynomial $v(x)$ will always be one degree less than $u(x)$.

Using this representation of elements in a Jacobian $J(\mathbb F_p)$, we can define an addition law in $J(\mathbb F_p)$, making it a group. 

Now we study how to define pairings on Jacobians of genus 2 curves. 

## A concrete example of a pairing an a genus 2 curve with embedding degree $k = 8$

