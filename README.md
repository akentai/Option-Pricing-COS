# Multivariate COS Method for Basket Option Pricing

This repository presents a mathematical derivation and Python implementation of the **Multivariate COS method**, extending the approach of Fang & Oosterlee [1] and Ruijter & Oosterlee [2], and applying it to the pricing of **basket options**.

> ⚠️ This repository includes an original derivation of the **general multivariate COS coefficient formula** (Equation 10), not found in the referenced literature.

---

## 🧭 Project Overview

This project is split into two parts:

1. **Mathematical Derivation**
   - Formalizes the characteristic function for basket payoffs.
   - Generalizes the COS method to multivariate settings.
   - Derives COS coefficients using the characteristic function alone.
   - Applies the method to pricing a **put basket option**.

2. **Numerical Implementation** (coming soon)
   - Implements the multivariate COS method in Python.
   - Compares it against Monte Carlo simulation.

---

## 📌 Problem Setup

![Problem](docs/Equation1.png)

---

## 🧠 Solution Overview

![Approximation](docs/Equation2.png)

This formula enables **efficient multidimensional density recovery** from the characteristic function.

---

### COS Pricing Formula

Using the COS approximation, we write the final pricing formula for the basket option as:

![COS Pricing](docs/Equation3.png)

- \( F_k \): Fourier-Cosine coefficients from the characteristic function \( \phi(t) \)
- \( V_k \): Fourier-Cosine coefficients of the **payoff function**

---

## 📈 Monte Carlo Reference

For comparison, we also implement:

![Monte Carlo](docs/Equation4.png)

---

## 📚 References

[1] Fang, F. and Oosterlee, C.W.,  
*A novel pricing method for European options based on Fourier cosine series expansions*,  
SIAM J. Sci. Comput., 2009.  
[https://doi.org/10.1137/080718061](https://doi.org/10.1137/080718061)

[2] Ruijter, M. and Oosterlee, C.W.,  
*Two-dimensional Fourier cosine series expansion method for pricing financial options*,  
SIAM J. Sci. Comput., 2012.  
[https://doi.org/10.1137/120862053](https://doi.org/10.1137/120862053)

---

## 🗂️ File Structure (planned)
👉 For full derivation with LaTeX-rendered equations, see [docs/Derivation.pdf](docs/Derivation.pdf)

👉 To run the full code, open [COS_Basket.ipynb](COS_Basket.ipynb) 
