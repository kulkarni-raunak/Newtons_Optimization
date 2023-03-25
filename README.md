# Newtons_Optimization
 ## Implement Newton's method in C/C++ and Python using below iteration step

$$
\begin{aligned}
J_{F}\left(x^{[i]}\right) \Delta x^{[i]} & =F\left(x^{[i]}\right) \\
x^{[i+1]} & =x^{[i]}-\Delta x^{[i]}
\end{aligned}
$$

## Using the following logic to generate Jacobian
$$
\left.\frac{\partial F_{i}(x)}{\partial x_{j}}\right|_{x=x x^{[i]}}=\frac{F_{i}\left(x^{[i]}+h e_{j}\right)-F_{i}\left(x^{[i]}\right)}{h}
$$

## The function that need to be minimised
$$
F\left(x_{1}, x_{2}, x_{3}\right)=\left(\begin{array}{c}
\frac{x_{1}}{x_{2}}+\frac{x_{3}}{x_{1}} \\
\frac{1}{2} x_{2}^{3}-250 x_{2} x_{3}-75000 x_{3}^{2} \\
e^{-x_{3}}+x_{3} \cdot e^{1}
\end{array}\right)
$$
