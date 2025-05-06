This repository provides implementations of the Shooting Method, a numerical technique used to solve Boundary Value Problems (BVPs) by transforming them into Initial Value Problems (IVPs). It explores the synergy between root-finding algorithms and IVP solvers.

Specifically, this repository includes implementations combining:

-   **Newton-Raphson Method:** A powerful iterative root-finding algorithm used here to refine the initial guess for the missing initial condition. This is coupled with:
    -   **Euler's Method:** A first-order numerical method for solving IVPs.
    -   **Runge-Kutta (RK) Methods:** A family of higher-order numerical methods (you can specify the order, e.g., RK4) for more accurate IVP solutions.

-   **Regula Falsi (False Position) Method:** Another root-finding algorithm, known for its guaranteed convergence (though potentially slower than Newton-Raphson). This is also coupled with:
    -   **Euler's Method:** For solving the resulting IVPs.
    -   **Runge-Kutta (RK) Methods:** For solving the resulting IVPs with higher accuracy.

This repository can be useful for understanding and comparing the performance of different combinations of root-finding and IVP solving techniques within the context of the Shooting Method.
