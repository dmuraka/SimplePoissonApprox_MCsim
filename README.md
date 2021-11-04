# SimplePoissonApprox_MCsim
R code for the Monte Carlo experiments used Murakami and Matsui (2021) Improved log-Gaussian approximation for over-dispersed Poisson regression: application to spatial analysis of COVID-19. Arxiv 2104.13588.

As explained in this paper, the proposed closed-form approximation is fast, practical, and free from the identification problem, which can be severe for the conventional 
Poisson regression for small samples with many zeros. The Monte Carlo experiments for the conevntional Poisson regression (PoissonApprox_MCsim_case1(basic)_code.R) and Poisson additive mixed model (PoissonApprox_MCsim_case2(spatial)_code.R) demonstrate accuracy of the proposed approximation. 
