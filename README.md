# btip_divsim
The scripts and code used for a diversity simulation project: "How low can we go?" Evaluating minimum sample size requirements for estimating nucleotide diversity via population genetic simulations.

Abstract:
Diversity simulations have become a recent technique for analyzing and modeling population-level genetic diversity and demographic history, eliminating the need for physical specimen collection. However, computational expense increases substantially as model parameters scale. We present an empirical workflow and statistical analysis that identify adequate sampling parameters for estimating nucleotide diversity in an isolated, non-migrating population model. This was achieved by testing combinations of sample size and locus count against a known baseline population, using the popgenART wrapper, fastsimcoal2, pixy, and ARTIllumina to simulate realistic sequencing read error. Two-way ANOVA confirms that individual sample size has a highly significant effect, while locus count and the interaction between the two factors showed no effect on the estimated nucleotide diversity. Descriptively, estimated diversity was observed to converge toward the true baseline value as both sample size and locus count increased, following an increasing asymptotic curve. These findings suggest an achievable balance between computational cost and sampling effort that future researchers may apply when designing diversity studies under similar resource constraints.

Co-interns: Carmel Saloma, Kyla Perocho, Marvin Bigay

Supervisors: Khylle Perdon, Dr. Ambrosio Matias

The creators used the popgenART wrapper, Art_illumina, fastsimcoal2, and pixy for this project. The user may user to these links and documentations to explore these tools:
1. https://manpages.debian.org/testing/art-nextgen-simulation-tools/art_illumina.1.en.html
2. https://github.com/demboc/popgenart/blob/main/tutorial.md
3. https://cmpg.unibe.ch/software/fastsimcoal2-25221/man/fastsimcoal25.pdf
4. https://pixy.readthedocs.io/en/latest/
