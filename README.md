# Environmental fluctuations and ecological interactions jointly constrain the adaptive potential of species

This repository contains the model, numerical simulations, and figures of the manuscript:

**Environmental fluctuations and ecological interactions jointly constrain the adaptive potential of species**

Leandro G. Cosmo and Jordi Bascompte. Department of Evolutionary Biology and Environmental Studies, University of Zurich.

The persistence of species under environmental change depends on their adaptive potential, the capacity of natural selection to increase or maintain the mean fitness of populations. Such adaptive potential can be depicted through the fitness landscape, which maps the mean trait value of a population into its mean fitness. For most species, these landscapes portray adaptations to both the abiotic environment and ecological interactions, such as mutualisms. However, fitness landscapes are not static. The abiotic environment fluctuates, and the traits of mutualistic partners change through reciprocal selection and coevolution. Here, we combine an analytical framework with 220 empirical mutualistic networks to show how environmental fluctuations and coevolution jointly constrain the adaptive potential of species.

Our results show three main ways through which environmental fluctuations and species interactions shape the long-term fitness of populations:

1. The long-term, geometric mean fitness (GMF) of a population depends on two components of its fitness landscape, the distance of the population to its adaptive peak and the height of this peak. The height of the peak sets the maximum mean fitness that the population can achieve.
2. When the traits of mutualistic partners fluctuate independently, generalists benefit through a portfolio effect, which increases the height of their adaptive peaks.
3. Coevolution correlates the trait fluctuations of mutualistic partners. This correlation constrains the portfolio effect and grows with indirect evolutionary effects, which are strongest in species-rich, nested networks.

Consequently, environmental fluctuations and coevolution may lower the ceiling of adaptation, even when natural selection keeps populations at their adaptive peaks.

## Modeling coevolution and fitness in mutualistic networks

Each species has a single continuous trait that mediates its mutualistic interactions, such as the depth of a floral tube or the length of a proboscis. The mean trait value of the population of species $`i`$, $`\overline{z}_i`$, affects adaptations through two sources in the environment. Mutualistic interactions contribute with a proportion $`m_i`$ and favor trait matching between species $`i`$ and each of its mutualistic partners $`j`$. The abiotic environment contributes with the remaining proportion, $`(1-m_i)`$, and favors an environmental optimum, $`\theta_i^{(t)}`$. Measuring fitness relative to its theoretical maximum, the mean population fitness of species $`i`$ at generation $`t`$ is (equation 1 of the main text):

$$
\overline{w}_i^{(t)} = \exp\left[\frac{m_i}{\alpha}\ln\left(\sum_j \frac{a_{ij}}{k_i}\, e^{-\alpha\left(\overline{z}_j^{(t)}-\overline{z}_i^{(t)}\right)^2}\right) - (1-m_i)\left(\theta_i^{(t)}-\overline{z}_i^{(t)}\right)^2\right],
$$

where $`a_{ij}=1`$ when species $`i`$ and $`j`$ interact and zero otherwise; $`k_i`$ is the number of mutualistic partners of species $`i`$; and $`\alpha`$ controls how strongly the benefits of mutualisms decay with mismatches in traits. Fitness decreases when the trait of species $`i`$ mismatches, on average, the traits of its mutualistic partners, and when it differs from its environmental optimum. We derive this function by integrating the selection gradient that drives coevolution in mutualisms across different studies (Supplementary Information, section S1).

## Mapping geometric mean fitness onto the fitness landscape

When the environment fluctuates, fitness changes from one generation to the next and accumulates multiplicatively. As such, the long-term growth of a population is governed by its geometric mean fitness, whose log equals the temporal average of log fitness, $`\ln G(\overline{w}_i)=\mathbb{E}_t[\ln\overline{w}_i^{(t)}]`$. Assuming that the benefits of mutualisms decay slowly with trait differences (small values of $`\alpha`$), the log of GMF can be written as a function of three components (equation 2 of the main text):

$$
\ln G(\overline{w}_i) \approx -\,\mathbb{E}_t\!\left[\left(\overline{z}_i^{(t)}-\overline{z}_{i,\mathrm{opt}}^{(t)}\right)^2\right] - \mathbb{E}_t\!\left[V_i^{(t)}\right] + \frac{m_i\alpha}{2}\,\mathbb{E}_t\!\left[V_{D,i}^{(t)}\right].
$$

The first term is the distance to the peak. It contains the squared distance between the mean trait of the population and its adaptive peak, $`\overline{z}_{i,\mathrm{opt}}^{(t)} = m_i \sum_j (a_{ij}/k_i)\,\overline{z}_j^{(t)} + (1-m_i)\,\theta_i^{(t)}`$, and depicts how maladapted the population is. The two remaining terms compose the height of the peak. The variance of selective targets, $`V_i^{(t)}`$, quantifies the trade-off in adaptation between mutualistic partners and the abiotic environment, and always decreases the height of the peak. In turn, the variance in trait matching between species $`i`$ and each of its partners, $`V_{D,i}^{(t)}`$, increases the height of the peak, because the well-matched partners compensate for the poorly matched ones.

This equation averages the fitness landscapes that a species experiences across generations and defines the geometric mean fitness landscape. It then follows that a fluctuating environment can hinder the long-term fitness of species in two distinct ways. On the one hand, a population may lag behind its moving peak and accumulate maladaptation. On the other hand, even a population that is at its adaptive peak in every generation can still achieve a smaller GMF if these adaptive peaks are, on average, low.

## The two scenarios

We evaluate the three components of GMF in two scenarios that differ in how the traits of mutualistic partners change over time:

| Scenario | Traits | Environmental optima | Script | Output table |
|---|---|---|---|---|
| Independent fluctuations | Partner traits drawn independently at every generation, $`\mathcal{N}(0,\,1.5^2)`$. The trait of the focal species remains fixed at $`\overline{z}_i=1`$ | Drawn independently at every generation, $`\mathcal{N}(0,\,1)`$ | `scripts/run_random_scenario.R` | `results/random_gmf_by_degree.csv` |
| Coevolution | The traits of all species reach their coevolutionary steady state at each time step | Drawn independently for every species at each time step, $`\mathcal{N}(0,\,2.88^2)`$ | `scripts/run_coevolution.R` | `results/coevolution_gmf.csv` |

Both scenarios use $`\alpha = 0.2`$ and $`m_i = 0.5`$, and the values above are the ones set in the scripts. The scenario with independent fluctuations depends on the network only through the number of partners of the focal species, so we store it once per number of partners. The coevolutionary scenario depends on the whole network, so we store it once per species.

### Independent fluctuations favor generalists through a portfolio effect

The scenario with independent fluctuations provides a baseline without coevolution, in which partners behave as independent components of a fluctuating environment. Under this scenario, GMF can be evaluated in a closed form (equation 3 of the main text):

$$
\ln G(\overline{w}_i) \approx -\left(\overline{z}_i - \mathbb{E}_t\!\left[\overline{z}_{i,\mathrm{opt}}\right]\right)^2 - \mathcal{V}_i + m_i\,\alpha\,\sigma_{Dz}^{2}\left(1-\frac{1}{k_i}\right).
$$

The first term depicts the mismatch between the trait of species $`i`$ and the average position of its moving peak. The term $`\mathcal{V}_i`$ is the long-run counterpart of the variance of selective targets. The final term is the portfolio effect, scaled by the temporal variance in the trait matching between species $`i`$ and its partners, $`\sigma_{Dz}^{2}`$. The factor $`(1-1/k_i)`$ equals zero for specialists with a single partner and increases toward one with the number of partners. When partner traits fluctuate independently, well-matched partners compensate for poorly matched ones in any given generation. The portfolio effect does not affect the distance between a population and its peak. Instead, it raises the adaptive peak, and generalists achieve a higher GMF than specialists because their maximum adaptive potential increases.

### Coevolution reshapes the adaptive potential of species

Coevolution couples trait dynamics and correlates the response of species to environmental changes. When species coevolve, their traits converge to a globally stable steady state (equation 4 of the main text):

$$
\mathbf{z}^{*} = (\mathbf{I}-\mathbf{Q})^{-1}\,\boldsymbol{\Psi}\,\boldsymbol{\theta} = \mathbf{T}\boldsymbol{\theta},
$$

where $`\mathbf{Q}`$ contains the direct evolutionary effects that mutualistic partners exert on each other, and $`\boldsymbol{\Psi}`$ is a diagonal matrix that weights the selective pressures of the abiotic environment. The matrix $`\mathbf{T}`$ accumulates the evolutionary effects that flow through all the direct and indirect coevolutionary pathways of the network. Each entry $`T_{ij}`$ measures how strongly the environmental optimum of species $`j`$ shapes the steady-state trait of species $`i`$. We assume that species traits reach this steady state before the environmental optima change. Coevolving species are then always near their adaptive peaks, and their GMF depends only on the temporal average of the height of these peaks (equation 5 of the main text):

$$
\ln G(\overline{w}_i) \approx -\,\sigma_\theta^{2}\left(S_i - m_i\,\bar{S}_j\right) + m_i\,\alpha\left(1-\frac{1}{k_i}\right)\bar{\nu}_i^{2}\left(1-r_i^{2}\right).
$$

The first term is the conflict of evolutionary effects, the coevolutionary counterpart of the variance of selective targets. Here, $`S_i = 1-\sum_k T_{ik}^2`$ is the Simpson diversity of the coevolutionary pathways that reach species $`i`$, and $`\bar{S}_j`$ is the average diversity of the pathways that reach its partners. When the pathways reaching species $`i`$ are more diverse than the ones reaching its partners, mutualistic and environmental selection favor different trait values, and the adaptive peaks of species $`i`$ become, on average, lower.

The second term is the portfolio effect of generalists under coevolution. The term $`\bar{\nu}_i`$ is the average temporal variance of the trait differences between species $`i`$ and its partners. In turn, $`r_i^{2}`$ is the temporal correlation between the squared trait differences of species $`i`$ with its partners, averaged across all pairs of partners. Thus, $`r_i^{2}`$ measures the extent to which the trait matching of species $`i`$ with its different partners fluctuates in synchrony across generations. When partners receive evolutionary effects from the same species through shared coevolutionary pathways, $`r_i^{2}`$ approaches one and the portfolio effect collapses, even for species with many partners. The Supplementary Information contains the complete derivation of these equations (sections S1 to S6).

## Indirect effects and the number of pathways connecting species

Both components of peak height under coevolution are built from the rows of the $`\mathbf{T}`$-matrix. Thus, they depend on how evolutionary effects propagate through the pathways of the network. We quantify the contribution of indirect effects of each network as the proportion of the total evolutionary effects in $`\mathbf{T}`$ that connects species that do not directly interact:

$$
F = \frac{\sum_i \sum_{j \neq i} (1-a_{ij})\,T_{ij}}{\sum_i \sum_{j \neq i} T_{ij}}.
$$

In turn, we measure the number of pathways connecting species in each network as $`1-1/\lambda`$, where $`\lambda`$ is the leading eigenvalue of the adjacency matrix. The number of pathways of length $`\ell`$ connecting any two species grows as $`\lambda^{\ell}`$. Consequently, the more pathways a network has, the easier it is for species that do not directly interact to affect each other. A high number of pathways is typical of the species-rich and nested networks that characterize pollination and seed-dispersal mutualisms. In these networks, indirect effects increase the conflict of evolutionary effects and the correlation among partners, and species with the same number of partners achieve a lower GMF.

## Empirical networks

We parameterize the model with the structure of empirical mutualistic networks from the Web of Life database (www.web-of-life.es). The repository contains 230 plant-by-animal incidence matrices in `src/empirical_networks/`. We set $`a_{ij}=1`$ whenever two species interact, drop species without interactions, and exclude the ten networks with 300 or more species. The 220 remaining networks compose the dataset of the analyses:

| Mutualism | Networks |
|---|---:|
| Plants and their pollinators | 165 |
| Plants and the frugivores that disperse their seeds | 33 |
| Anemones that protect anemonefishes | 17 |
| Plants with extrafloral nectaries and their protective ants | 5 |

Figure 4a,b isolates the effects of network structure from those of the number of partners using species with three mutualistic partners. Species with three partners are present in 203 of the 220 networks.

## Manuscript figures

| Figure | Result | Script | Output |
|---|---|---|---|
| Figure 1 | Fluctuations in the traits of mutualistic partners and in the environment reshape the fitness landscape of a population and determine its long-term fitness | `scripts/main_fig1.R` | `figures/main_fig1.png` |
| Figure 2 | Independent fluctuations of partner traits and the environment create a portfolio effect that increases the maximum adaptive potential of generalists | `scripts/main_fig2.R` | `figures/main_fig2.png` |
| Figure 3 | Coevolution reshapes the geometric mean fitness of species through the height of their adaptive peaks | `scripts/main_fig3.R` | `figures/main_fig3.png` |
| Figure 4 | Indirect evolutionary effects link the structure of mutualistic networks to the adaptive potential of species under environmental fluctuations | `scripts/main_fig4.R` | `figures/main_fig4.png` |

Each script also saves an SVG version of its figure. Because coevolving species track their moving adaptive peaks, the landscapes of Figure 3 are centered on the optimum trait value of the focal species at each time step to improve visualization.

## Repository structure

```text
Cosmo_Bascompte_coevo_adapt_potential/
|-- src/
|   |-- aux_functions.R                # Network reading, T-matrix and analytical GMF components
|   |-- main_functions.R               # Exact log fitness, coevolutionary model, numerical simulations
|   `-- empirical_networks/            # 230 empirical incidence matrices (Web of Life)
|-- scripts/
|   |-- packages_dependencies.R        # Installs the R packages used by the project
|   |-- run_coevolution.R              # Coevolutionary scenario for every empirical network
|   |-- run_random_scenario.R          # Scenario with independent fluctuations, one run per number of partners
|   |-- main_fig1.R                    # Main Figure 1
|   |-- main_fig2.R                    # Main Figure 2
|   |-- main_fig3.R                    # Main Figure 3
|   `-- main_fig4.R                    # Main Figure 4
|-- results/
|   |-- coevolution_gmf.csv            # Coevolutionary scenario, one row per species
|   |-- network_subset.csv             # Analysed networks; the row order sets the random seeds
|   |-- random_gmf_by_degree.csv       # Independent fluctuations, one row per number of partners
|   `-- fig4_indirecteffects_fits.csv  # Linear regressions of Figure 4
`-- figures/                           # Manuscript figures and the PNG icons used by Figures 1 and 2
```

The table `results/coevolution_gmf.csv` stores the quantities of equation 5 of the main text under the following column names:

| Column | Quantity |
|---|---|
| `gmf_exact`, `ln_gmf_exact` | GMF from numerical simulations of the exact coevolutionary model |
| `gmf_approx`, `ln_gmf_approx` | Analytical GMF (equation 5 of the main text) |
| `var_selective_target` | Conflict of evolutionary effects, $`\sigma_\theta^{2}(S_i - m_i\,\bar{S}_j)`$ |
| `generalist_benefit` | Portfolio effect of generalists under coevolution |
| `R_i`, `R_bar` | Diversity of coevolutionary pathways, $`S_i`$ and $`\bar{S}_j`$ |
| `nu_bar` | Average temporal variance of the trait differences, $`\bar{\nu}_i`$ |
| `r2_partner` | Correlation among partners, $`r_i^{2}`$ |
| `indirect_network` | Contribution of indirect effects, $`F`$ |
| `leading_eigenvalue` | Leading eigenvalue of the adjacency matrix, $`\lambda`$ |

## Running the analyses

The analyses are written in R (version 4.1.0 or later). The figure scripts need ggplot2 3.5.0, patchwork 1.2.0 and dplyr 1.1.0, or later versions. All scripts locate files with the `here` package, so run them from the project root or from the `coevo_gmf.Rproj` project. To install the required packages, run:

```bash
Rscript scripts/packages_dependencies.R
```

The repository already contains the result tables needed to reproduce the manuscript figures. To generate the figures, run the following scripts:

```r
source("scripts/main_fig1.R")
source("scripts/main_fig2.R")
source("scripts/main_fig3.R")
source("scripts/main_fig4.R")
```

Figure 1 is illustrative and needs no result table. Figure 2 reads the empirical networks to obtain the number of partners of each species, and runs its own numerical simulations. Figures 3 and 4 read `results/coevolution_gmf.csv`, and Figure 3 also reads `results/network_subset.csv`.

To rebuild the result tables from the empirical networks, run the scripts of the two scenarios:

```r
source("scripts/run_coevolution.R")
source("scripts/run_random_scenario.R")
```

The coevolutionary scenario simulates 5,000 time steps of environmental change for each of the 220 networks and takes hours, even in parallel. Each network receives the seed `42 + net_id`, where `net_id` is the row index of the network in `results/network_subset.csv`. Figure 3 simulates again the coevolution of its two focal species from these seeds, and stops if the resulting GMF differs from the stored one.