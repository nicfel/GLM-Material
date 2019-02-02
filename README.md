# Inferring time-dependent migration and coalescence patterns from genetic sequence and predictor data in structured populations

Nicola F. Müller<sup>1,2</sup>, Gytis Dudas<sup>3</sup>, Tanja Stadler<sup>1,2</sup>

<sup>1</sup>ETH Zurich, Department of Biosystems Science and Engineering, 4058 Basel, Switzerland

<sup>2</sup>Swiss Institute of Bioinformatics (SIB), Switzerland

<sup>3</sup>Vaccine and Infectious Disease Division, Fred Hutchinson Cancer Research Center, Seattle, United States



## Abstract
Population dynamics can be inferred from genetic sequence data by using phylodynamic methods. These methods typically quantify the dynamics in unstructured populations or assume migration rates and effective population sizes to be constant through time in structured populations. When considering varying rates through time in structured populations, the number of parameters to infer increases rapidly and the available data might not be sufficient to inform these. Additionally, it is often of interest to know what predicts these parameters rather than knowing the parameters themselves. Here, we introduce a method to jointly infer parameters of and predictors for time-varying migration rates and effective population sizes by using a GLM (generalized linear model) approach under the marginal approximation of the structured coalescent. Using simulations, we show that our approach is able to reliably infer the model parameters and its predictors when phylogenetic trees are known. Further, when simulating trees under the structured coalescent, we show that our new approach outperforms the discrete trait GLM model. We then apply our framework to a previously described Ebola virus dataset, where we infer the parameters and predictors from genome sequences while accounting for phylogenetic uncertainty. We infer weekly cases to be the strongest predictor for effective population size and geographic distance the strongest predictor for migration. This approach is implemented as part of the BEAST2 package MASCOT, which allows to jointly infer population dynamics, i.e. the parameters and predictors, within structured populations, the phylogenetic tree, and evolutionary parameters.

## License

The content of this project itself is licensed under the [Creative Commons Attribution 3.0 license](http://creativecommons.org/licenses/by/3.0/us/deed.en_US), and the java source code of **esco** is licensed under the [GNU General Public License](http://www.gnu.org/licenses/).
