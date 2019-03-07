# PTNUTS
Parallel tempered No-U-Turn Sampler on the PESTO framework. 
This repository provides the inclusion of the No-U-Turn Sampler into parallel tempering. The No-U-Turn sampler is based on 

Hoffman, M. D. and Gelman, A. (2014). The No-U-Turn Sampler:  Adaptively setting path lengths in Hamiltonian Monte Carlo. Journal of Machine Learning Research, 15(1):1593–1623.


The Parallel tempering Pesto framework is based on PESTO toolbox provided by

Stapor, P., Weindl, D., Ballnus, B., Hug, S., Loos, C., Fiedler, A., Krause, S., Hross, S., Fröhlich, F., Hasenauer, J. (2018). PESTO: Parameter EStimation TOolbox. Bioinformatics, 34(4), 705-707. doi: 10.1093/bioinformatics/btx676.


This repository features a No-U-Turn sampler, the Parallel tempered No-U-Turn Sampler and two example models (mRNA transfection and JAK/STAT signaling pathway). To perform sampling on the JAK/STAT model the AMICI toolbox is required.

Fröhlich, F., Kaltenbacher, B., Theis, F. J., & Hasenauer, J. (2017). Scalable Parameter Estimation for Genome-Scale Biochemical Reaction Networks. Plos Computational Biology, 13(1), e1005331. doi: 10.1371/journal.pcbi.1005331

Fröhlich, F., Theis, F. J., Rädler, J. O., & Hasenauer, J. (2017). Parameter estimation for dynamical systems with discrete events and logical operations. Bioinformatics, 33(7), 1049-1056. doi: 10.1093/bioinformatics/btw764


Further are tools to benchmark the sampling methods included, which are provided by 

Ballnus, B., Hug,  S., Hatz, K., Görlitz, L., Hasenauer, J.,  and Theis, F. (2017). Comprehensive benchmarking of markov chain monte carlo methods for dynamical systems. BMC Systems Biology, 63.
