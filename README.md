# SINGER
SINGER stands for Sampling and inference of genealogies with recombination, and it is a Bayesian method to do posterior sampling of Ancestral Recombination Graph under Sequentially Markovian Coalescent. SINGER works by iterative threading one haplotype to the partially-built ARG, until the ARG for all haplotypes have been built. After initialization, MCMC will be performed to update the ARG to explore the posterior distribution. For a full description and cite our method, use:


Here we maintained the version which is under active development, but you can still direclty download the binary files for all past versions. 

## Requirements

If you want to compile the source files, then C++17 and cmake are required. Otherwise you can also used the pre-compiled binary files on various platforms. 

## Installations

## Input and output

SINGER takes **.vcf(gz)** file and outputs a **.trees** file in tskit format. The mutations are already mapped to the branches, but non-polymorphic, multi-allelic sites and structral variants are excluded from inference. The branch length should be interpreted with units of generations, for example, for homo sapiens, you would need multiply that by 28 to convert to units of years. There will also be a **.log** file for you to check the argument you ran, and the summary statistic in MCMC iterations. 

## Documentation

To sample ARGs with SINGER, you can run command line like shown below. 

**IMPORTANT**:if you wish to get ARG for: **(1) a long chromosome or (2) a series of regions**, we have provided more support to help you (see the [next section: Tools](#Tools)). If you think there are other specific job pipeline which many people might want to use, please contact us and we might add it! 

```
path_to_singer/bin/singer -fast -Ne 1e4 -m 1.25e-8 -r 1.25e-8
-input prefix_of_vcf_file -output prefix_of_output_file
-start 0 -end 1e6
-n 1000 -thin 5
```

We specify the details of the arguments here:

|flag|required?|details|  
|-------------------|-----|---|  
|**-fast**|optional|you will run fast-SINGER with this flag, otherwise regular full SINGER|
|**-Ne**|required|the diploid effective population size, which means the haploid effective population size will be **2*Ne**|
|**-m**|required|per base pair per generation mutation rate|
|**-r**|required|per base pair per generation recombination rate|
|**-input**|required|the prefix of the input .vcf(.gz) file name|
|**-output**|required|the prefix of the output .trees file name| 
|**-n**|optional|MCMC iterations to run, default at 0, only getting initialization|
|**-thin**|optional|we sample the ARG from MCMC every this number of iterations, default at 1|
|**-seed**|optional|the seed for random number generator, default at 5498u in C++ standard|
|**-bandwith**|optional|the bandwidth of the search in fast-SINGER, default at INT_MAX|
|**-polar**|optional|the probability of incorrectly phased probability, default at 0.01|

Here are some other parameters of the software which we **DON'T** recommend changing, unless you really understand what they do in the algorithms and then feel necesary to do so:

|flag|required?|details|  
|-----|-----|--------------|  
|-penalty|optional|extra penalty for violation of infinite sites model, default at 0.01|
|-hmm_epsilon|optional|the precison parameter in branch-HMM, default at 0.01|
|-psmc_bins|optional|the number of PSMC time bins, default at 20|

## Tools

### Running SINGER for a long chromosome

Often people would like to run the ARG inference method for the entire chromosome (or even the entire genome), and we have provided a python script `parallelize_singer.py` to facilitate you to this end. It automatically handles parallelization for you and runs SINGER multi-threaded. 

```
python parallelize_singer.py -Ne 2e4 -m 1.2e-8 
```
This script will:

1. Cut the genome into windows (default at 1Mb)
2. Remove the windows of unsequenced regions (<5 variants in the window)
3. Automatically parallelize running SINGER on these windows
4. Convert the output to `.trees` files with `tskit` format


### Running SINGER for a series of regions

Sometimes it is of interest to only look at certain regions on the genome (e.g. selection signals), and we have provided support for this with the python script `multiple_windows_singer.py`. It will automatically parallelize running SINGER on the regions you specify with a given `.bed` file. 

```
python multiple_windows_singer.py
```

This script will:

1. Index the vcf file for these specified windows
2. Automatically parallelize running SINGER on these windows
3. Convert the output to `.trees` files with `tskit` format

### Convert to `tskit` format

We have our internal data structure, but these files are only useful for resume MCMC runs from there. When you finishing running, please use this script to convert it to `tskit` format for your own analysis:

```
python convert_to_tskit.py
```

## Suggestions from developer

1. As a Bayesian sampling method, SINGER works best when you sample some ARGs from posterior, **only using one single sample is NOT ideal**. To this point, we highly encourage specifying **-n, -thin** flags. You can find how we run SINGER on real datasets on:
2. To decide whether to run SINGER or fast-SINGER, the best way is to run both on test data (for example, subsampled individuals on a small genomic region), and compare the inference results from them. If fast-SINGER agrees with SINGER well then it is good to go with fast-SINGER to save computational time. 
3. It is of importance to carefully choose the parameters, such as -Ne, -m, and -r. We recommend first choosing the mutation rate m, and then based on average pairwise diversity \($\pi=4\cdot N_e \cdot m\$), you can decide the Ne parameter. If you are not super sure about the recombination rate, you can simply use the same rate as mutation rate. 
4. Unfortunately for now we only support phased, high-quality genomes, and polymorphic sites with missingness will be excluded. We are working on incorporating missingness and unphased data in the near future. ARGweaver has better support in these regards.
