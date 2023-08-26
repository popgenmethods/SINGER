# SINGER
SINGER stands for Sampling and inference of genealogies with recombination, and it is a Bayesian method to do posterior sampling of Ancestral Recombination Graph under Sequentially Markovian Coalescent. 


Here we maintained the version which is under active development, but you can still direclty download the binary files for all past versions. 

## Requirements

If you want to compile the source files, then C++17 and cmake are required. Otherwise you can also used the pre-compiled binary files on various platforms. 

## Installations

## Input and output

SINGER takes **.vcf(gz)** file and outputs a **.trees** file in tskit format. The mutations are already mapped to the branches, but non-polymorphic, multi-allelic sites and structral variants are excluded from inference. The branch length should be interpreted with units of generations, for example, for homo sapiens, you would need multiply that by 28 to convert to units of years. There will also be a **.log** file for you to check the argument you ran, and the summary statistic in MCMC iterations. 

## Documentation

To sample ARGs with SINGER, you can run command line like this:

```
path_to_singer/bin/singer -fast -Ne 1e4 -m 1.25e-8 -r 1.25e-8
-input prefix_of_vcf_file -output prefix_of_output_file
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

Here are some other parameters of the software which we **DON'T** recommend changing, unless you really understand what they do in the algorithms and then feel necesary to do so:

|flag|required?|details|  
|-----|-----|--------------|  
|-penalty|optional|extra penaly for violation of infinite sites model, default at 0.01|
|-polar_penalty|optional|the proportion of incorrectly polarized variants, default at 0.01|  
|-hmm_epsilon|optional|the precison parameter in branch-HMM, default at 0.01|
|-psmc_bins|optional|the number of PSMC time bins, default at 20|

## Suggestions from developer

1. As a Bayesian sampling method, SINGER works best when you sample some ARGs from posterior, **only using one single sample is NOT ideal**. To this point, we highly encourage specifying **-n, -thin** flags. You can find how we run SINGER on real datasets on:
2. To decide whether to run SINGER or fast-SINGER, the best way is to run both on test data (for example, subsampled individuals on a small genomic region), and compare the inference results from them. If fast-SINGER agrees with SINGER well then it is good to go with fast-SINGER to save computational time. 
3. It is of importance to carefully choose the parameters, such as -Ne and -r, you can get an estimate of Ne simply by average pairwise differences. If you are not so sure about the recombination rate, use a value from the lower side and let the algorithm tell you more.
4. Unfortunately for now we only support phased, high-quality genomes, and polymorphic sites with missingness will be excluded. We are working on incorporating missingness and unphased data in the future. ARGweaver has better support in these regards.
