# Population Structure and Correlations Among Loci

*INDIVIDUALS RARELY MATE COMPLETELY AT RANDOM*; your parents weren't two Bilateria plucked at random from the tree of life. Even within species, there's often geographically-restricted mating among individuals. Individuals tend to mate with individuals from the same, or closely related sets of populations. This form of non-random mating is called population structure and can have profound effects on the distribution of genetic variation within and among natural populations.

Populations can often differ in their allele frequencies, either due to genetic drift or selection driving differentiation among populations. In this chapter we'll talk through some ways to summarize and visualize population genetic structure. Population differentiation is also a major driver of correlations in allelic state among loci, and we'll start our discussion of these correlations at the end of this chapter. One reason for talking about population structure so early in the book is that summarizing population structure is often a key initial stage in population genomic analyses. Thus you'll often encounter summaries and visualizations of population structure when we read research papers, so it's good to have some understanding of what they represent.

## Inbreeding as a summary of population structure

Our statements about inbreeding, and inbreeding coefficients, represent one natural way to summarize population structure. In the previous chapter, we defined inbreeding as having parents that are more closely related to each other than two individuals drawn at random from some reference population. The question that naturally arises is: Which reference population should we use? While I might not look inbred in comparison to allele frequencies in the United Kingdom (UK), where I am from, my parents certainly are not two individuals drawn at random from the world-wide population. If we estimated my inbreeding coefficient $F$ using allele frequencies within the UK, it would be close to zero, but would likely be larger if we used world-wide frequencies. This is because there is a somewhat lower level of expected heterozygosity within the UK than in the human population across the world as a whole.

:::{margin}
```{figure} ../figures/Pop_struct/FST_hierarchy.pdf
---
name: figure-3.1
align: left
---

\- The hierarchical nature of F-statistics. The two dots within an individual represent the two alleles at a locus for an individual $I$. We can compare the heterozygosity in individuals ($H_I$), to that found by randomly drawing alleles from the sub-population (S), to that found in the total population (T).
```
:::

Building on this idea of inbreeding coefficients estimated at various levels, {cite:t}`Wright:43` developed a set of 'F-statistics' (also called 'fixation indices') that formalize the idea of inbreeding with respect to different levels of population structure {cite:p}`Wright:43,Wright:49`. See {numref}`figure-3.1` for a schematic diagram. Wright defined $F_{XY}$ as the correlation between random gametes, drawn from the same level $X$, relative to level $Y$. We will return to why $F$-statistics are statements about correlations between alleles in just a moment. One commonly used $F$-statistic is $F_{IS}$, which is the inbreeding coefficient between an individual ($I$) and the subpopulation ($S$). Consider a single locus, where in a subpopulation ($S$) a fraction $H_I=f_{12}$ of individuals are heterozygous. In this subpopulation, let the frequency of allele $A_1$ be $p_S$, such that the expected heterozygosity under random mating is $H_S = 2 p_S (1 - p_S)$. We will write $F_{IS}$ as

:::{math}
:label: eq-3.1

F_{IS} = 1-\frac{H_I}{H_S}= 1-\frac{f_{12}}{2p_Sq_S},
:::

a direct analog of Equation {eq}`eq-2.12`. Hence, $F_{IS}$ is the relative difference between observed and expected heterozygosity due to a deviation from random mating within the subpopulation. We could also compare the observed heterozygosity in individuals ($H_I$) to that expected in the total population, $H_T$. If the frequency of allele $A_1$ in the total population is $p_T$, then we can write $F_{IT}$ as

:::{math}
:label: eq-3.2

F_{IT} =1-\frac{H_I}{H_T}= 1-\frac{f_{12}}{2p_Tq_T},
:::

which compares heterozygosity in individuals to that expected in the total population. As a simple extension of this, we could imagine comparing the expected heterozygosity in the subpopulation ($H_S$) to that expected in the total population $H_T$, via $F_{ST}$:

:::{math}
:label: eq-3.3

F_{ST} = 1-\frac{H_S}{H_T}=1-\frac{2p_Sq_S}{2p_Tq_T}.
:::

We can relate the three $F$-statistics to each other as

:::{math}
:label: eq-3.4

(1-F_{IT}) =\frac{H_I}{H_S} \frac{H_S}{H_T}=(1-F_{IS})(1-F_{ST}).
:::

Hence, the reduction in heterozygosity within individuals compared to that expected in the total population can be decomposed to the reduction in heterozygosity of individuals compared to the subpopulation, and the reduction in heterozygosity from the total population to that in the subpopulation.

If we want a summary of population structure across multiple subpopulations, we can average $H_I$ and/or $H_S$ across populations, and use a $p_T$ calculated by averaging $p_S$ across subpopulations (or our samples from sub-populations). For example, the average $\bar{F_{ST}}$ across $K$ subpopulations (sampled with equal effort) is

:::{math}
:label: eq-3.5

\bar{F_{ST}} = 1 - \frac{\bar{H}_{S}}{H_T},
:::

:::{margin}
Averaging heterozygosity across loci first, then calculating $F_{ST}$, rather than calculating $F_{ST}$ for each locus individually and then taking the average, has better statistical properties as statistical noise in the denominator is averaged out.
:::

where $\bar{H}_S = \frac{1}{K} \sum_{i = 1}^{K} H_{S}^{(i)}$, and $H_{S}^{(i)} = 2 p_{i} q_{i}$ is the expected heterozygosity in subpopulation $i$. It follows that the average heterozygosity of the sub-populations $\bar{H}_S  \leq H_T$, and so $\bar{F_{ST}} \geq 0$ and $\bar{F_{IS}} \leq \bar{F_{IT}}$. This observation that the average heterozygosity of the sub-populations must be less than of equal to that of the total population is called the Wahlund effect {cite:p}`wahlund1928zusammensetzung`. Furthermore, if we have multiple sites, we can replace $H_I$, $H_S$, and $H_T$ with their averages across loci (as above).

:::{admonition} Question 1
:name: question-3.1

In a species of lemurs, you estimate the allele frequency to be $20\%$. In a particular population, you estimate that the allele frequency is $10\%$. In this population, only $9\%$ of individuals are heterozygote. What is $F_{IT}$, $F_{ST}$, and $F_{IS}$ for this population?   
:::

As an example of comparing a genome-wide estimate of $F_{ST}$ to that at individual loci we can look at some data from blue- and golden-winged warblers (*Vermivora cyanoptera* and *V. chrysoptera* 1-2 & 5-6 in {numref}`figure-3.2`).

:::{margin}
```{figure} ../illustration_images/alleles_genotypes/blue_golden_winged_warblers/The_warblers_of_North_America_6309257188.jpg
---
name: figure-3.2
align: left
---

\- Blue-, golden-winged, and Lawrence's warblers (*Vermivora*). <span style="font-size: smaller;">The warblers of North America. Chapman, F.M. 1907. Image from the [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/page/9165714#page/101/mode/1up). Contributed by American Museum of Natural History Library. Not in copyright.</span>
```
:::

These two species are spread across eastern Northern America, with the golden-winged warbler having a smaller, more northernly range. They're quite different in terms of plumage, but have long been known to have similar songs and ecologies. The two species hybridize readily in the wild; in fact two other previously-recognized species, Brewster's and Lawrence's warbler (4 & 3 in {numref}`figure-3.2`), are actually found to just be hybrids between theses two species. The golden-winged warbler is listed as 'threatened' under the Canadian endangered species act as its habitat is under pressure from human activity and and due to increasing hybridization with the blue-winged warbler, which is moving north into its range. {cite:t}`Toews:16` investigated the population genomics of these warblers, sequencing ten golden- and ten blue-winged warblers. They found very low divergence among these species, with a genome-wide $F_{ST}=0.0045$. In {numref}`figure-3.3`, per SNP $F_{ST}$ is averaged in $2000$bp windows moving along the genome. The average is very low, but some regions of very high $F_{ST}$ stand out. Nearly all of these regions correspond to large allele frequency differences at loci in, or close, to genes known to be involved in plumage colouration differences in other birds.

```{figure} ../Journal_figs/alleles_genotypes/blue_golden_winged_warblers/GW_FST_warblers.png
---
name: figure-3.3
align: left
---

\- FST between blue- and golden-winged warbler population samples at SNPs across the genome. Each dot is a SNP, and SNPs are coloured alternating by scaffold. Thanks to David Toews for the figure.
```

To illustrate these frequency differences @Toews:16 genotyped a SNP in each of these high-$F_{ST}$ regions. Here's their genotyping counts from the SNP, segregating for an allele 1 and 2, in the *Wnt* region, a key regulatory gene involved in feather development:

:::{table}
| Species-Genotypes | 11  | 12  | 22  |
| :---------------: | :-: | :-: | :-: |
| Blue-winged       | 2   | 21  | 31  |
| Golden-winged     | 48  | 12  | 1   |
:::

:::{admonition} Question 2
:name: question-3.2

With reference to the table of *Wnt*-allele counts:

- **A)** Calculate $F_{IS}$ in blue-winged warblers.
- **B)** Calculate $F_{ST}$ for the sub-population of blue-winged warblers compared to the combined sample.
- **C)** Calculate mean $F_{ST}$ across both sub-populations.
:::

#### Interpretations of F-statistics

Let's now return to Wright's definition of the $F$-statistics as correlations between random gametes, drawn from the same level $X$, relative to level $Y$. Without loss of generality, we may think about $X$ as individuals and $S$ as the subpopulation. Rewriting $F_{IS}$ in terms of the observed homozygote frequencies ($f_{11}$, $f_{22}$) and expected homozygosities ($p_{S}^2$, $q_{S}^2$) we find

:::{math}
:label: eq-3.6

F_{IS} = \frac{2p_Sq_S - f_{12}}{2p_Sq_S} = \frac{f_{11}+f_{22} - p_S^2 - q_S^2}{2p_Sq_S},
:::

:::{margin}
To see why the numerator of Equation {eq}`eq-3.6` is the covariance of a discrete random variable see Appendix Equation {eq}`eq-A.41`, where we imagine that the random variable is $1$ if the alleles drawn from the population are the same and $0$ if not. The denominator is the binomial variance of a sample of two, and so our equation is a covariance divided by a variance and so interpretable as a correlation (see Appendix Equation {eq}`eq-A.43`).
:::

using the fact that $p^2+2pq+q^2=1$, and $f_{12} = 1 - f_{11} - f_{12}$. The form of Equation {eq}`eq-3.6` reveals that $F_{IS}$ is the covariance between pairs of alleles found in an individual, divided by the expected variance under binomial sampling. Thus, $F$-statistics can be understood as the correlation between alleles drawn from a population (or an individual) above that expected by chance (i.e. drawing alleles sampled at random from some broader population).

We can also interpret $F$-statistics as proportions of variance explained by different levels of population structure. To see this, let us think about $F_{ST}$ averaged over $K$ subpopulations, whose frequencies are $p_1,\dots,p_K$. The frequency in the total population is $p_T=\bar{p} = \frac{1}{K} \sum_{i=1}^K p_i$. Then, we can write

:::{math}
:label: eq-3.7

F_{ST} &= \frac{2 \bar{p}\bar{q} - \frac{1}{K}\sum_{i=1}^K 2p_iq_i }{2 \bar{p}\bar{q}} = \frac{ \left(\frac{1}{K} \sum_{i=1}^K p_i^2 + \frac{1}{K} \sum_{i=1}^K q_i^2 \right) -  \bar{p}^2-\bar{q}^2 }{2\bar{p}\bar{q}}\\
    &= \frac{\mathrm{Var}(p_1,\dots,p_K)}{\mathrm{Var}(\bar{p})},
:::

:::{margin}
This follows because the numerator, in the middle step of Equation {eq}`eq-3.7`, is the averaged squared frequency minus the squared frequency, i.e. the variance (see Appendix Equation {eq}`eq-A.23`).
:::

which shows that $F_{ST}$ is the proportion of the variance explained by the subpopulation labels.

## Other approaches to population structure

There is a broad spectrum of methods to describe patterns of population structure in population genetic datasets. We'll briefly discuss two broad-classes of methods that appear often in the literature: assignment methods and principal components analysis.

## Assignment Methods

Here we'll describe a simple probabilistic assignment to find the probability that an individual of unknown population comes from one of $K$ predefined populations. For example, there are three broad populations of common chimpanzee (*Pan troglodytes*) in Africa: western, central, and eastern. Imagine that we have a chimpanzee whose population of origin is unknown (e.g. it's from an illegal private collection). If we have genotyped a set of unlinked markers from a panel of individuals representative of these populations, we can calculate the probability that our chimp comes from each of these populations.

We'll then briefly explain how to extend this idea to cluster a set of individuals into $K$ initially unknown populations. This method is a simplified version of what population genetics clustering algorithms such as STRUCTURE and ADMIXTURE do {cite:p}`pritchard:00,alexander:09`.

#### A simple assignment method

We have genotype data from unlinked $S$ biallelic loci for $K$ populations. The allele frequency of allele $A_1$ at locus $l$ in population $k$ is denoted by $p_{k,l}$, so that the allele frequencies in population 1 are $p_{1,1},\dots, p_{1,L}$ and population 2 are $p_{2,1},\dots, p_{2,L}$ and so on.

You genotype a new individual from an unknown population at these $L$ loci. This individual's genotype at locus $l$ is $g_l$, where $g_l$ denotes the number of copies of allele $A_1$ this individual carries at this locus ($g_l=0,1,2$).

The probability of this individual's genotype at locus $l$ conditional on coming from population $k$, i.e. their alleles being a random HW draw from population $k$, is

:::{math}
:label: eq-3.8

P(g_l | \textrm{pop k}) =
    \begin{cases}
        (1-p_{k,l})^2  & g_l=0 \\
        2 p_{k,l} (1-p_{k,l}) & g_l=1\\
        p_{k,l}^2  & g_l=2
    \end{cases}
:::

Assuming that the loci are independent, the probability of the individual's genotype across all S loci, conditional on the individual coming from population $k$, is

:::{math}
:label: eq-3.9

P(\textrm{ind.} | \textrm{pop k})  = \prod_{l=1}^S P(g_l | \textrm{pop k})
:::

We wish to know the probability that this new individual comes from population $k$, i.e. $P(\textrm{pop k} | \textrm{ind.})$. We can obtain this through Bayes' rule

:::{math}
:label: eq-3.10

P(\textrm{pop k} | \textrm{ind.})  = \frac{P(\textrm{ind.} | \textrm{pop k}) P(\textrm{pop k})}{P(\textrm{ind.})}
:::

where

:::{math}
:label: eq-3.11

P(\textrm{ind.}) = \sum_{k=1}^K  P(\textrm{ind.} | \textrm{pop k}) P(\textrm{pop k})
:::

:::{margin}
See the Appendix Equation {eq}`eq-A.16` for more on Bayes' Rule
:::

is the normalizing constant. We can interpret $P(\textrm{pop k})$ as the prior probability of the individual coming from population $k$, and unless we have some other prior knowledge we will assume that the new individual has a equal probability of coming from each population $P(\textrm{pop k})=\frac{1}{K}$.

We interpret

:::{math}
:label: eq-3.12

P(\textrm{pop k} | \textrm{ind.})
:::

as the posterior probability that our new individual comes from each of our $1,\dots, K$ populations.

More sophisticated versions of this are now used to allow for hybrids, e.g, we can have a proportion $q_k$ of our individual's genome come from population $k$ and estimate the set of $q_k$'s.

:::{admonition} Question 3
:name: question-3.3

Returning to our chimp example, imagine that we have genotyped a set of individuals from the Western and Eastern populations at two SNPs (we'll ignore the central population to keep things simpler). The frequency of the capital allele at two SNPs ($A/a$ and $B/b$) is given by

<center>
<div style="width: 50%;">

| Population | locus A | locus B |
| :--------: | :-----: | :-----: |
| Western    | 0.1     | 0.85    |
| Eastern    | 0.95    | 0.2     |

</div>
</center>

- **A)** Our individual, whose origin is unknown, has the genotype $AA$ at the first locus and $bb$ at the second. What is the posterior probability that our individual comes from the Western population versus Eastern chimp population?
- **B)** (Trickier) Lets assume that our individual from part A is a hybrid (not necessarily an F1). At each locus, with probability $q_W$ our individual draws an allele from the Western population and with probability $q_E=1-q_W$ they draw an allele from the Eastern population. What is the probability of our individual's genotype given $q_W$?
- **Optional** You could plot this probability as a function of $q_W$. How does your plot change if our individual is heterozygous at both loci?
:::

:::{margin}
```{figure} ../Journal_figs/alleles_genotypes/chimp/chimp.png
:name: figure-3.4
:align: left

\- Chimpanzee <span style="font-size: smaller;">Archives du Muséum d’Histoire Naturelle, Paris. (tome X, 1856) Image from the [Biodiversity Heritage Library](https://www.flickr.com/photos/biodivlibrary/19930229848/in/photolist-XGMmmV-wnaCYy-wDMHFP-dpNAc7-bw9qDC-bw9orQ-bJPuCF-eV4etD-d9s1By-cbcysm-c5izpA-bvSkgm-bu7ij8-azL7vS-ayF5zC-atDTpg-atDTuz-akH6mL-ag15Rf-ag15Lb). Contributed by Natural History Museum Library, London. Licensed under CC BY-2.0.</span>
```
:::

#### Clustering based on assignment methods

While it is great to be able to assign our individuals to a particular population, these ideas can be pushed to learn about how best to describe our genotype data in terms of discrete populations without assigning any of our individuals to populations *a priori*. We wish to cluster our individuals into $K$ unknown populations. We begin by assigning our individuals at random to these $K$ populations.

1.  Given these assignments we estimate the allele frequencies at all of our loci in each population.

2.  Given these allele frequencies we chose to reassign each individual to a population $k$ with a probability given by Equation {eq}`eq-3.9`.

We iterate steps 1 and 2 for many iterations (technically, this approach is known as *Gibbs Sampling*). If the data is sufficiently informative, the assignments and allele frequencies will quickly converge on a set of likely population assignments and allele frequencies for these populations.

```{figure} ../figures/Becquet_et_al_STRUCTURE_journal_pgen_0030066_g001.png
:name: figure-3.5
:align: left

\- {cite:t}`becquet:07` genotyped 78 common chimpanzee and 6 bonobo at over 300 polymorphic markers (in this case microsatellites). They ran STRUCTURE to cluster the individuals using these data into $K=4$ populations. In {cite:t}`becquet:07` above figure they show each individual as a vertical bar divided into four colours depicting the estimate of the fraction of ancestry that each individual draws from each of the four estimated populations (licensed under CC BY 4.0). We can see that these four colours/populations correspond to: Red, central; blue, eastern; green, western; yellow, bonobo.
```

To do this in a full Bayesian scheme we need to place priors on the allele frequencies (for example, one could use a beta distribution prior). Technically we are using the joint posterior of our allele frequencies and assignments. Programs like STRUCTURE, use this type of algorithm to cluster the individuals in an "unsupervised" manner (i.e. they work out how to assign individuals to an unknown set of populations). See {numref}`figure-3.5` for an example of {cite:t}`becquet:07` using STRUCTURE to determine the population structure of chimpanzees.

STRUCTURE-like methods have proven incredible popular and useful in examining population structure within species. However, the results of these methods are open to misinterpretation; see {cite:t}`lawson:18` for a recent discussion. Two common mistakes are 1) taking the results of STRUCTURE-like approaches for some particular value of K and taking this to represent the best way to describe population-genetic variation. 2) Thinking that these clusters represent 'pure' ancestral populations.

There is no right choice of K, the number of clusters to partition into. There are methods of judging the 'best' K by some statistical measure given some particular dataset, but that is not the same as saying this is the most meaningful level on which to summarize population structure in data. For example, running STRUCTURE on world-wide human populations for low value of K will result in population clusters that roughly align with continental populations {cite:p}`rosenberg:02`. However, that does not tell us that assigning ancestry at the level of continents is a particularly meaningful way of partitioning individuals. Running the same data for higher value of K, or within continental regions, will result in much finer-scale partitioning of continental groups {cite:p}`rosenberg:02,li:08`. No one of these layers of population structure identified is privileged as being more meaningful than another.

It is tempting to think of these clusters as representing ancestral populations, which themselves are not the result of admixture. However, that is not the case, for example, running STRUCTURE on world-wide human data identifies a cluster that contains many European individuals, however, on the basis of ancient DNA we know that modern Europeans are a mixture of distinct ancestral groups.

## Principal components analysis

Principal component analysis (PCA) is a common statistical approach to visualize high dimensional data, and used by many fields. The idea of PCA is to give a location to each individual data-point on each of a small number principal component axes. These PC axes are chosen to reflect major axes of variation in the data, with the first PC being that which explains largest variance, the second the second most, and so on. The use of PCA in population genetics was pioneered by Cavalli-Sforza and colleagues and now with large genotyping datasets, PCA has made a comeback {cite:p}`menozzi:78,patterson:06`.

Consider a dataset consisting of N individuals at $S$ biallelic SNPs. The $i^{th}$ individual's genotype data at locus $\ell$ takes a value $g_{i,\ell}=0,1,\; \text{or} \; 2$ (corresponding to the number of copies of allele $A_1$ an individual carries at this SNP). We can think of this as a $N \times S$ matrix (where usually $N \ll S$).

Denoting the sample mean allele frequency at SNP $\ell$ by $p_{\ell}$, it's common to standardize the genotype in the following way

:::{math}
:label: eq-3.13

\frac{g_{i,\ell} - 2 p_{\ell}}{\sqrt{2 p_{\ell}(1-p_{\ell})}}
:::

i.e. at each SNP we center the genotypes by subtracting the mean genotype ($2p_{\ell}$) and divide through by the square root of the expected variance assuming that alleles are sampled binomially from the mean frequency ($\sqrt{2 p_{\ell}(1-p_{\ell})}$). Doing this to all of our genotypes, we form a data matrix (of dimension $N \times S$). We can then perform principal component analysis of this data matrix to uncover the major axes of genotype variance in our sample. {numref}`figure-3.6` shows a PCA from {cite:t}`becquet:07` using the same chimpanzee data as in {numref}`figure-3.5`.