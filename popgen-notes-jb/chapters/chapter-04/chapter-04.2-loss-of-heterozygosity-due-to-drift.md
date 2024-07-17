(sec-4.1)=
# Loss of heterozygosity due to drift

Genetic drift will, in the absence of new mutations, slowly purge our population of neutral genetic diversity, as alleles slowly drift to high or low frequencies and are lost or fixed over time.

Imagine a randomly mating population of a constant size $N$ diploid individuals, and that we are examining a locus segregating for two alleles that are neutral with respect to each other. This population is randomly mating with respect to the alleles at this locus. See {numref}`figure-4.1` and {numref}`figure-4.2` to see how genetic drift proceeds,
by tracking alleles within a small population.

In generation $t$ our current level of heterozygosity is $H_t$, i.e. the probability that two randomly sampled alleles in generation $t$ are non-identical is $H_t$. Assuming that the mutation rate is zero (or vanishingly small), what is our level of heterozygosity in generation $t+1$?

```{figure} ../../figures/Loss_of_he_col_two_alleles.png
:name: figure-4.1
:align: left
\- Loss of heterozygosity over time, in the absence of new mutations. A
diploid population of 5 individuals over the generations, with lines
showing transmission. In the first generation every individual is a
heterozygote.
```

```{figure} ../../figures/Loss_of_het_2_many_alleles.png
:name: figure-4.2
:align: left
\- Loss of heterozygosity over time, in the absence of new mutations. A
diploid population of 5 individuals. In the first generation I colour
every allele a different colour so we can track their descendants.
```

In the next generation ($t+1$) we are looking at the alleles in the offspring of generation $t$. If we randomly sample two alleles in generation $t+1$ which had different parental alleles in generation $t$, that is just like drawing two random alleles from generation $t$. So the probability that these two alleles in generation $t+1$, that have different parental alleles in generation $t$, are non-identical is $H_t$.

Conversely, if the two alleles in our pair had the same parental allele in the proceeding generation (i.e. the alleles are identical by descent one generation back) then these two alleles must be identical (as we are
not allowing for any mutation).

In a diploid population of size $N$ individuals there are $2N$ alleles. The probability that our two alleles have the same parental allele in the proceeding generation is $\frac{1}{(2N)}$ and the probability that they have different parental alleles is is $1-\frac{1}{(2N)}$. So by the above argument, the expected heterozygosity in generation $t+1$ is

:::{math}
:label: eq-4.1

H_{t+1} = \frac{1}{2N} \times 0 + \left(1-\frac{1}{2N} \right)H_t
:::

Thus, if the heterozygosity in generation $0$ is $H_0$, our expected heterozygosity in generation $t$ is

:::{math}
:label: eq-4.2

H_t = \left(1-\frac{1}{2N} \right)^tH_0
:::

i.e. the expected heterozygosity within our population is decaying geometrically with each passing generation. If we assume that $\frac{1}{(2N)}\ll 1$ then we can approximate this geometric decay by an exponential decay (see [Question 4.1](#question-4.1), such that

:::{math}
:label: eq-4.3

H_t =H_{0} e^{ - \frac{t}{(2N)} }
:::

i.e. heterozygosity decays exponentially at a rate $\frac{1}{(2N)}$.

In {numref}`figure-4.3`, we show trajectories through time for 40 independently simulated loci drifting in a population of 50 individuals. Each population was started from a frequency of $30\%$. Some drift up and some drift down, eventually being lost or fixed from the population, but, on average across simulations, the allele frequency doesn't change. We also track heterozygosity, you can see that heterozygosity sometimes goes up, and sometimes goes down, but on average we are losing heterozygosity, and this rate of loss is well predicted by equation {eq}`eq-4.2`.

:::{margin}
```{figure} ../../illustration_images/Genetic_drift/smelt/20497452375_9be855d9ff_z.jpg
:name: figure-4.4
:align: left
\- Pond smelt (*Hypomesus olidus*), a close relative of delta smelt. <span style="font-size: smaller;">Bulletin of the United States Fish Commission. 1906.</span> <span style='font-size: bigger; color: red;'>MISSING</span>
```
:::

```{figure} ../../figures/WF_loss_het/WF_loss_het_N50.pdf
:name: figure-4.3
:align: left
\- Change in allele frequency and loss of heterozygosity over time for 40
replicates. Simulations of genetic drift in a diploid population of 50
individuals, in the absence of new mutations. We start 40 independent,
biallelic loci each with an initial allele at 30% frequency. The left
panel shows the allele frequency over time and the right panel shows the
heterozygosity over time, with the mean decay matching equation {eq}`eq-4.2`. Code [here](https://github.com/cooplab/popgen-notes/blob/master/Rcode/Genetic_drift/WF_loss_of_het.R).
```

:::{admonition} Question 1
:name: question-4.1

You are in charge of maintaining a population of delta smelt in the Sacramento River delta. Using a large set of microsatellites you estimate that the mean level of heterozygosity in this population is 0.005. You set yourself a goal of maintaining a level of heterozygosity of at least 0.0049 for the next two hundred years. Assuming that the smelt have a generation time of 3 years, and that only genetic drift affects these loci, what is the smallest fully outbreeding population that you would need to maintain to meet this goal?
:::

Note how this picture of decreasing heterozygosity stands in contrast to the consistency of Hardy-Weinberg equilibrium from the previous chapter. However, our Hardy-Weinberg *proportions* still hold in forming each new generation. As the offspring genotypes in the next generation ($t+1$) represent a random draw from the previous generation ($t$), if the parental frequency is $p_t$, we *expect* a proportion $2p_t(1-p_t)$ of our offspring to be heterozygotes (and HW proportions for our homozygotes). However, because population size is finite, the observed genotype frequencies in the offspring will (likely) not match exactly with our expectations. As our genotype frequencies likely change slightly due to sampling, biologically this reflects random variation in family size and Mendelian segregation, the allele frequency will changed. Therefore, while each generation represents a sample from Hardy-Weinberg proportions based on the generation before, our genotype proportions are not at an equilibrium (an unchanging state) as the underlying allele frequency changes over the generations. We'll develop some mathematical models for these allele frequency changes later on. For now, we'll simply note that under our simple model of drift (formally the Wright-Fisher model), our allele count in the $t+1^{th}$ generation represents a binomial sample (of size $2N$) from the population frequency $p_t$ in the previous generation. If you've read to here, please email Prof Coop a picture of JBS Haldane in a striped suit with the title "I'm reading the chapter 3 notes". (It's well worth googling JBS Haldane and to read more about his life; he's a true character and one of the last great polymaths.)

```{figure} ../../Journal_figs/genetic_drift/black_footed_ferrets/black_footed_ferrets_He.pdf
:name: figure-4.5
:align: left
\- Loss of heterozygosity in the Black-footed Ferrets in their declining
population. Numbers in brackets give estimated number of individuals
alive at that time. Data from {cite:t}`Wisely:02`. <span style='font-size: bigger; color: red;'>MISSING</span>
```

:::{margin}
```{figure} ../../illustration_images/Genetic_drift/Black_footed_ferrets/Black_footed_ferret.pdf
:name: figure-4.6
:align: left
\- The black-footed ferret (*M. nigripes*). <span style="font-size: smaller;">Wild animals of North America, The National geographical society, 1918.</span> <span style='font-size: bigger; color: red;'>MISSING</span>
```
:::

To see how a decline in population size can affect levels of heterozygosity, let's consider the case of black-footed ferrets (*Mustela nigripes*). The black-footed ferret population has declined dramatically through the twentieth century due to destruction of their habitat and sylvatic plague. In 1979, when the last known black-footed ferret died in captivity, they were thought to be extinct. In 1981, a very small wild population was rediscovered ($40$ individuals), but in 1985 this population suffered a number of disease outbreaks.

At that point of the $18$ remaining wild individuals were brought into captivity, 7 of which reproduced. Thanks to intense captive breeding efforts and conservation work, a wild population of over 300 individuals has been established since. However, because all of these individuals are descended from those 7 individuals who survived the bottleneck, diversity levels remain low. {cite:t}`Wisely:02` measured heterozygosity at a number of microsatellites in individuals from museum collections, showing the sharp drop in diversity as population sizes crashed (see {numref}`figure-4.5`).

:::{admonition} Question 2
:name: question-4.2

In mathematical population genetics, a commonly used approximation is $(1-x) \approx e^{-x}$ for $x << 1$ (formally, this follows from the Taylor series expansion of $\exp(-x)$, ignoring second order and higher terms of $x$, see {eq}`eq-A.4`). This approximation is especially useful for approximating a geometric decay process by an exponential decay process, e.g. $(1 - x)^t \approx e^{-xt}$. Using your calculator, or R, check how well this expression approximates the exact expression for two values of $x$, $x = 0.1$, and $0.01$, across two different values of t, $t=5$ and $t=50$. Briefly comment on your results.
:::

(sec-4.1.1)=
## Levels of diversity maintained by a balance between mutation and drift

Next we're going to consider the amount of neutral polymorphism that can be maintained in a population as a balance between genetic drift removing variation and mutation introducing new neutral variation, see {numref}`figure-4.7` for an example. Note in our example, how no single allele is maintained at a stable equilibrium, rather an equilibrium level of polymorphism is maintained by a constantly shifting set of alleles.

```{figure} ../../figures/Mut_drift_balance.png
:name: figure-4.7
:align: left
\- Mutation-drift balance. A diploid population of 5 individuals. In the
first generation everyone has the same allele (black). Each generation
the transmitted allele can mutate and we generate a new colour. In the
bottom plot, I trace the frequency of alleles in our population over
time. The mutation rate we use is very high, simply to maintain
diversity in this small population. <span style='font-size: bigger; color: red;'>MISSING</span>
```

### The neutral mutation rate

We'll first want to consider the rate at which neutral mutations arise in the population.Thinking back to our discussion of the neutral theory of molecular evolution, let's suppose that there are only two classes of mutation that can arise in our genomic region of interest: neutral mutations and highly deleterious mutations. The total mutation rate at our locus is $\mu$ per generation, i.e. per transmission from parent to child. A fraction $C$ of our mutations are new alleles that are highly deleterious and so quickly removed from the population. We'll call this $C$ parameter the constraint, and it will differ according to the genomic region we consider. The remaining fraction $(1-C)$ are our neutral mutations, such that our neutral mutation rate is $(1-C)\mu$. This is the per generation rate. In the rest of the chapter for simplicity we'll assume that $C=0$ and use a neutral mutation rate of $\mu$. However, we'll return to this discussion of constraint when we discuss molecular divergence in a subsequent chapter.

:::{admonition} Question 3
:name: question-4.3

It's worth taking a minute to get familiar with both how rare, and how common, mutation is. The per base pair mutation rate in humans is around $1.5$ $\times$ $10^{-8}$ per generation. That means, on average, we have to monitor a site for $\sim 66.6$ million transmissions from parent to child to see a mutation. Yet populations and genomes are big places, so mutations are common at these levels.

**A)** Your autosomal genome is $\sim$ 3 billion base pairs long ($3 \times 10^9$). You have two copies, the one you received from your mum and one from your dad. What is the average (i.e. the expected) number of mutations that occurred in the transmission from your mum and your dad to you?

**B)** The current human population size is $\sim$7 billion individuals. How many times, at the level of the entire human population, is a single base-pair mutated in the transmission from one generation to the next?
:::

### Levels of heterozygosity maintained as a balance between mutation and drift

Looking backwards in time from one generation to the previous generation, we are going to say that two alleles which have the same parental allele (i.e. find their common ancestor) in the preceding generation have *coalesced*, and refer to this event as a *coalescent event*. If our pairs of alleles are to be different from each other in the present day, a mutation must have occured more recently on one or other lineage before they found a common ancestor.

The probability that our pair of randomly sampled alleles have coalesced in the preceding generation is $\frac{1}{(2N)}$, and the probability that our pair of alleles fail to coalesce is $1-\frac{1}{(2N)}$.

The probability that a mutation changes the identity of the transmitted allele is $\mu$ per generation. So the probability of no mutation occurring is $(1-\mu)$. We'll assume that when a mutation occurs it creates some new allelic type which is not present in the population. This assumption (commonly called the infinitely-many-alleles model) makes the math slightly cleaner, and also is not too bad an assumption biologically. See {numref}`figure-4.7` for a depiction of mutation-drift balance in this model over the generations.

This model lets us calculate when our two alleles last shared a common ancestor and whether these alleles are identical as a result of failing to mutate since this shared ancestor. For example, we can work out the probability that our two randomly sampled alleles coalesce $2$ generations in the past (i.e. they fail to coalesce in generation $1$ and then coalesce in generation $2$), and that they are identical as

:::{math}
:label: eq-4.5
    \left(1- \frac{1}{2N} \right) \frac{1}{2N} (1-\mu)^4
:::

Note the power of $4$ is because our two alleles have to have failed to mutate through $2$ meioses each.

More generally, the probability that our alleles coalesce in generation $t+1$ (counting backwards in time) and are identical due to no mutation to either allele in the subsequent generations is

:::{math}
:label: eq-4.6
    P(\textrm{coal. in t+1 \& no mutations}) =  \frac{1}{2N} \left(1- \frac{1}{2N} \right)^t \left(1-\mu \right)^{2(t+1)}
:::

To make this slightly easier on ourselves let's further assume that $t \approx t+1$ and so rewrite this as:

:::{math}
:label: eq-4.7
    P(\textrm{coal. in t+1 \& no mutations}) \approx \frac{1}{2N} \left(1- \frac{1}{2N} \right)^t \left(1-\mu \right)^{2t}
:::

This gives us the approximate probability that two alleles will coalesce in the $(t+1)^\text{th}$ generation. In general, we may not know when two alleles may coalesce: they could coalesce in generation $t=1, t=2, \ldots$, and so on. Thus, to calculate the probability that two alleles coalesce in *any* generation before mutating, we can write:

:::{math}
:label: eq-4.8
  P(\textrm{coal. in any generation \& no mutations}) \approx & P(\textrm{coal. in} \; t=1 \; \textrm{\& no mutations}) \; + \nonumber\\
&  P(\textrm{coal. in} \; t=2 \; \textrm{\& no mutations}) + \ldots \nonumber\\
  %P(\textrm{coal. in} \; t=3 \; \textrm{\& no mutations})  +\ldots \nonumber\\
  = & \sum_{t=1}^\infty P(\textrm{coal. in } \; t \; \textrm{generations \& no mutation})
:::

an example of using the Law of Total Probability, see equation {eq}`eq-A.12`, combined with the fact that coalescing in a particular generation is mutually exclusive with coalescing in a different generation.

While we could calculate a value for this sum given $N$ and $\mu$, it's difficult to get a sense of what's going on with such a complicated expression. Here, we turn to a common approximation in population genetics (and all applied mathematics), where we assume that $\frac{1}{(2N)} \ll 1$ and $\mu \ll 1$. This allows us to approximate the geometric decay as an exponential decay (see Appendix equation {eq}`eq-A.13`). Then, the probability two alleles coalesce in generation $t+1$ and don't mutate can be written as:

:::{math}
:label: eq-4.9
    P(\textrm{coal. in t+1 \& no mutations}) &\approx \frac{1}{2N}
    \left(1- \frac{1}{2N} \right)^t \left(1-\mu \right)^{2t} \\
    & \approx \frac{1}{2N} e^{-t/(2N)} e^{-2\mu t } \\
    &=\frac{1}{2N} e^{-t(2\mu+1/(2N))}
:::

Then we can approximate the summation by an integral, giving us:

:::{math}
:label: eq-4.12
    \frac{1}{2N} \int_0^{\infty} e^{-t(2\mu+1/(2N))} dt = \frac{1/(2N)}{1/(2N)+2\mu} \label{eqn:coal_no_mut}
:::

:::{margin}
We can use a very similar argument for a haploid population and replace $\theta=4N\mu$ with $\theta=2N \mu$. Haploids can't be heterozygous, but we interpret `heterozygosity' as the probability that two alleles paired at random in our population differ from each other.
:::

The equation above gives us the probability that our two alleles coalesce at some point in time, and do not mutate before reaching their common ancestor. Equivalently, this can be thought of as the probability our two alleles coalesce *before* mutating, i.e. that they are homozygous.

Then, the complementary probability that our pair of alleles are non-identical (or heterozygous) is simply one minus this. The following equation gives the equilibrium heterozygosity in a population at equilibrium between mutation and drift:

:::{margin}
This result was derived by {cite:t}`kimura1964number` and {cite:t}`malecot:48` (see {cite:t}`malecot:69` for an English translation, the lack of earlier translation meant this result was missed). Technically we're assuming that every new mutation creates a new allele, the so-called "infinitely many alleles" model, otherwise our pair of sequences could be identical due to repeat or back mutation. See this GENETICS [blog post](http://genestogenomes.org/kimura-crow-infinite-alleles/) and {cite:t}`ewens2016motoo` for a nice discussion of the history.
:::

:::{math}
:label: eq-4.13
    H = \frac{2\mu}{1/(2N)+2\mu}  = \frac{4N\mu}{1+4N\mu} \label{eqn:hetero}
:::

The compound parameter $4N\mu$, the population-scaled mutation rate,
will come up a number of times so we'll give it its own name:


:::{math}
:label: eq-4.14
    \theta = 4N\mu
:::

:::{margin}
See Math Appendix equation {eq}`eq-A.9` for more background on conditional probabilities.
:::

What's the intuition of our equation {eq}`eq-4.13`, well the probability that any event happens in a particular generation is $P(\textrm{mutation or coalescence}) \approx \frac{1}{(2N)}+2\mu$, so conditional on an event happening the probability that it is a mutation is $P(\textrm{mutation} \mid \textrm{mutation or coalescence}) =  \frac{2\mu}{\left(\frac{1}{(2N)}+2\mu \right)}$.

So all else being equal, species with larger population sizes should
have proportionally higher levels of neutral polymorphism. Indeed,
populations of animals, e.g. birds, on small islands have lower levels
of diversity than closely related species on the mainland with larger
ranges. More generally, we do see higher levels of heterozygosity in
larger census population sizes across animals {numref}`figure-4.8`.
However, while census population sizes vary over many orders of
magnitude, levels of diversity vary much less than that. So, if levels
of diversity in natural populations represent a balance between genetic
drift and mutation, levels of genetic drift in large populations must be
a lot faster than their census population size suggests. In the next
section we'll talk about some possible reasons why.


```{figure} ../../Journal_figs/genetic_drift/Allozyme_pop_size/bird_allozyme_pop_size.pdf
---
name: figure-4.8
align: left
---
\- Average basepair heterozygosity plotted against the log of range size for endemic island and mainland bird populations {cite:p}`leroy2020endemic`. Average allozyme heterozygosity plotted against the log of census population size (N) for animals. Data from {cite:t}`soule1976allozyme` and {cite:t}`frankham1996relationship`. Code [here](https://github.com/cooplab/popgen-notes/blob/master/Journal_figs/genetic_drift/Allozyme_pop_size/allozyme_pop_size.R).
```

## The effective population size

:::{margin}
The effective population size ($N_e$) is the population size that
would result in the same rate of drift in an idealized population of constant size (following our modeling
assumptions)
as that observed in our true population.
:::

In practice, populations rarely conform to our assumptions of being
constant in size with low variance in reproductive success. Real
populations experience dramatic fluctuations in size, and there is often
high variance in reproductive success. Thus rates of drift in natural
populations are often a lot higher than the census population size would
imply. See {numref}`figure-4.9` for a depiction of a repeatedly
bottlenecked population losing diversity at a fast rate.

```{figure} ../../figures/Loss_of_he_col_alleles_varying_pop_dark.png
---
name: figure-4.9
align: left
---
\- Loss of heterozygosity over time in a bottlenecking population. A diploid population of 10 individuals, that bottlenecks down to three individuals repeatedly. In the first generation, I colour every allele a different colour so we can track their descendants. There are no new mutations. <span style='font-size: bigger; color: red;'>MISSING</span>
```

To cope with this discrepancy, population geneticists often invoke the
concept of an *effective population size* ($N_e$). In many situations
(but not all), departures from model assumptions can be captured by
substituting $N_e$ for $N$.

If population sizes vary rapidly in size, we can (if certain conditions
are met) replace our population size by the harmonic mean population
size. Consider a diploid population of variable size, whose size is
$N_t$ $t$ generations into the past. The probability our pairs of
alleles have not coalesced by generation $t$ is given by

```{margin}
To see this, note that if $1/(N_i)$ is
small, then we can approximate equation {eq}`eq-4.15` using the exponential approximation:

:::{math}
:label: eq-4.17
    \prod_{i=1}^{t} \exp \left( -\frac{1}{2N_i} \right)   =
\exp \left(- \sum_{i=1}^{t} \frac{1}{2N_i} \right) .
:::

When we put the product inside the exponent, it becomes a sum.  We can also write the probability of not coalescing by generation $t$ in a population of constant size ($N_e$) as an exponential, so that it takes the same form as the expression above on the right. Comparing the exponent in the two cases, we see

:::{math}
:label: eq-4.18
    \frac{t}{2N_e} = \sum_{i=1}^{t} \frac{1}{(2N_i)}
:::

So that if we want a constant effective population size ($N_e$) that has the same
rate of loss of heterozygosity as our variable population, we need to rearrange and solve this equation to give equation {eq}`eq-4.16`.
```

:::{math}
:label: eq-4.15
    \prod_{i=1}^{t} \left(1-\frac{1}{2N_i} \right)
:::

Note that this simply collapses to our original expression
$\left(1-\frac{1}{2N } \right)^t$ if $N_i$ is constant. Under this
model, the rate of loss of heterozygosity in this population is
equivalent to a population of effective size

:::{math}
:label: eq-4.16
    N_e =\frac{1}{\frac{1}{t} \sum_{i=1}^{t} \frac{1}{N_i} }.
:::

This is the harmonic mean of the varying population size.

Thus our effective population size, the size of an idealized constant
population which matches the rate of genetic drift, is the harmonic mean
true population size over time. The harmonic mean is very strongly
affected by small values, such that if our population size is one
million $99\%$ of the time but drops to $1000$ every hundred or so
generations, $N_e$ will be much closer to $1000$ than a million.

```{figure} ../../figures/Loss_of_he_col_alleles_varying_RS.png
---
name: figure-4.10
align: left
---
\- High variance on reproductive success increases the rate of genetic drift. A diploid population of 10 individuals, where the circled individuals have much higher reproductive success. In the first generation I colour every allele a different colour so we can track their descendants, there are no new mutations. <span style='font-size: bigger; color: red;'>MISSING</span>
```

Variance in reproductive success will also affect our effective
population size. Even if our population has a large constant size $N$
individuals, if only small proportion of them get to reproduce, then the
rate of drift will reflect this much smaller number of reproducing
individuals. See {numref}`figure-4.10` for a depiction of the higher rate
of drift in a population where there is high variance in reproductive
success.

To see one example of this, consider the case where $N_F$ of females get
to reproduce and $N_M$ males get reproduce. While every individual has a
biological mother and father, not every individual gets to be a parent.
In practice, in many animal species far more females get to reproduce
than males, i.e. $N_M <N_F$, as a few males get many mating
opportunities and many males get no/few mating opportunities (see
{cite:t}`janicke:16` for a broad analysis, and note that there a certainly many
exceptions to this general pattern). When our two alleles pick an
ancestor, $25\%$ of the time our alleles were both in a female ancestor,
in which case they are IBD with probability $1/(2N_F)$, and $25\%$ of
the time they are both in a male ancestor, in which case they coalesce
with probability $1/(2N_M)$. The remaining $50\%$ of the time, our
alleles trace back to two individuals of different sexes in the prior
generation and so cannot coalesce. Therefore, our probability of
coalescence in the preceding generation is

:::{margin}
```{figure} ../../illustration_images/Genetic_drift/Hamadryas_baboon/Hamadryas_baboon.pdf
:name: figure-4.11
:align: left
\- Male Hamadryas baboons. Up to ten females live in a harem with a single male. <span style="font-size=smaller;">Brehm's Tierleben (Brehm's animal life). Brehm,
  A.E. 1893.</span> <span style='font-size: bigger; color: red;'>MISSING</span>
```
:::

:::{math}
:label: eq-4.19
    \frac{1}{4}\left(\frac{1}{2N_M} \right)+\frac{1}{4}\left(\frac{1}{2N_F} \right)
:::

i.e. the rate of coalescence is the
harmonic mean of the two sexes' population sizes, equating this to
$\frac{1}{2N_e}$ we find

:::{math}
:label: eq-4.20
    N_e = \frac{4N_FN_M}{N_F+N_M}
:::

Thus if reproductive success is very skewed in one sex (e.g. $N_M \ll N/2$), our autosomal effective population size will be much reduced as a result. For more on how different evolutionary forces affect the rate of genetic drift, and their impact on the effective population size, see {cite:t}`charlesworth:09`.

:::{admonition} Question 4
:name: question-4.4
You are studying a population of 500 male and 500 female Hamadryas baboons. Assume that all of the females but only 1/10 of the males get to mate. What is the effective population size for the autosome?
:::

Variance in male and female reproductive success can have very different effects on chromosomes with differing modes of inheritance such as the X chromosome, mitochondria, and Y chromosome. The mitochondria (mtDNA) and Y chromosome are haploid and only inherited through the females and males respectively, so they have a haploid effective population sizes of $N_M$ and $N_F$.

```{figure} ../../Journal_figs/genetic_drift/Scythian_horses/Scythian_horses.pdf
:name: figure-4.12
:align: left

\- Levels of genome-wide diversity in Scythian horses from $2300$ year
old Scythian horses and Modern horses (Nordic). The numbers next to each
column given the fraction of diversity remaining in the present day,
Data from {cite:t}`librado2017ancient`. Code [here](https://github.com/cooplab/popgen-notes/blob/master/Journal_figs/genetic_drift/Scythian_horses/Scythian_horses.r).
```

:::{margin}
```{figure} ../../illustration_images/Genetic_drift/Scythian_rider/616px-Szkíta.jpg
:name: figure-4.13
:align: left
\- A gold plaque showing Scythian rider found in a burial mound in eastern Crimera (c400–350 BC). <span style="font-size: smaller;">Photograph: V Terebenin/State Hermitage Museum. Image from [wikimedia](https://commons.wikimedia.org/wiki/File:Szk\%C3\%ADta.jpg). This is a faithful photographic reproduction of a two-dimensional, public domain work of art.</span>
```
:::

{cite:t}`librado2017ancient` sequenced ancient DNA from 13 sacrificed stallions
from an $2300$ year old Scythian burial mound in Kazakhstan. The
Scythian were a nomadic people whose Russian Steppe empire stretched
from the Black Sea to the borders of China. They were among the first
people to master horseback warfare with both men and women riding armed
with short bows. 

By comparing these data to modern horses, {cite:t}`librado2017ancient` found that
levels of diversity had been substantially reduced on the autosomes and
greatly reduced on the Y chromosome. This contrasts with the mtDNA
where levels of diversity have decreased only slightly. This pattern
likely reflects the fact that much of modern horse breeding relies on a
breeding a small number of stallions to a large number of mares, and so
the effective population size of the Y chromosome has been much smaller
than the mtDNA leading to a much higher rate of loss of diversity on the
Y than on other chromosomes.

:::{admonition} Question 5
:name: question-4.5
Using the data on the reduction in horse genetic diversity in {numref}`figure-4.12`:

**A)** Estimate the effective number of stallions and mares contributing to the horse population using the mtDNA and Y chromosome data.

**B)** Predict what the reduction in diversity over the $2300$ years should be on the autosomes using these numbers?

Assume a horse generation time of $8$ years. Assume no new mutations
during this time interval.
:::

:::{margin}
```{figure} ../../illustration_images/Genetic_drift/splitgill_mushroom/splitgill_mushrooms.png
:name: figure-4.14
:align: left
\- Split-gill fungus (*Schizophyllum commune*). <span style="font-size=smaller;">长江三角洲及邻近地区孢子植物志 (Spore Flora of the Yangtze River Delta and Adjacent Areas) 1989. 上海自然博物馆. Image from the [Biodiversity Heritage Library](https://www.biodiversitylibrary.org/item/112321\#page/142/mode/1up). Contributed by Institute of Botany, Chinese Academy of Sciences. No known copyright restrictions.</span>
```
:::

:::{admonition} Question 6
:name: question-4.6

One of the highest levels of genetic diversity is seen in the diploid
split-gill fungus, *Schizophyllum commune*. Populations in the USA have
a sequence-level heterozygosity of $0.13$ per synonymous base
{cite:p}`baranova2015extraordinary`. {cite:t}`baranova2015extraordinary` sequenced
parents and multiple offspring to estimate that
$\mu= 2 \times 10^{-8} bp^{-1}$ per generation. What is your estimate of
the effective population size of *S. commune*?
:::