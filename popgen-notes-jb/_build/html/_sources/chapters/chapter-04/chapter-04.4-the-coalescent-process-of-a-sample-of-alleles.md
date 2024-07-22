# The coalescent process of a sample of alleles

Usually we are not just interested in pairs of alleles, or the average
pairwise diversity. Generally we are interested in the properties of
diversity in samples of a number of alleles drawn from the population.
Instead of just following a pair of lineages back until they coalesce,
we can follow the history of a sample of alleles back through the
population.

Consider first sampling three alleles at random from the population. The
probability that all three alleles choose exactly the same ancestral
allele one generation back is $\frac{1}{(2N)^2}$. If $N$ is
reasonably large, then this is a very small probability. As such, it is
very unlikely that our three alleles coalesce all at once, and in a
moment we'll see that it is safe to ignore such unlikely events.

```{figure} ../../figures/Coalescent/Coal_three_lineages.pdf
:name: figure-4.18
:align: left

\- A simple simulation of the coalescent process for three lineages. We
track the ancestry of three modern-day alleles, the first pair (red and
purple) coalesce four generations back, after which there are only two
independent lineages we are tracking. This pair then coalesces twelve
generations in the past. Note that different random realizations of this
process will differ from each other a lot. The $T_{MRCA}$ is $T_3+T_2$.
The total time in the tree is $T_{tot}=3T_3 + 2T_2= 25$ generations. <span style='font-size: bigger; color: red;'>MISSING</span>
```

The probability that a specific pair of alleles find a common ancestor
in the preceding generation is still $\frac{1}{(2N)}$. There are
three possible pairs of alleles, so the probability that no pair finds a
common ancestor in the preceding generation is

:::{math}
:label: eq-4.28
    \left(1-\frac{1}{2N} \right)^3 \approx \left( 1- \frac{3}{2N} \right)
:::

In making this approximation we are multiplying out the right hand-side
and ignoring terms of $1/N^2$ and higher (a Taylor approximation, see
Appendix equation {eq}`eq-A.2`). See {numref}`figure-4.18` for a random realization of
this process.

:::{margin}
said as "i choose 2"
:::

More generally, when we sample $i$ alleles there are ${i \choose 2}$
pairs, i.e. $i(i-1)/2$ pairs. Thus, the probability that no pair of
alleles in a sample of size $i$ coalesces in the preceding generation is

:::{math}
:label: eq-4.29
    \left(1-\frac{1}{(2N)} \right)^{{i \choose
 2}} \approx \left( 1- \frac{{i \choose
 2}}{2N}\right)
:::

while the probability any pair coalesces is
$\approx \frac{{i \choose
 2}}{2N}$, again using Appendix equation {eq}`eq-A.2`.

We can ignore the possibility that more than pairs of alleles (e.g.
tripletons) simultaneously coalesce at once as terms of
$\frac{1}{N^2}$ and higher can be ignored as they are vanishingly
rare. Obviously in reasonable sample sizes there are many more triples
(${i \choose 3}$) and higher order combinations than there are pairs
(${i \choose 2}$), but if $i \ll N$ then we are safe to ignore these
terms.

When there are $i$ alleles, the probability that we wait until the $t+1$
generation before any pair of alleles coalesces is

:::{math}
:label: eq-4.30
    P(T_i =t+1) = \frac{{i \choose
 2}}{2N}\left( 1- \frac{{i \choose
 2}}{2N}\right)^{t}
:::

:::{margin}
see Appendix equation {eq}`eq-A.30`.
:::

Thus the waiting time to the first
coalescent event while there are $i$ lineages is a geometrically
distributed random variable with probability of success
$p=\frac{{i \choose 2}}{2N}$, which we denote by

:::{math}
:label: eq-4.31
    T_i \sim \text{Geo}
\left(  \frac{{i \choose
      2}}{2N} \right).
:::

The mean waiting time till any of pair within
our sample coalesces is

```{margin}
To see the continuous time  version of this, note that equation {eq}`eq-4.30` is

:::{math}
:label: eq-label
    \approx  \frac{{i \choose
 2}}{2N} \exp \left( - \frac{{i \choose
 2}}{2N} t \right)
:::

The waiting time $T_i$ to the first coalescent event in a sample
of $i$ alleles is thus exponentially distributed with rate $\frac{{i \choose
 2}}{2N}$, i.e. $T_i \sim \text{Exp}\left(\frac{{i \choose
   2}}{2N} \right)$.
```

:::{math}
:label: eq-4.32
    E( T_i) = \frac{2N}{{i \choose  2}}
:::

which again
follows from the mean of a geometric random variable being
$\frac{1}{p}$.

After a pair of alleles first finds a common ancestral allele some
number of generations back in the past, we only have to keep track of
that common ancestral allele for the pair when looking further into the
past. In our example coalescent genealogy for our 3 alleles, shown in {numref}`figure-4.18`, we start by tracking the 3
lineages, then by chance the blue and purple coalesce in the four
generations back. Then we're tracking just two lineages, the red lineage
and the ancestral lineage of the blue and purple alleles; then those two
coalesce and we've found our most recent common ancestor of our sample.
Another example with four tips is shown in {numref}`figure-4.19`; we're track four lineages, then a pair
coalesce, then we tracking three lineages, then a pair coalesce, then
we're tracking two lineages, then this final pair coalesce and we've
found the most recent common ancestor of our sample (fin, end scene).

More generally, when a pair of alleles in our sample of $i$ alleles
coalesces, we then switch to having to follow $i-1$ alleles back in
time. Then when a pair of these $i-1$ alleles coalesce, we then only
have to follow $i-2$ alleles back. This process continues until we
coalesce back to a sample of two, and from there to a single most recent
common ancestor (MRCA).

### Simulating a coalescent genealogy

To simulate a coalescent genealogy at a locus for a sample of $n$
alleles we therefore simply follow the following algorithm:

1.  Set $i=n$.

2.  Simulate a random variable to be the time $T_i$ to the next
    coalescent event from $T_i \sim
      \text{Exp}\left(\frac{{i \choose
     2}}{2N} \right)$

3.  Choose a pair of alleles to coalesce at random from all possible
    pairs.

4.  Set $i=i-1$

5.  Continue looping steps 2-4 until $i=1$, i.e. the most recent common
    ancestor of the sample is found.

By following this algorithm we are generating realizations of the
genealogy of our sample.

## Expected properties of coalescent genealogies and mutations

```{figure} ../../figures/Coalescent/Coal_w_muts.pdf
:name: figure-4.19
:align: left

\- A simple coalescent tree from a single coalescent simulation, tracing
the genealogy of 4 alleles with mutational changes marked with dashes
showing transitions away from the MRCA sequence (AGTTT) . The $T_{MRCA}$
is $T_4+T_3+T_2$. The total time in the tree is
$T_{tot}=4 T_4+3T_3 + 2T_2= 54$ generations. <span style='font-size: bigger; color: red;'>MISSING</span>
```

### The expected time to the most recent common ancestor

We will first consider the time to the most recent common ancestor of
the entire sample ($T_{MRCA}$). This is 

:::{math}
:label: eq-4.34
    T_{MRCA} = \sum_{i=n}^2 T_i
:::


generations back, where we are summing from $i=n$ alleles counting
backwards to $i=2$ alleles (see {numref}`figure-4.19` for example). As our coalescent times for
different $i$ are independent, the expected time to the most recent
common ancestor is

:::{math}
:label: eq-4.35
    E(T_{MRCA}) = \sum_{i=n}^2 E(T_i) = \sum_{i=n}^2  2N/{i \choose
 2}
:::
 
Using the fact that $\frac{1}{i(i-1)}=\frac{1}{i-1} - \frac{1}{i}$
and a bit of rearrangement, we can rewrite this as

:::{math}
:label: eq-4.36
    E(T_{MRCA}) = 4N\left(1- \frac{1}{n} \right)
:::

So the average $T_{MRCA}$ scales linearly with population size $N$.
Interestingly, as we move to larger and larger samples (i.e. $n \gg 1$),
the average time to the most recent common ancestor converges on $4N$.
What's happening here is that in large samples our lineages typically
coalesce rapidly at the start and very soon coalesce down to a much
smaller number of lineages.

:::{admonition} Question 9
:name: question-4.9

Assume an autosomal effective population of 10,000 individuals (roughly
the long-term human estimate) and a generation time of 30 years. What is
the expected time to the most recent common ancestor of a sample of 20
people? What is this time for a sample of 500 people?
:::

### The expected total time in a genealogy and the number of segregating sites

Mutations fall on specific lineages of the coalescent genealogy and are
transmitted to all descendants of their lineage. Furthermore, under the
infinitely-many-sites assumption, each mutation creates a new
segregating site. The mutation process is a *Poisson process*, and the
longer a particular lineage, i.e. the more generations of meioses it
represents, the more mutations that can accumulate on it. The total
number of segregating sites in a sample is thus a function of the
*total* amount of time in the genealogy of the sample, or the sum of all
the branch lengths on the genealogical tree, $T_{tot}$. Our total amount
of time in the genealogy is

:::{math}
:label: eq-4.37
    T_{tot} = \sum_{i=n}^2 iT_i
:::

as when there are $i$ lineages, each
contributes a time $T_i$ to the total time (see {numref}`figure-4.19` for an example). Taking the expectation of
the total time in the genealogy,

:::{math}
:label: eq-4.38
    E(T_{tot}) = \sum_{i=n}^2 i \frac{2N}{{i \choose
 2} } = \sum_{i=n}^2 \frac{4N}{i -1} =\sum_{i=n-1}^1 \frac{4N}{i}
:::

we see that our expected total amount of time in the genealogy scales
linearly with our population size $N$. Our expected total amount of time
is also increasing with sample size $n$, but is doing so very slowly.
This again follows from the fact that in large samples, the initial
coalescence usually happens very rapidly, so that extra samples add
little to the total amount of time in the genealogical tree.\
We saw above that the number of mutational differences between a pair of
alleles that coalescence $T_2$ generations ago was Poisson with a mean
of $2 \mu T_2$, where $2T_{2}$ is the total branch length in this simple
2-sample genealogical tree. A mutation that occurs on any branch of our
genealogy will cause a segregating polymorphism in the sample (meeting
our infinitely-many-sites assumption). Thus, if the total time in the
genealogy is $T_{tot}$, there are $T_{tot}$ generations for mutations.
So the total number of mutations segregating in our sample ($S$) is
Poisson with mean $\mu T_{tot}$. Thus the expected number of segregating
sites in a sample of size $n$ is

:::{math}
:label: eq-4.39
    E(S) = \mu \E(T_{tot}) = \sum_{i=n-1}^1 \frac{4N\mu }{i} = \theta
\sum_{i=n-1}^1 \frac{1}{i}
:::

Note that this is
growing with the sample size $n$, albeit very slowly (roughly at the
rate of the $\log$ of the sample size). We can use this formula to
derive another estimate of the population scaled mutation rate $\theta$,
by setting our observed number of segregating sites in a sample ($S$)
equal to this expectation. We'll call this estimator
$\widehat{\theta}_W$:

:::{math}
:label: eq-4.40
    \widehat{\theta}_W =\frac{ S}{\sum_{i=n-1}^1 \frac{1}{i}}
:::

This estimator of $\theta$ was devised by {cite:t}`watterson:75`, hence the $W$.

### The neutral site-frequency spectrum

We can use our coalescent process to find the expected number of derived
alleles present $i$ times out of a sample size $n$, e.g. how many
singletons ($i = 1$) do we expect to find in our sample? For example, in
{numref}`figure-4.19` in our sample of four sequences, there are
3 singletons and 2 doubletons. The number of sites with these different
allele frequencies depends on the lengths of specific genealogical
branches. A mutation that falls on a branch with $i$ descendants will
create a derived allele with frequency $i$. For example, in our example
tree in {numref}`figure-4.19`, the total number of generations where a
mutation could arise and be a doubleton is $T_3+2T_2$, the total length
of the branch ancestral to just the orange and red allele $(T_3+T_2)$
plus the branch ancestral to just the blue and purple allele $(T_2)$.

To see how we could go about working this out, let's start by
considering the simple coalescent tree, shown in {numref}`figure-4.20`, for sample of $3$ alleles drawn from a
population. Mutations that fall on the branches coloured in black will
be derived singletons, while mutations that fall along the orange branch
will be doubletons in the sample. The total number of generations where
a singleton mutation could arise is $3 T_3 + T_2$. Note that we only
count the time where there are two lineages $(T_{2})$ once. So our
expected number of singletons, using equation {eq}`eq-4.32`, is

:::{math}
:label: eq-4.41
    E(S_i) = \mu \left( 3E(T_3) +  E(T_2) \right) = \mu \left( 3
  \frac{2N}{3}+ 2N \right) = \theta
:::


By similar logic, the time where
doubletons could arise is $T_2$ and our expected number of doubletons is
$E(S_i)
=\theta/2$. Thus, there are on average half as many doubletons as
singletons.

Extending this logic to larger samples might be doable, but is tedious
(I mean really tedious: for 10 alleles there are thousands of possible
tree shapes and the task quickly gets impossible even computationally).
A nice, relatively simple proof of the neutral site frequency spectrum
is given by {cite:t}`Hudson:15`, but we won't give this here. The general form
is:

:::{math}
:label: eq-4.42
    E(S_i) = \frac{\theta }{i}
:::

i.e.
there are twice as many singletons as doubletons, three times as many
singletons as tripletons, and so on. The other thing that will be
helpful for us to know is that neutral alleles at intermediate frequency
tend to be old, and those that are rare in the sample are on average
young. We expect to see a lot more rare alleles in our sample than
common alleles.

:::{admonition} Question 10
:name: question-4.10
There are two possible tree shapes that could relate four samples. Draw
both of them and separately colour (or otherwise mark) the branches by
where singletons, doubletons, and tripleton derived alleles could arise.
:::

We can also ask the probability of observing a derived allele
segregating at frequency $i/n$ given that the site is polymorphic in our
sample of size $n$ (i.e. given that $0<i<n$ ). We can obtain this
probability by dividing the expected number of sites segregating for an
allele at frequency $i$ by the expected number segregating at all of the
possible allele frequencies for polymorphisms in our sample

:::{math}
:label: eq-4.43
    \begin{aligned}
    P(i |0<i<n) &=\frac{E(S_i)}{\sum_{j=1}^{n-1} E(S_j)} = \frac{\frac{1}{i}}{\sum_{j=1}^{n-1} \frac{1}{j}}.\end{aligned}
:::

We can interpret this probability as the fraction of polymorphic sites
we expect to find at a frequency $i/n$.

### Tests based on the site frequency spectrum

Population geneticists have proposed a variety of ways to test whether
an observed site frequency spectrum conforms to its neutral,
constant-size expectations. These tests are useful for detecting
population size changes using data across many loci, or for detecting
the signal of selection at individual loci. One of the first tests was
proposed by {cite:t}`tajima:89`, and is called Tajima's $D$. Tajima's $D$ is

:::{math}
:label: eq-4.44
    D = \frac{\hat{\theta}_{\pi}-\hat{\theta}_{W}}{C}
:::

where the numerator is the difference between the estimate of $\theta$
based on pairwise differences and that based on segregating sites. As
these two estimators both have expectation $\theta$ under the neutral,
constant-size model, the expectation of $D$ is zero. The denominator $C$
is a positive constant; it's the square-root of an estimator of the
variance of this difference under the constant population size, neutral
model. This constant was chosen for $D$ to have mean zero and variance
$1$ under the null model, so we can test for departures from this simple
null model.

An excess of rare alleles compared to the constant-size, neutral model
will result in a negative Tajima's $D$, because each additional rare
allele increases the number of segregating sites by $1$, but only has a
small effect on the number of pairwise differences between samples. In
contrast, a positive Tajima's $D$ reflects an excess of intermediate
frequency alleles relative to the constant-size, neutral expectation.
Alleles at intermediate-frequency increase pairwise diversity more per
segregating site than typical, thus increasing $\theta_{\pi}$ more than
$\theta_{W}$. In the next section we'll see how long-term changes in
population size systematically change the site frequency spectrum and so
are detectable by statistics such as Tajima's $D$.

## Demography and the coalescent

We've already seen how changes in population size can change the rate at
which heterozygosity is lost from the population (see the discussion
around equation {eq}`eq-4.15`). If the population size in generation $i$
is $N_i$, the probability that a pair of lineages coalesce is
$\frac{1}{(2N_i)}$; this conforms to our intuition that if the
population size is small, the rate at which pairs of lineages find their
common ancestor is faster. We can potentially accommodate rapid random
fluctuations in population size by simply using the effective population
size $N_e$ in place of $N$. However, longer-term, more systematic
changes in population size will distort the coalescent genealogies, and
hence patterns of diversity, in more systematic ways.

We can see how demography potentially distorts the observed frequency
spectrum away from the neutral expectation in a very large sample of
humans shown in {numref}`figure-4.21`. For comparison, the neutral frequency
spectrum, equation {eq}`eq-4.42`, is shown as a red line. There are
vastly more rare alleles than expected under our neutral,
constant-size-size model, but the neutral prediction and reality agree
somewhat more for alleles that are more common.

Why is this? Well, these patterns are likely the result of the very
recent explosive growth in human populations. If the population has
grown rapidly, then the pairwise-coalescent rate in the past may be much
higher than the coalescent rate closer to the present. (see {numref}`figure-4.22`).

One consequence of a recent population expansion is that there is much
less genetic diversity in the population than you'd predict using the
census population size. Humans are one example of this effect; there are
$7$ billion of us alive today, but this is due to very rapid population
growth over the past thousand to tens of thousands of years. Our level
of genetic diversity is very much lower than you'd predict given our
census size, reflecting our much smaller ancestral population. A second
consequence of recent population expansion is that the deeper coalescent
branches are much more squished together in time compared to those in a
constant-sized population. Mutations on deeper branches are the source
of alleles at more intermediate frequencies, and so there are even fewer
intermediate-frequency alleles in growing populations. That's why there
are so many rare alleles, especially singletons, in this large sample of
Europeans.

Another common demographic scenario is a population bottleneck. In a
bottleneck, the population size crashes dramatically, and subsequently
recovers. For example, our population may have had size
$N_{\textrm{Big}}$ and crashed down to $N_{\textrm{Small}}$. One example
of a bottleneck is shown in {numref}`figure-4.23`. Looking at a sample of lineages drawn from the population today, if the
bottleneck was somewhat recent ($\ll N_{\textrm{Big}}$ generations in
the past) many of our lineages will not have coalesced before reaching
the bottleneck, moving backward in time. But during the bottleneck our
lineages coalesce at a much higher rate, such that many of our lineages
will coalesce if the bottleneck lasts long enough
($\sim N_{\textrm{Small}}$ generations). If the bottleneck is very
strong, then all of our lineages will coalesce during the bottleneck,
and the resulting site frequency spectrum may look very much like our
population growth model (i.e. an excess of rare alleles). However, if
some pairs of lineages escape coalescing during the bottleneck, they
will coalesce much more deeply in time (e.g. the blue and orange
ancestral lineages in {numref}`figure-4.23`).

```{figure} ../../figures/Genetic_drift/Demography/Mimulus_coalescent_times.pdf
:name: figure-4.24
:align: left
\- Diversity along a region of the Mimulus genome. Black dots give $\pi$
in 1kb windows between chromosomes sampled from two individuals, the red
line is a moving average (data from {cite:t}`brandvain:14`). Pairwise coalescent
times ($t$) estimated assuming $t= \frac{\pi}{2 \mu}$ using
$\mu_{BP}=10^{-9}$.
```

An example of this is shown {numref}`figure-4.24`, data from {cite:t}`brandvain:14`. *Mimulus
nasutus* is a selfing species that arose recently from an out-crossing
progenitor *M. guttatus*, and experienced a strong bottleneck. *M.
guttatus* has very high levels of genetic diversity ($\pi=4\%$ at
synonymous sites), but *M. nasutus* has lost much of this diversity
($\pi =1\%$). Looking along the genome, between a pair of *M. guttatus*
chromosomes, levels of diversity are fairly uniformly high.

But in comparing two *M. nasutus* chromosomes, diversity is low because
the pair of lineages generally coalesce recently. Yet in a few places we
see levels of diversity comparable to *M. guttatus*; these regions
correspond to genomic sites where our pair of lineages fail to coalesce
during the bottleneck and subsequently coalesce much more deeply in the
ancestral *M. guttatus* population.

```{figure} ../../Journal_figs/genetic_drift/Maize_bottleneck/Wright_Tajima_D.pdf
:name: figure-4.26
:align: left
\- Data for polymorphism from Maize and Teosinite: 774 loci from
{cite:t}`Wright:05`. **Left)** Genetic diversity levels in maize and and teosinte
samples at each of these loci. Note how diversity levels are lower in
maize than teosinte, i.e. most points are below the red $x=y$ line.
**Right)** The distribution of Tajima's D in maize and teosinte, see how
the maize distribution is shifted towards positive values.
```

Mutations that arise on deeper lineages will be at intermediate
frequency in our sample, and so mild bottlenecks can lead to an excess
of intermediate frequency alleles compared to the standard constant-size
model. This can skew Tajima's D (see equation {eq}`eq-4.44`) towards positive values and away from its
expectation of zero. One example of this skew is shown in {numref}`figure-4.26`. Maize (*Zea mays* subsp. *mays*) was
domesticated from its wild progenitor teosinte (*Zea mays subsp.
parviglumis*) roughly ten thousand years ago. We can see how the
bottleneck associated with domestication has resulted in a loss of
genetic diversity in maize compared to teosinte, and the polymorphism
that remains is somewhat skewed towards intermediate frequencies
resulting in more positive values of Tajima's D.

:::{admonition} Question 11
:name: question-4.11
{cite:t}`voight2005interrogating` sequenced 40 autosomal regions from 15 diploid
samples of Hausa people from Yaounde, Cameroon. The average length of
locus they sequenced for each region was $2365$bp. They found that the
average number of segregating sites per locus was $S= 11.1$ and the
average $\pi = 0.0011$ per base over the loci. Is Tajima's D positive or
negative? Is a demographic model with a bottleneck or growth more
consistent with this result?
:::