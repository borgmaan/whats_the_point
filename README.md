Genetics of Pointing Instinct
===============

Collection of Python and R modules to investigate the genetic mechanisms driving inherited behavioral traits (like pointing) in dogs. Initial round of analyses will be conducted using Illumina Canine HD 173K SNP array.

![fst-summary](https://raw.github.com/borgmaan/whats_the_point/master/images/fst_figure_2.png)


---

## Background 

A pointing breed is a type of gundog typically used in finding game. The name pointer comes from the dog's instinct to point, by stopping and aiming its muzzle towards game. Pointers were selectively bred for dogs who had abundant pointing and backing instinct. They typically start to acquire their hunting instincts at about 2 months of age. We hope to interrogate the canine genome using [high-density SNP genotyping arrays](http://res.illumina.com/documents/products/datasheets/datasheet_caninehd.pdf) to find genetic regions under selective pressure that may be associated with pointing instinct.

---

## Analysis Approach

1. Compute pairwise genetic distances between each breed and every other breed
	* E.g. for `n` total breeds, compute `n - 1` pairwise *Fst* measures using 173k Arrays  
2. Calculate breed specific normalized Fst measures (*d(i)* from [Akey et al 2010](http://www.pnas.org/content/early/2010/01/06/0909918107)) for each pointing breed using two approaches:
	1. Using all pointer vs. non-pointer comparisons. Differentiated markers in this set presumably include regions under selection for pointing instinct and regions under selection for other breed-specific traits.
	2. Using all pointer vs. pointer comparisons. Regions associated with pointing instinct should be undifferentiated in these comparisons while regions corresponding to breed specific traits should still show differentiation.
3. To remove breed-specific signatures of selection, the *d(i)* statistics for pointer vs. pointer comparisons will be subtracted from the *d(i)* statistics resulting from the comparisons with non-pointing breeds to produce a new set of adjusted *d(i)* values.

---

## Signatures of Interest

### Selective Sweeps
The elimination of standing variation in regions linked to a recently fixed beneficial mutation is known as a “selective sweep”.

### Extended Regions with High Fst
Easier to locate. Use these to narrow search for selective sweeps. 

## TODO

* Document better
* Add Makefile to make pipeline reproducible 
* Mine NGS data to look for overlapping regions of homozygosity