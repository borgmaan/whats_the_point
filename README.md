Genetics of Pointing Instinct
===============

Collection of Python and R modules to investigate the genetic mechanisms driving inherited behavioral traits (like pointing) in dogs. Initial round of analyses will be conducted using Illumina Canine HD 173K SNP array.


Background 
---------------
A pointing breed is a type of gundog typically used in finding game. Gundogs are traditionally divided into three classes: retrievers, flushing dogs, and pointing breeds. The name pointer comes from the dog's instinct to point, by stopping and aiming its muzzle towards game. This demonstrates to the hunter the location of his or her quarry and allows them to move into gun range. Pointers were selectively bred for dogs who had abundant pointing and backing instinct. They typically start to acquire their hunting instincts at about 2 months of age.

Pointing breeds are thought to have originated in about the 1650s.

Available Data
---------------
SNP Data for Pointing Families:

Strong pointers:
Gordon Setter
Irish Setter
Pointer (English)
English Setter
Irish Red and White Setter

Versatile Breeds:
---------------
Brittany
Vizsla
German Shorthaired Pointer
Bracco Italiano
Braque Du Bourbonnais
Braque Francais
German Longhaired Pointer
German Wirehaired Pointer
Large Munsterlander
Portuguese Pointer
Pudelpointer
Spinone Italiano
Weimaraner
Wirehaired Pointing Griffon
Braque d'Auvergne
Small Munsterlander

Non-Pointing Breeds
42 Total

Steps Moving Forward:
1. Compute pairwise Fst for all breed pairs.
2. Calculate breed specific normalized Fst measures (d(i) from Akey et al 2010) for each pointing breed using two approaches. 
	a. Using all pointer vs. non-pointer comparisons. Differentiated markers in this set presumably include regions under selection for pointing instinct and regions under selection for other breed-specific traits.
	b. Using all pointer vs. pointer comparisons. Regions associated with pointing instinct should be undifferentiated in these comparisons while regions corresponding to breed specific traits should still show differentiation.
3. To remove breed-specific signatures of selection, the di statistics for pointer-vs.-pointer comparisons will be subtracted from the di statistics resulting from the comparisons with non-pointing breeds to produce a new set of adjusted di values.



Signatures of Interest

Selective Sweeps
The elimination of standing variation in regions linked to a recently fixed beneficial mutation is known as a “selective sweep”.

Difficulties:
First, the effects of selection are confounded by the effects of demographic factors. The neutral null hypothesis (absence of selective sweeps) is a composite hypothesis that also makes as- sumptions regarding the demography of the populations inves- tigated. Typically, it is assumed that the population is in equilib- rium at constant size and with no population subdivision or gene-flow with other populations.

The second challenge faced in genomic scans of selective sweeps is that much of the available data consist of SNP genotypes that had been initially identified using an ascertainment (or SNP discovery) process. Because these data deviate from a random sample of fully identified genotypes, standard population genetic methods cannot be applied without taking this “ascertainment bias” into account.


Extended Regions with High Fst

