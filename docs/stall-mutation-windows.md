# Stall / mutation windows vs transcription

**v1 scope.** The pipeline v1 is wild-type versus RNase H2 (`rnh201`) rNMP mapping.
This note is not a v1 result. Hydroxyurea comparisons below are Future work / v2
(Meroni et al. 2019, PMID 30847483) and are not a claim that hydroxyurea changed
mutation load or polymerase use. v1 does not assign replicative polymerases
(Daigaku et al. 2015, PMID 25664722, is the v2 reference for that question).

Spatial overlay on sacCer3. HydEn hotspots are Zhou 2019 *rnh201*
(SRR8503261 NT, SRR8503270 RNase HII). WGS T→C/A→G candidates are
GSE234409 SRR24871232 (HU) and SRR24871236 (C0). These are different experiments.
The overlay is co-occurrence, not a within-sample stall measurement.

## Definition

A continued-stall window is a 500 bp bin that is a hotspot in **both**
HydEn libraries, contains at least one WGS T→C or A→G candidate, and
places that candidate in a gene body or a 500 bp strand-aware promoter.
HII-only bins are not called continued stall (NT=0 ratio class).
Koh ribose-seq hotspots are not in this table. `rnmp_mutations.tsv`
is a type/count summary and is not a site list.

WGS positions are 1-based VCF. Overlap uses position − 1 so a call at
VCF 500 stays in bin [0, 500).

## Counts

- HydEn NT bins: 2525
- HydEn HII bins: 2947
- Continued hotspot bins (both libraries): 2014
- Bins with a WGS transition inside a HydEn hotspot: 309
- Continued-stall windows: 145
- of which touch a protein-coding feature: 143
- within 20 kb of a chromosome end, or on mtDNA: 37
- interior nuclear continued-stall windows: 108
- WGS candidates: HU 3794 / C0 3832 (raw ratio 0.99)

## Transcription classes

A site can be `both` when it sits in one gene's body and another gene's promoter.

| Set | n | gene_body | promoter | both | intergenic |
|---|---:|---:|---:|---:|---:|
| All WGS | 7626 | 4363 | 1486 | 447 | 1330 |
| Any HydEn hotspot | 1342 | 775 | 214 | 45 | 308 |
| Continued hotspot | 846 | 505 | 113 | 29 | 199 |

## Enrichment vs nuclear length

Fold is (site fraction) / (feature fraction of nuclear bp, mito excluded).
Hotspot threshold is the pipeline's arbitrary ≥2× median. A fold near 1
means the transitions sit in hotspots about as often as bin length predicts.

- Continued hotspot, HU: 1.445
- Continued hotspot, C0: 1.431
- Any HydEn hotspot, HU: 1.303
- Any HydEn hotspot, C0: 1.286
- Gene body, all WGS: 0.876
- Promoter 500 bp, all WGS: 1.038

## Future work / v2

Hydroxyurea versus the matched untreated library is not a v1 result
(Meroni et al. 2019, Nucleic Acids Research, PMID 30847483). Counts below are
not depth-adjusted. Prior mapped-read ratio was ~1.15–1.21×. n=1 vs 1, no p.
A window `hu_over_c0` above 1 is a count ratio with a pseudocount, not an effect.
Shared coordinates are consistent with the same strain background.

## Continued-stall windows (top 20 by min HydEn fold)

The top of this list is chromosome-end and mtDNA occupancy. `end` means
the bin is within 20 kb of a chromosome end, or on mtDNA. That flag does
not remove the window from the continued-stall call.

| chrom | start | end | NT fold | HII fold | HU | C0 | tx | end | genes |
|---|---:|---:|---:|---:|---:|---:|---|---:|---|
| chrXIII | 924000 | 924500 | 88.42 | 275.66 | 0 | 2 | promoter | 1 | YMR326C|YMR326C |
| chrX | 745500 | 746000 | 82.04 | 257.95 | 7 | 2 | promoter | 1 | YJR162C|YJR162C |
| chrXV | 0 | 500 | 74.36 | 222.02 | 5 | 3 | promoter | 1 | YOL166W-A|YOL166W-A |
| chrX | 715000 | 715500 | 63.28 | 94.1 | 1 | 1 | gene_body | 0 | YJR151C|DAN4 |
| chrXII | 5500 | 6000 | 55.37 | 126.06 | 1 | 1 | gene_body | 1 | YLL066W-B|YLL066W-B |
| chrIV | 1525000 | 1525500 | 41.05 | 124.27 | 0 | 1 | promoter | 1 | YDR543C|YDR543C |
| chrV | 570000 | 570500 | 28.94 | 56.05 | 1 | 3 | promoter | 1 | YER188C-A|YER188C-A |
| chrmt | 29000 | 29500 | 22.7 | 59.85 | 5 | 5 | gene_body | 1 | Q0085|ATP6 |
| chrXIV | 783500 | 784000 | 21.86 | 49.17 | 0 | 2 | promoter | 1 | YNR077C|YNR077C |
| chrI | 0 | 500 | 19.88 | 55.15 | 1 | 2 | both | 1 | YAL068W-A|YAL068W-A;YAL069W|YAL069W |
| chrI | 13000 | 13500 | 15.42 | 32.1 | 1 | 1 | gene_body | 1 | YAL064C-A|TDA8 |
| chrIII | 316000 | 316500 | 13.11 | 31.84 | 1 | 1 | gene_body | 1 | YCR108C|YCR108C |
| chrIV | 385000 | 385500 | 13.05 | 15.92 | 2 | 2 | gene_body | 0 | YDL037C|BSC1 |
| chrX | 745000 | 745500 | 12.73 | 23.18 | 1 | 2 | mixed | 1 | YJR162C|YJR162C |
| chrI | 25000 | 25500 | 10.53 | 18.59 | 2 | 2 | gene_body | 0 | YAL063C|FLO9 |
| chrVIII | 527000 | 527500 | 10.14 | 19.23 | 6 | 5 | gene_body | 0 | YHR211W|FLO5 |
| chrXII | 11500 | 12000 | 9.47 | 22.27 | 1 | 1 | promoter | 1 | YLL065W|YLL065W |
| chrIV | 384500 | 385000 | 7.4 | 11.65 | 3 | 3 | mixed | 0 | YDL037C|BSC1;YDL039C|PRM7 |
| chrV | 6000 | 6500 | 7.26 | 8.51 | 2 | 0 | gene_body | 1 | YEL074W|YEL074W |
| chrVII | 983000 | 983500 | 7.04 | 6.99 | 4 | 4 | gene_body | 0 | YGR246C|BRF1 |

## Not licensed

- A hydroxyurea effect on rNMP-proximal mutations or stalled transcription. That question is Future work / v2 (PMID 30847483), not a result of this table.
- A gene in this table is an FTD, tau, or patient-drug locus.
- HII-only bins are failed excision. Those stay in the window table as `hii_only` and are not continued-stall calls.
- Rank order inside mtDNA, rDNA, or subtelomeres is copy-number / transcription occupancy, not a disease rank.
