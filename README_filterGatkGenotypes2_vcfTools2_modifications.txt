Scripts Modified
1. filterGatkGenotypes.py

Added helper functions and updated genotype checks to recognize triploid reference and missing genotypes:

    Added is_ref_genotype() helper function (checks for 0/0/0, 0|0|0)

    Added is_missing_genotype() helper function (checks for ././., .|.|.)

    Updated all genotype checks to use these helpers instead of hardcoded diploid checks

2. vcfTools.py - Multiple Functions Modified
For SNP Alignment (vcf_to_fasta.py):

    get_alt_degenerate() - Extended to handle:

        Diploid non-reference heterozygotes (1/2)

        Triploid heterozygotes (0/1/2, 1/2/3, etc.)

        Returns proper IUPAC codes for all ploidy combinations

        Filters homozygous reference triploids (0/0/0)

For Genotype Filtering (filterGatkGenotypes.py):

    get_genotype() - Changed ploidy detection from string length to allele counting, properly returns ././. for filtered triploids

    get_GQ_index() - Added triploid reference genotypes (0/0/0, 0|0|0) to RGQ/GQ lookup logic

    get_percent_AD() - Completely rewritten to be ploidy-agnostic, calculates minimum allele frequency across any number of alleles

    get_total_DP() - No changes (already ploidy-agnostic)

    get_DP_index() - No changes (already ploidy-agnostic)

What This Enables

✅ SNP alignments with proper IUPAC degenerate codes for mixed-ploidy data
✅ Quality filtering that works for both diploid and triploid genotypes
✅ Conservative allele balance filtering ensuring all alleles in heterozygous calls have sufficient support
✅ Proper handling of triploid reference (0/0/0) and missing (././.) genotypes throughout the pipeline

#########################3
# Changes to make vcfToFasta-Edited2.py work:

Summary of Changes to get_variant_type()
Original Function Issues:

    Only checked for diploid reference genotypes (0/0, 0|0)

    Only checked for diploid missing genotypes (./., .|.)

    Only examined the last allele in multi-allelic genotypes

    Didn't handle triploid genotypes properly

Changes Made:
1. Added Triploid Genotype Recognition

    Added 0/0/0 and 0|0|0 to reference genotype checks

    Added ././. and .|.|. to missing genotype checks

2. Multi-Allelic Genotype Handling

    Changed from checking only the last allele to examining all alleles in the genotype

    Parse genotype using re.split(r"[/|]", genotype) to get all alleles

    Extract unique non-reference alleles for checking

3. Optimized SNP Detection Logic

    Calculate ref_len once upfront instead of repeatedly in loops

    Early exit: If ref_len != 1, immediately return 'INDEL' (can't be SNP)

    Only check ALT allele lengths after confirming REF is length 1

    Return 'SNP' only if all ALT alleles present are also length 1

4. Conservative Filtering Approach

    If any allele in a heterozygous genotype is an INDEL, the entire position is classified as 'INDEL'

    Mixed SNP/INDEL positions (e.g., ALT=C,AT) are excluded from SNP alignments

    Ensures only pure SNP positions make it into phylogenetic analyses

Final Result:

The function is now ploidy-agnostic and correctly identifies SNP vs INDEL variants for diploid, triploid, and higher-ploidy genotypes, with efficient early-exit logic that avoids unnecessary checking.
