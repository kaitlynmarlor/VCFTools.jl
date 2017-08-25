"""
    gtstats(vcffile, [out])

Calculate genotype statistics for each marker in a VCF file with GT field data.
It is similar to the gsstats.jar utility  <https://faculty.washington.edu/browning/beagle_utilities/utilities.html#gtstats>.

# Input
- `vcffile`: VCF file, ending with .vcf or .vcf.gz
- `out`: IOStream for output. Defulat is the `STDOUT`
One line with 15 tab-delimited fiels is written per marker:
    - 1-8)  VCF fixed fields (CHROM, POS, ID, REF, ALT, QUAL, FILT, INFO)
    -   9)  Missing genotype count
    -  10)  Missing genotype frequency
    -  11)  Non-REF allele count
    -  12)  Non-REF allele frequency
    -  13)  Minor allele count             (REF allele vs Non-REF alleles)
    -  14)  Minor allele frequency         (REF allele vs Non-REF alleles)
    -  15)  HWE P-value                    (REF allele vs Non-REF alleles)

# Output
- `records`: number of records in the input VCF file
- `samples`: number of individuals in the input VCF file
- `lines`  : number of lines written to `out`
"""
function gtstats(vcffile::AbstractString, out::IO=STDOUT)
    # open VCF file
    vcflines = endswith(vcffile, ".gz")? countgzlines(vcffile) : countlines(vcffile)
    reader = endswith(vcffile, ".gz") ?
        VCF.Reader(GzipDecompressionStream(open(vcffile, "r"))) :
        VCF.Reader(open(vcffile, "r"))
    # set up progress bar
    records = vcflines - length(VCF.header(reader)) - 1
    out == STDOUT || (pbar = ProgressMeter.Progress(records, 1))
    # allocate ouput arrays
    samples = length(VCF.header(reader).sampleID)
    missings_by_sample = zeros(Int, samples)
    missings_by_record = Int[]
    maf_by_record = Float64[]
    minorallele_by_record = Int[]
    sizehint!(maf_by_record, records)
    sizehint!(missings_by_record, records)
    sizehint!(minorallele_by_record, records)
    # loop over records
    records = lines = 0
    for record in reader
        records += 1
        # if no "GT" field, skip this record
        VCF.findgenokey(record, "GT") == 0 && continue
        # calcuate summary statistics
        lines += 1
        n00, n01, n11, n0, n1, altfreq, reffreq, missings,
            minorallele, maf, hwepval = gtstats(record, missings_by_sample)
        missfreq = missings / (n0 + n1)
        altfreq  = n0 / (n0 + n1)
        minoralleles = minorallele == 0 ? n0 : n1
        push!(missings_by_record, missings)
        push!(maf_by_record, maf)
        push!(minorallele_by_record, minorallele)
        # output
        nbytes = record.format[1][1] - 2
        unsafe_write(out, pointer(record.data), nbytes)
        print(out, '\t', missings, '\t', missfreq, '\t', n0, '\t',
        altfreq, '\t', minoralleles, '\t', maf, '\t', hwepval, '\n')
        # update progress bar
        out == STDOUT || ProgressMeter.update!(pbar, records)
    end
    close(out); close(reader)
    return records, samples, lines, missings_by_sample, missings_by_record,
        maf_by_record, minorallele_by_record
end

"""
    gtstats(vcffile, out)

Calculate genotype statistics for each marker in a VCF file with GT field data.
Output is written to the file specified by `out`.
"""
function gtstats(vcffile::AbstractString, out::AbstractString)
    ofile = endswith(out, ".gz") ?
        GzipCompressionStream(open(out, "w")) :
        open(out, "w")
    gtstats(vcffile, ofile)
end

"""
    gtstats(record[, missings_by_sample])

Calculate genotype statistics for a VCF record with GT field.

# Input
- `record`: a VCF record
- `missings_by_sample`: accumulator of misisngs by sample, `missings_by_sample[i]`
is incremented by 1 if `i`-th individual has missing genotype in this record

# Output
- `n00`: number of homozygote ALT/ALT or ALT|ALT
- `n01`: number of heterozygote REF/ALT or REF|ALT
- `n11`: number of homozygote REF/REF or REF|REF
- `n0`: number of ALT alleles
- `n1`: number of REF alleles
- `altfreq`: proportion of ALT alleles
- `reffreq`: proportion of REF alleles
- `missings`: number of missing genotypes
- `minorallele`: minor allele, 0 (ALT allele) or 1 (REF allele)
- `maf`: minor allele frequency
- `hwepval`: Hardy-Weinberg p-value
"""
function gtstats(
    record::VCF.Record,
    missings_by_sample::Union{Vector,Void}=nothing
    )
    # n11: number of homozygote REF/REF or REF|REF
    # n00: number of homozygote ALT/ALT or ALT|ALT
    # n01: number of heterozygote REF/ALT or REF|ALT
    missings = n00 = n10 = n11 = 0
    samples = length(record.genotype)
    gtkey = VCF.findgenokey(record, "GT")
    for i in 1:samples
        geno = record.genotype[i]
        if gtkey > endof(geno)
            missings += 1
            (missings_by_sample == nothing) || (missings_by_sample[i] += 1)
        else
            # "0" => 0x30, "1" => 0x31
            if record.data[geno[gtkey][1]] == 0x30
                n00 += record.data[geno[gtkey][3]] == 0x30? 1 : 0
            elseif record.data[geno[gtkey][1]] == 0x31
                n11 += record.data[geno[gtkey][3]] == 0x31? 1 : 0
            end
        end
    end
    n01         = length(record.genotype) - missings - n00 - n11
    n0          = n01 + 2n00
    n1          = n01 + 2n11
    altfreq     = n0 / (n0 + n1)
    reffreq     = n1 / (n0 + n1)
    minorallele = n0 < n1? 0 : 1
    maf         = n0 < n1? altfreq : reffreq
    hwepval     = hwe(n00, n01, n11)
    return n00, n01, n11, n0, n1, altfreq, reffreq, missings,
        minorallele, maf, hwepval
end

"""
    hwe(n00, n01, n11)

Hardy-Weinberg equilibrium test.
"""
function hwe(n00::Integer, n01::Integer, n11::Integer)
    n = n00 + n01 + n11
    n == 0 && return 1.0
    p0 = (n01 + 2n00) / 2n
    (p0 ≤ 0.0 || p0 ≥ 1.0) && return 1.0
    p1 = 1 - p0
    # Pearson's Chi-squared test
    e00 = n * p0 * p0
    e01 = 2n * p0 * p1
    e11 = n * p1 * p1
    ts = (n00 - e00)^2 / e00 + (n01 - e01)^2 / e01 + (n11 - e11)^2 / e11
    pval = ccdf(Chi(1), ts)
    # TODO Fisher exact test
    return pval
end

function countgzlines(gzfile::AbstractString)
    endswith(gzfile, "gz") || throw(ArgumentError("File name should end with .gz"))
    countlines(GzipDecompressionStream(open(gzfile, "r")))
end
