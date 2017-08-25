using CodecZlib, GeneticVariation

@testset "filter_genotype" begin
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    #@code_warntype filter_genotype(record, ["GT"])
    record_out = filter_genotype(record, ["GT"])
    @test VCF.format(record_out) == ["GT"]
    @test VCF.genotype(record_out) == [["0|0"], ["1|0"]]
end

@testset "flip_gt_allele" begin
    record = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    #@code_warntype flip_gt_allele(record)
    record_out = flip_gt_allele(record)
    #@show record_out
    @test VCF.ref(record_out) == "A"
    @test VCF.alt(record_out) == ["G"]
    @test VCF.format(record_out) == ["GT"]
    @test VCF.genotype(record_out) == [["1|1"], ["0|1"]]
end

@testset "match_gt_allele" begin
    record1 = VCF.Record("20\t14370\trs6054257\tG\tA\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t0|0:48:1:51,51\t1|0:48:8:51,51")
    record2 = VCF.Record("20\t14370\trs6054257\tA\tG\t29\tPASS\tNS=3;DP=14;AF=0.5;DB;H2\tGT:GQ:DP:HQ\t1|1:48:1:51,51\t0|1:48:8:51,51")
    #@code_warntype match_gt_allele(record1, record2)
    #@inferred match_gt_allele(record1, record2)
    record1_out, record2_out = match_gt_allele(record1, record2)
    #@show record1_out
    #@show record2_out
    @test VCF.ref(record1_out) == VCF.ref(record2_out)
    @test VCF.alt(record1_out) == VCF.alt(record2_out)
    @test VCF.genotype(record1_out) == VCF.genotype(record2_out)
end

@testset "binomial_proportion_test" begin
    # large sample test
    # <http://www.itl.nist.gov/div898/software/dataplot/refman1/auxillar/binotest.htm>
    pval = VCFTools.binomial_proportion_test(32, 6, 39, 5)
    @test pval ≈ 0.55760 atol=1e-4
    # Fisher exact test
    # Lady tasting tea: <https://en.wikipedia.org/wiki/Fisher%27s_exact_test>
    pval = VCFTools.binomial_proportion_test(1, 9, 11, 3)
    @test pval ≈ 0.002759456 atol=1e-6
end

@testset "conformgt_by_id" begin
    refvcf = "chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz"
    tgtvcf = "test.08Jun17.d8b.vcf.gz"
    outvcf = "conformgt.matched"
    isfile(refvcf) || download("http://bochet.gcc.biostat.washington.edu/beagle/1000_Genomes_phase1_vcf/chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz", joinpath(Pkg.dir("VCFTools"), "test/chr22.1kg.ref.phase1_release_v3.20101123.vcf.gz"))
    isfile(tgtvcf) || download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        joinpath(Pkg.dir("VCFTools"), "test/test.08Jun17.d8b.vcf.gz"))
    #@code_warntype conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    #@test @inferred conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    @time lines = conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    @test lines == 766
    reader_ref = VCF.Reader(GzipDecompressionStream(open(join([outvcf, ".ref.vcf.gz"]), "r")))
    reader_tgt = VCF.Reader(GzipDecompressionStream(open(join([outvcf, ".tgt.vcf.gz"]), "r")))
    state_ref, state_tgt = start(reader_ref), start(reader_tgt)
    while !done(reader_ref, state_ref)
        record_ref, state_ref = next(reader_ref, state_ref)
        record_tgt, state_tgt = next(reader_tgt, state_tgt)
        @test VCF.chrom(record_ref) == VCF.chrom(record_tgt)
        @test VCF.pos(record_ref) == VCF.pos(record_tgt)
        @test VCF.id(record_ref) == VCF.id(record_tgt)
    end
    # Profile.clear()
    # @profile conformgt_by_id(refvcf, tgtvcf, outvcf, "22", 20000086:20099941, false)
    # Profile.print(format=:flat)
end