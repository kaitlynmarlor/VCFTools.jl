using GeneticVariation

@testset "convert_gt" begin
    #@code_warntype convert_gt(Float64, (false, false), true, :additive)
    @inferred convert_gt(Float64, (false, false), true, :additive)
    # REF is the minor allele
    @test convert_gt(Float64, (false, false), true, :additive) == 0.0
    @test convert_gt(Float64, (false, true), true, :additive) == 1.0
    @test convert_gt(Float64, (true, false), true, :additive) == 1.0
    @test convert_gt(Float64, (true, true), true, :additive) == 2.0
    @test convert_gt(Float64, (false, false), true, :dominant) == 0.0
    @test convert_gt(Float64, (false, true), true, :dominant) == 1.0
    @test convert_gt(Float64, (true, false), true, :dominant) == 1.0
    @test convert_gt(Float64, (true, true), true, :dominant) == 1.0
    @test convert_gt(Float64, (false, false), true, :recessive) == 0.0
    @test convert_gt(Float64, (false, true), true, :recessive) == 0.0
    @test convert_gt(Float64, (true, false), true, :recessive) == 0.0
    @test convert_gt(Float64, (true, true), true, :recessive) == 1.0
    # ALT is the minor allele
    @test convert_gt(Float64, (false, false), false, :additive) == 2.0
    @test convert_gt(Float64, (false, true), false, :additive) == 1.0
    @test convert_gt(Float64, (true, false), false, :additive) == 1.0
    @test convert_gt(Float64, (true, true), false, :additive) == 0.0
    @test convert_gt(Float64, (false, false), false, :dominant) == 1.0
    @test convert_gt(Float64, (false, true), false, :dominant) == 1.0
    @test convert_gt(Float64, (true, false), false, :dominant) == 1.0
    @test convert_gt(Float64, (true, true), false, :dominant) == 0.0
    @test convert_gt(Float64, (false, false), false, :recessive) == 1.0
    @test convert_gt(Float64, (true, false), false, :recessive) == 0.0
    @test convert_gt(Float64, (false, true), false, :recessive) == 0.0
    @test convert_gt(Float64, (true, true), false, :recessive) == 0.0
end

@testset "convert_gt(vcfile)" begin
    vcffile = "test.08Jun17.d8b.vcf.gz"
    if !isfile(vcffile) 
        download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
        abspath(joinpath(dirname(pathof(VCFTools)), "..", "test/$vcffile"))) 
        run(`gunzip -k $vcffile`)
    end
    #@code_warntype convert_gt(UInt8, vcffile; impute = false, center = false, scale = false)
    @inferred convert_gt(UInt8, vcffile; impute = false, center = false, scale = false)
    # convert to a matrix of UInt8, no impute/center/scale (default)
    @time A = convert_gt(UInt8, vcffile)
    @test size(A) == (191, 1356)
    @test eltype(A) == Union{UInt8, Missing}
    # convert to a matrix of Float64, impute = center = scale = true
    @time B = convert_gt(Float64, vcffile; impute = true, center = true, scale = true)
    @test size(B) == (191, 1356)
    @test eltype(B) == Union{Float64, Missing}
    # read the first record to a nullable UInt8 vector, no impute/center/scale (default)
    reader = VCF.Reader(openvcf(vcffile, "r"))
    A1 = Vector{Union{UInt8, Missing}}(undef, nsamples(vcffile))
    @time copy_gt!(A1, reader)
    @test all(A1[findall(!ismissing, A1)] .== A[findall(!ismissing, A[:, 1])])
    @test all(findall(ismissing, A1) .== findall(ismissing, A[:, 1]))
    # convert next 5 records to a nullable Float64 matrix, impute = center = scale = true
    B5 = Matrix{Union{Float64, Missing}}(undef, nsamples(vcffile), 5)
    @time copy_gt!(B5, reader; impute = true, center = true, scale = true)
    B_copy = copy(B[:, 2:6])
    @test all(B5[findall(!ismissing, B5)] .== B_copy[findall(!ismissing, B[:, 2:6])])
    @test all(B5[findall(ismissing, B5)] .== B_copy[findall(ismissing, B_copy)])
    close(reader)
end
