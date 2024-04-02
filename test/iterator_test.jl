# Test the iteration behavior of the VCFIterator

# vcf file 
# read the data and print 
# create the VCFIterator 
# expected rows vs rows 

"""
function test_iteration()
    
    # Create a temporary VCF file for testing
    isfile("test.08Jun17.d8b.vcf.gz") || Downloads.download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
    joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))

    vcf_file = joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz")
        
    # Create the VCFIterator
    iterator = vcfiterator(vcf_file)

    # Iterate over the variants and collect the rows
    next_state = 1
    final_state = 2
    row = VCFRow("", 0, [""], "", [""], 0)
    while next_state !== final_state
        row, next_state = iterate(iterator, next_state)
    end

   #@test row == VCFRow("22", 20000086, ["rs138720731"], "T", ["C"], 100.0)
   @test string(row.CHROM) == string("22")
   @test Int(row.POS) == Int(20000086)
   @test string(row.ID) == string(["rs138720731"])
   @test string(row.REF) == string("T")
   @test string(row.ALT) == string(["C"])
   @test float(row.QUAL) == float(100.0)
   @test next_state == 2

   next_state = 1
    final_state = 3
    row = VCFRow("", 0, [""], "", [""], 0)
    while next_state !== final_state
        if next_state == 2
            row = VCFRow("", 0, [""], "", [""], 0)
        end
        row, next_state = iterate(iterator, next_state)
    end
   
    # return a tuple not an array 

    #@test row == VCFRow("22", 20000146, ["rs73387790"], "G", ["A"], 100.0)
    @test string(row.CHROM) == string("22")
    @test Int(row.POS) == Int(20000146)
    @test string(row.ID) == string(["rs73387790"])
    @test string(row.REF) == string("G")
    @test string(row.ALT) == string(["A"])
    @test float(row.QUAL) == float(100.0)
    @test next_state == 3

    # Clean up the temporary VCF file
    rm(vcf_file)
end

@testset "iterator(vcf)" begin
    test_iteration()
end 

# put genetic variant base tests here 
"""

 # Create a temporary VCF file for testing
 isfile("test.08Jun17.d8b.vcf.gz") || Downloads.download("http://faculty.washington.edu/browning/beagle/test.08Jun17.d8b.vcf.gz",
 joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz"))

 vcf_file = joinpath(dirname(pathof(VCFTools)), "..", "test/test.08Jun17.d8b.vcf.gz")

#  vcf_iter = VCFIterator(vcf_file)
#  vcf, state = iterate(vcf_iter,1299)
#  print(vcf)

 # at record 292 id is missing 
 # at record 429 qual is missing 
 # at record 1299 id is missing

#  @testset "VCF file tests" begin
#     # Initialize a counter to keep track of passed tests
    
#     # Create an iterator for your VCF file
#     vcf_iter = VCFIterator(vcf_file)

#     vcf, state = iterate(vcf_iter,1356)
#     print(vcf)
    
    # # Iterate over the length of the VCF file
    # for i in 1:nrecords(vcf_file)
    #     # Obtain the VCFRow object at the current iteration
    #     vcf_row, _ = iterate(vcf_iter, i)
            
    #     # Test chrom function
    #     @test chrom(VCFData(vcf_file), vcf_row) == vcf_row.CHROM
    #     # Test pos function
    #     @test pos(VCFData(vcf_file), vcf_row) == vcf_row.POS
    #     # Test rsid function
    #     @test rsid(VCFData(vcf_file), vcf_row) == vcf_row.ID
    #     # Test alleles function
    #     @test alleles(VCFData(vcf_file), vcf_row) == (vcf_row.REF, vcf_row.ALT)
    #     # Test alt_allele function
    #     @test alt_allele(VCFData(vcf_file), vcf_row) == vcf_row.ALT
    #     # Test ref_allele function
    #     @test ref_allele(VCFData(vcf_file), vcf_row) == vcf_row.REF
    # end

# end

# @testset "MAF Tests" begin
#     vcf_iter = VCFIterator(vcf_file)
#     num = nrecords(vcf_file)
#     vcf_data = VCFData(vcf_file)
#     records, samples, lines, missing_by_sample, missings_by_record, maf_by_record, minor_allele_by_record, hwe_by_record = gtstats(vcf_file)

#     for i in 1:nrecords(vcf_file)
#         # Obtain the VCFRow object at the current iteration
#         vcf_row, _ = iterate(vcf_iter, i)
#         @test maf(vcf_data,vcf_row) == maf_by_record[i]
#     end
# end

# @testset "HWE Tests" begin
#     vcf_iter = VCFIterator(vcf_file)
#     num = nrecords(vcf_file)
#     vcf_data = VCFData(vcf_file)
#     records, samples, lines, missing_by_sample, missings_by_record, maf_by_record, minor_allele_by_record, hwe_by_record = gtstats(vcf_file)

#     for i in 1:nrecords(vcf_file)
#         # Obtain the VCFRow object at the current iteration
#         vcf_row, _ = iterate(vcf_iter, i)
#         @test hwepval(vcf_data,vcf_row) == hwe_by_record[i]
#     end
# end

# vcf_iter = VCFIterator(vcf_file)
# vcf, state = iterate(vcf_iter,1355)
# print(vcf)

# @testset "Alternate Dosages" begin
#     vcf_iter = VCFIterator(vcf_file)
#     vcf_data = VCFData(vcf_file)

#     for i in 1:nrecords(vcf_file)
#         vcf_row, _ = iterate(vcf_iter, i)
#         ds = vcf_row.DOSAGES
#         gt = vcf_row.GENOTYPE
#     end
# end

