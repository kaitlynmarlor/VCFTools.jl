abstract type VariantIterator end
abstract type GeneticData end 

mutable struct VCFIterator <: VariantIterator
    vcffile::AbstractString
    vcf::VCF.Reader
    sample_num::Int
end

function VCFIterator(vcffile::AbstractString; sample_num::Int=VCFTools.nsamples(VCF.Reader(openvcf(vcffile, "r"))))
    vcf = VCF.Reader(openvcf(vcffile, "r")) 
    return VCFIterator(vcffile, vcf, sample_num) #then returns an iostream 
end

mutable struct VCFData <: GeneticData # feed VCFData VCFIterator function output output of openvcf function 
    file_name::AbstractString
end
# placeholder because other genotype files need that 

mutable struct VCFRow <: GeneticVariantBase.Variant
    CHROM::String
    POS::Int64
    ID::Vector{String}
    REF::String
    ALT::Vector{String}
    QUAL::Float64
    GENOTYPE::Vector{Union{Float64, Missing}}
    DOSAGES::Union{Nothing, Vector{Union{Missing, Float64}}}
    INDEX::Int64
end

# for dosages there is a key DS 
# copy_gt! and copy_ds!
# for genotypes taking values from GT

# return VCFRow item from the iterate function 

@inline function Base.eltype(::Type{<:VariantIterator}) # vector of string 
    VCFRow
end

function Base.iterate(itr::VCFIterator, state=1)
    rows = nrecords(itr.vcffile)

    if state <= 0 || state > rows
        return nothing
    else
        reader = VCF.Reader(openvcf(itr.vcffile, "r")) 
        count = 0

        for record in reader
            count += 1
            if count == state
                chr = ""
                pos = 0
                ids = Vector{String}()
                ref = ""
                alt = Vector{String}()
                qual = 0.0
            
                chr = VCF.chrom(record)
                pos = VCF.pos(record) 
                ref = VCF.ref(record)
                alt = VCF.alt(record)
                geno = gt_key(record, impute=true)
                ds = ds_key(record, impute=true)

                # Handle ID and QUAL using try...catch
                try
                    ids = VCF.id(record)
                catch e
                    println("Missing ID at record $count")
                    ids = Vector{String}()
                end
                
                try
                    qual = VCF.qual(record)
                catch e
                    println("Missing QUAL at record $count")
                    qual = 0.0
                end

                index = count 
                vcf_row = VCFRow(chr, pos, ids, ref, alt, qual, geno, ds, index)

                close(reader)
                return(vcf_row, state+1)
            end
        end
        close(reader)
        return nothing 
    end
end
        
            # out = zeros(Union{Missing, Float64}, nsamples(itr.vcffile))
            # geno = copy_gt!(out, reader)
            # return a tuple not an array 


@inline function Base.length(itr::VCFIterator)
    return nrecords(itr.vcffile)
end

@inline function vcfiterator(vcffile::AbstractString)
    VCFIterator(vcffile)
end

function ds_key(record::VCF.Record; key::String = "DS", impute::Bool = false, center::Bool = false, scale::Bool = false)
    A = Vector{Union{Missing, Float64}}(undef, length(record.genotype))  
    dskey = VariantCallFormat.findgenokey(record, key)

    for i in 1:length(record.genotype)
        if dskey === nothing
            break
        end
        geno = record.genotype[i]
        # Missing genotype: dropped field or "." => 0x2e
        if dskey > lastindex(geno) || record.data[geno[dskey]] == [0x2e]
            A[i] = missing # (impute ? ct : missing)
        else # not missing
            A[i] = parse(Float64, String(record.data[geno[dskey]]))
        end

    end

    if center || scale || impute
        total_ds = zero(Float64)
        cnt = 0

        for i in eachindex(A)
            if A[i] !== missing
                total_ds += A[i]
                cnt += 1
            end
            ct = total_ds / cnt
            wt = ct == 0 ? 1.0 : 1.0 / √(ct * (1 - ct/2))
            # center and scale if asked
            center && !ismissing(A[i]) && (A[i] -= ct)
            scale && !ismissing(A[i]) && (A[i] *= wt)
        end

        if impute
            for i in eachindex(A)
                if A[i] === missing
                    A[i] = ct
                end
            end
        end
    end 

    return A
end

function gt_key(record::VCF.Record; model::Symbol = :additive, impute::Bool = false, center::Bool = false, scale::Bool = false)
    A = Vector{Union{Missing, Float64}}(undef, length(record.genotype)) 
    gtkey = VariantCallFormat.findgenokey(record, "GT")

    _, _, _, _, _, alt_freq, _, _, _, _, _ = gtstats(record, nothing)
    ct = 2alt_freq
    wt = alt_freq == 0 ? 1.0 : 1.0 / √(2alt_freq * (1 - alt_freq))

    for i in eachindex(record.genotype)
        geno = record.genotype[i]
        # Missing genotype: dropped field or when either haplotype contains "."
        if gtkey > lastindex(geno) || geno_ismissing(record, geno[gtkey])
            if impute
                a1, a2 = rand() ≤ alt_freq, rand() ≤ alt_freq
                A[i] = convert_gt(T, (a1, a2), model)
            else
                A[i] = missing
            end
        else # not missing
            # "0" (REF) => 0x30, "1" (ALT) => 0x31
            a1 = record.data[geno[gtkey][1]] == 0x31
            a2 = record.data[geno[gtkey][3]] == 0x31
            A[i] = convert_gt(Float64, (a1, a2), model)
        end
        # center and scale if asked
        center && !ismissing(A[i]) && (A[i] -= ct)
        scale && !ismissing(A[i]) && (A[i] *= wt)
    end

    return A
end

function GeneticVariantBase.chrom(data::VCFData, row::VCFRow)::String
    return row.CHROM
end

function GeneticVariantBase.pos(data::VCFData, row::VCFRow)::Int
    return row.POS
end

function GeneticVariantBase.rsid(data::VCFData, row::VCFRow)::Vector{String}
    return row.ID
end

function GeneticVariantBase.alleles(data::VCFData, row::VCFRow)::Tuple{String, Vector{String}}
    return (row.REF, row.ALT)
end

function GeneticVariantBase.alt_allele(data::VCFData, row::VCFRow)::Vector{String}
    return row.ALT
end

function GeneticVariantBase.ref_allele(data::VCFData, row::VCFRow)::String
    return row.REF
end

# use GeneticData as argument for genetic variant base functions 
# GeneticData, VCFRow

function GeneticVariantBase.maf(data::VCFData, row::VCFRow)

    records, samples, lines, missing_by_sample, missings_by_record, maf_by_record, minor_allele_by_record, hwe_by_record = gtstats(data.file_name)
    maf = maf_by_record[row.INDEX]
    return maf

end

#copyto! function
#copygt! copydt
#reads in genotypes 0 1 2 average divided by two allele frequency of alternate allele 
#dosages for each snp 0-2

function GeneticVariantBase.hwepval(data::VCFData, row::VCFRow)
    
    records, samples, lines, missing_by_sample, missings_by_record, maf_by_record, minor_allele_by_record, hwe_by_record = gtstats(data.file_name)
    p_value = hwe_by_record[row.INDEX]
    return p_value 

end

function GeneticVariantBase.n_samples(data::VCFData)
    return nsamples(data.file_name)
end 

function GeneticVariantBase.n_records(data::VCFData)
    return nrecords(data.file_name)
end 


#copy_gt has keyword argument impute in VCFTools.jl 
# we need to run it with that option 
# impute, scale, center keyword arguments for this function 
# impute will impute the missing values 
# scale will divide by standard dev 
# center will subtract mean
# can do similar thing for other file formats 

function GeneticVariantBase.alt_dosages!(arr::AbstractArray{T}, row::VCFRow; use_genotype::Bool = false) where T <: Real


    if use_genotype
        genotypes = row.GENOTYPE
        @assert length(arr) == length(genotypes) "Array size does not match genotype size"
        for i in 1:length(genotypes)
            if genotypes[i] === missing 
                arr[i] = convert_gt(T, (a1, a2), model)
            else
                arr[i] = genotypes[i] 
            end
        end
        return arr 
    else
        dosages = row.DOSAGES
        @assert length(arr) == length(dosages) "Array size does not match dosages size"
        for i in 1:length(dosages)
            if dosages[i] === missing
                arr[i] = convert_gt(T, (a1, a2), model)
            else
                arr[i] = dosages[i] 
            end
        end
        return arr
    end
end



function GeneticVariantBase.alt_genotypes!(arr::AbstractArray{T}, row::VCFRow) where T <: Real
    genotypes = row.GENOTYPE
    @assert length(arr) == length(genotypes) "Array size does not match genotype size"
    for i in 1:length(genotypes)
        if genotypes[i] === missing 
            arr[i] = convert_gt(T, (a1, a2), model)
        else
            arr[i] = genotypes[i] 
        end
    end
    return arr 
end