using BioinformaticsBISC195
using Test

@testset "BioinformaticsBISC195" begin
    
@testset "Using Strings" begin
    
    @testset "normalizeDNA" begin
        @test normalizeDNA("aatgn") == "AATGN"
        @test_throws Exception normalizeDNA("ZCA")
        @test_throws Exception normalizeDNA(42)
        c = normalizeDNA('C') 
        @test c == "C"
        @test typeof(c) == String
    end # normalizeDNA

    @testset "composition" begin
        seq = rand(['A','T','G','C','N'], 20) |> join
        bc = composition(seq)
        @test bc isa Dict

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)

        bc = composition(lowercase(seq))

        @test bc['A'] == count(x-> x == 'A', seq)
        @test bc['C'] == count(x-> x == 'C', seq)
        @test bc['G'] == count(x-> x == 'G', seq)
        @test bc['T'] == count(x-> x == 'T', seq)
        @test bc['N'] == count(x-> x == 'N', seq)
    end # composition

    @testset "gc_content" begin
        @test gc_content("ANTG") == 1/3
        @test gc_content("cccggg") * 100 == 100.0
        @test gc_content("ATta") == 0.0
        @test_throws Exception gc_content("ATtq")
    end # gc_content

    @testset "complement" begin
        @test complement("ATTAN") == "TAATN"
        @test complement("gcta") == "CGAT"
        @test complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception complement("AxC")
    end # complement

    @testset "reverse_complement" begin
        @test reverse_complement("ATTAN") == "NTAAT"
        @test reverse_complement("gcta") == "TAGC"
        @test reverse_complement("nnnnnnn") == "NNNNNNN"
        @test_throws Exception reverse_complement("AxC")
    end # reverse_complement

    @testset "parse_fasta" begin
        testpath = normpath(joinpath(@__DIR__, "..", "data"))
        genomes = joinpath(testpath, "cov2_genomes.fasta")
        ex1_path = joinpath(testpath, "ex1.fasta")
        ex2_path = joinpath(testpath, "ex2.fasta")

        ex1 = parse_fasta(ex1_path)
        @test ex1 isa Tuple
        @test all(x-> x isa String, ex1[1])
        @test all(x-> x isa String, ex1[2])

        @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
        @test ex1[2] == ["AATTATAGC", "CGCCCCCCAGTCGGATT"]

        @test_throws Exception parse_fasta(ex2_path)

        cov2 = parse_fasta(genomes)
        @test length(cov2[1]) == 8
        @test length(cov2[2]) == 8
    end #parse_fasta

    @testset "getKmers" begin
        @test getKmers("ggg", 3) == Set(["GGG"])
        @test getKmers("ATATATATA", 4) == Set(["TATA", "ATAT"])
        @test getKmers("GTAGAGCTGT", 6) == Set(["GTAGAG", "TAGAGC", "AGAGCT", "GAGCTG", "AGCTGT"])
        @test getKmers("GTAGAGCTGT", 8) == Set(["GTAGAGCT", "TAGAGCTG", "AGAGCTGT"])
        @test_throws Exception getKmers("A", 2)
        @test_throws Exception getKmers("CXG", 2)
    end # getKmers
    
    @testset "getKmerDist" begin
        kmer1, kmer2, kmer3 = ["GGG", "CCAT"], ["ATAT", "GGG", "ATCG"], ["TAG"]
        @test getKmerDist(kmer1, kmer2) == 0.75
        @test getKmerDist(kmer1, kmer3) == 1
        @test getKmerDist(kmer2, kmer2) == 0
    end # getKmerDist

    @testset "getKmerCount" begin
        @test getKmers("ggg", 3) == Dict("GGG" => 1)
        @test getKmers("ATATATATA", 4) == Dict("TATA" => 3, "ATAT" => 3)
        @test getKmers("GTAGAGCTGT", 6) == Dict("GTAGAG" => 1, "TAGAGC" => 1, "AGAGCT" => 1, "GAGCTG" => 1, "AGCTGT" => 1)
        @test getKmers("GTAGAGCTGT", 8) == Dict("GTAGAGCT" => 1, "TAGAGCTG" => 1, "AGAGCTGT" => 1)
        @test_throws Exception getKmerCount("A", 2)
        # Assumes data is normalized, so no exception tests besides the one above
        # for our purposes, will be running on parse_fasta-processed data so it is guaranteed pre-normalized
    end # getKmerCount

    @testset "getHeaderAttrib" begin
        @test getHeaderAttrib("NC_045512.2 |Severe acute respiratory syndrome-related coronavirus||China|2019-12", "|", [4,5]) == ["China", "2019-12"]
        @test getHeaderAttrib("MZ413975.1 |Severe acute respiratory syndrome-related coronavirus|oronasopharynx|Bangladesh: Sylhet|2020-12-31", "|", [4,5]) == ["Bangladesh: Sylhet", "2020-12-31"]
        @test getHeaderAttrib("MZ042984.1 |Severe acute respiratory syndrome-related coronavirus|oronasopharynx|Egypt|2021-04-09", "|", [1,3,5]) == ["MZ042984.1", "oronasopharynx", "2021-04-09"]
        @test_throws Exception getHeaderAttrib("A|B|C", 6)
    end # getHeaderAttrib

    @testset "getCountryFromLocation" begin
        @test getCountryFromLocation("Argentina") == "Argentina"
        @test getCountryFromLocation("Pakistan: Gilgit Baltistan") == "Pakistan"
        @test getCountryFromLocation("Japan: Kanagawa, Sagamihara") == "Japan"
        @test getCountryFromLocation("") == ""
    end # getCountryFromLocation

end # strings

# @testset "Using BioSequences" begin
    
#     @testset "normalizeDNA" begin
#         @test normalizeDNA("aatgn") == dna"AATGN"
#         @test_throws Exception normalizeDNA("ZCA")
#         @test_throws Exception normalizeDNA(42)
#         c = normalizeDNA('C') 
#         @test c == dna"c"
#         @test c isa LongSequence
#     end #  normalizeDNA

#     @testset "gc_content" begin
#         @test gc_content(dna"ANTG") == 0.25
#         @test gc_content(dna"cccggg") * 100 == 100.0
#         @test gc_content(dna"ATta") == 0.0
#     end #  composition

#     @testset "composition" begin
#         seq = rand(['A','T','G','C','N'], 20) |> join |> LongDNASeq
#         bc = composition(seq)

#         @test bc[DNA_A] == count(==(DNA_A), collect(seq))
#         @test bc[DNA_C] == count(==(DNA_C), collect(seq))
#         @test bc[DNA_G] == count(==(DNA_G), collect(seq))
#         @test bc[DNA_T] == count(==(DNA_T), collect(seq))
#         @test bc[DNA_N] == count(==(DNA_N), collect(seq))
#     end #  gc_content

#     @testset "complement" begin
#         @test complement(dna"ATTAN") == dna"TAATN"
#         @test complement(dna"gcta") == dna"CGAT"
#         @test complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  complement

#     @testset "reverse_complement" begin
#         @test reverse_complement(dna"ATTAN") == dna"NTAAT"
#         @test reverse_complement(dna"gcta") == dna"TAGC"
#         @test reverse_complement(dna"nnnnnnn") == dna"NNNNNNN"
#     end #  reverse_complement

#     @testset "parse_fasta" begin
#         testpath = normpath(joinpath(@__DIR__, "..", "data"))
#         genomes = joinpath(testpath, "cov2_genomes.fasta")
#         ex1_path = joinpath(testpath, "ex1.fasta")
#         ex2_path = joinpath(testpath, "ex2.fasta")

#         ex1 = parse_fasta(ex1_path)
#         @test ex1 isa Tuple
#         @test all(x-> x isa String, ex1[1])
#         @test all(x-> x isa LongSequence, ex1[2])

#         @test ex1[1] == ["ex1.1 | easy", "ex1.2 | multiline"]
#         @test ex1[2] == [dna"AATTATAGC", dna"CGCCCCCCAGTCGGATT"]

#         @test_throws Exception parse_fasta(ex2_path)

#         cov2 = parse_fasta(genomes)
#         @test length(cov2[1]) == 8
#         @test length(cov2[2]) == 8
#     end # parse_fasta

# end # BioSequences

end # BioinformaticsBISC195
