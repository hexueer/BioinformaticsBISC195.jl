module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta,
       getKmers,
       getKmerDist

# uncomment the following line if you intend to use BioSequences types
# using BioSequences
# import BioSequences: composition, gc_content, complement, reverse_complement

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases). All invalid bases are replaced with N's.
Returns a the normalized DNA string
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    newSeq = ""
    for base in seq
        # note: `N` indicates an unknown base
        if occursin(base, "AGCTN")
            newSeq = newSeq * base
        elseif occursin(base, "RYSWKMBDHV")
            newSeq = newSeq * 'N'
        else
            error("Invalid base $base encountered")
        end
    end
    return newSeq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

"""
    basecomposition(sequence::AbstractString)

Counts the number of each type of base (A, C, G, T, N)
in a DNA sequence and returns a dictionary of those counts
"""
function composition(seq)
    counts = Dict('A' => 0, 'C' => 0, 'G' => 0, 'T' => 0, 'N' => 0)
    for base in normalizeDNA(seq)
        counts[base] += 1
    end
    return counts
end

"""
    gc_content(seq::AbstractString)

Calculates and returns the GC ratio of a DNA sequence.
The GC ratio is the total number of G and C bases divided by the total length of the sequence.
"""
# function gc_content(seq::LongDNASeq)
#     seq = composition(seq)        
#     return (seq[DNA_C] + seq[DNA_G]) / (seq[DNA_A] + seq[DNA_C] + seq[DNA_G] + seq[DNA_T] + seq[DNA_N])
# end

function gc_content(seq::AbstractString)
    counts = composition(seq)        
    return (counts['C'] + counts['G']) / (counts['A'] + counts['C'] + counts['G'] + counts['T'])
end

"""
    complement(base::Char)
    complement(seq::AbstractString)

Takes a base and returns the complementary base.
Also takes sequences and returns a sequence populated with the complementary bases.
"""
function complement(base::Char)
    comp = Dict('A'=>'T',
                'T'=>'A',
                'G'=>'C',
                'C'=>'G',
                'N'=>'N')
    return comp[base]
end

function complement(seq::AbstractString)
    return map(complement, normalizeDNA(seq))
end

"""
    reverse_complement(seq::AbstractString)

Takes a sequence and returns the reverse complement of that sequence.
"""
function reverse_complement(seq::AbstractString)
    return reverse(complement(seq))
end

"""
    function parse_fasta(path)

Reads a fasta-formated file and returns 2 vectors,
one containing the headers as 'Strings' withe leading ">" removed,
the other containing the sequences as `Strings`.

Note: function does validate DNA sequences for correctness (including Ns).

"""
function parse_fasta(path)
    headers, sequences, tempSeq = [], [], []
	for line in eachline(path)
		if startswith(line, '>')
			push!(headers, line[2:end])
			isempty(tempSeq) || push!(sequences, join(tempSeq))
			tempSeq = []
		else
			push!(tempSeq, normalizeDNA(line))
		end
	end
	push!(sequences, join(tempSeq)) 
	return headers, sequences
end

"""
    getKmers(sequence, k)

Takes a sequence and integer k and returns a set of all unique kmers of length k.

Examples
≡≡≡≡≡≡≡≡≡≡

    julia> getKmers("ggg", 3)
    Set{Any} with 1 element:
    "GGG"

    julia> getKmers("ATATATATA", 4)
    Set{Any} with 2 elements:
    "TATA"
    "ATAT"

    julia> getKmers("A", 2)
    ERROR: k must be a positive integer less than the length of the sequence
"""
function  getKmers(sequence, k)
    1 <= k <= length(sequence) || error("k must be a positive integer less than the length of the sequence")
    kmers = []
    
    stopindex = length(sequence)-k+1

    for i in 1:stopindex
        kmer = normalizeDNA(sequence[i:i+k-1]) 
        !occursin('N', kmer) && push!(kmers, kmer)
    end

    return Set(kmers)
end

"""
    getKmerDist(kmer1, kmer2)

Takes two kmer sets and calculates the distance between them.
**** Assumes that kmer sets are normalized.

The distance will always be a positive number between 0 and 1 
where identical things have a distance of 0 
and things with absolutely no similarity have a distance of 1
The distance metric is defined as 1 - (length of intersection / length of union), 
or alternatively, as ((length of set1 - set2) + (length of set2 - set1)) / length of union.
"""
function getKmerDist(kmer1, kmer2)
    return 1 - (length(intersect(kmer1, kmer2)) / length(union(kmer1, kmer2)))
end

end # module BioinformaticsBISC195