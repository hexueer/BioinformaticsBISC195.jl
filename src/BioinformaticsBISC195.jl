module BioinformaticsBISC195

export normalizeDNA,
       composition,
       gc_content,
       complement,
       reverse_complement,
       parse_fasta

# uncomment the following line if you intend to use BioSequences types
# using BioSequences
# import BioSequences: composition, gc_content, complement, reverse_complement

"""
    normalizeDNA(::AbstractString)

Ensures that a sequence only contains valid bases
(or `'N'` for unknown bases).
Returns a LongDNASeq.
"""
function normalizeDNA(seq)
    seq = uppercase(string(seq))
    for base in seq
        # note: `N` indicates an unknown base
        occursin(base, "AGCTN") || error("invalid base $base")
    end
    return seq # change to `return LongDNASeq(seq)` if you want to try to use BioSequences types
end

"""
    basecomposition(sequence)

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

# """
#     gc_content(::LongDNASeq) or gc_content(seq::AbstractString)

# Calculates and returns the GC ratio of a DNA sequence.
# The GC ratio is the total number of G and C bases divided by the total length of the sequence.
# """
# function gc_content(seq::LongDNASeq)
#     seq = composition(seq)        
#     return (seq[DNA_C] + seq[DNA_G]) / (seq[DNA_A] + seq[DNA_C] + seq[DNA_G] + seq[DNA_T] + seq[DNA_N])
# end

function gc_content(seq::AbstractString)
    counts = composition(seq)        
    return (counts['C'] + counts['G']) / (counts['A'] + counts['C'] + counts['G'] + counts['T'] + counts['N'])
end

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

end # module BioinformaticsBISC195