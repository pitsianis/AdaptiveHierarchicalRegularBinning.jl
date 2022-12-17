module BitInterleaveTests
# NOTE: Use SafeTestSets: https://github.com/YingboMa/SafeTestsets.jl

using AdaptiveHierarchicalRegularBinning: bit_interleave, InterleaveMethod, Brute, Pdep
using Test

function bit_interleave_equivalence_test(n)
    ground_truth_method = Brute

    @testset "$(bit_interleave) equivalence" begin
        @testset "$M" for M in (Brute, Pdep)
            @testset "$T" for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
                W = rand(T, n)

                expected = bit_interleave(ground_truth_method, W)
                actual   = bit_interleave(M, W)

                @test expected == actual
            end
        end
    end
end


function bit_interleave_test(n)
    @testset "$(bit_interleave)" begin
        @testset "$T" for T in (UInt8, UInt16, UInt32, UInt64, UInt128)
            shift_reg(b, reg) = (reg<<1) | (b&0x1)
            toBits(wVec) = foldr(shift_reg, wVec; init=zero(T))

            # Create random
            WMat = rand(Bool, sizeof(T)*8, n)

            # BoolVec interleave
            R = Vector{Bool}(undef, prod(size(WMat)))
            for i in 1:n
                R[i:n:end] .= WMat[:, i]
            end
            R = R[1:size(WMat, 1)]
            expected = toBits(R)

            # Bit interleave
            W = map(toBits, eachslice(WMat, dims=2))
            actual = bit_interleave(W)

            @test actual == expected
        end
    end

end


function runtest(n)
    @testset "All Tests [$n]" begin
        bit_interleave_equivalence_test(n)
        bit_interleave_test(n)
    end
end

@testset "Bit Interleave" begin
    for n in 1:128
        runtest(n)
    end
end

end