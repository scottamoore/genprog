include("Chap1.jl")
include("Chap3.jl")

module test
using Main.Chap1, Main.Chap3, Random

@time sga(1)
@time sga(2, 100, 30, 50, 0.95, 0.001)
@time sga(3, 100, 30, 100, 0.95, 0.001)
@time sga(4, 100, 30, 100, 0.6, 0.001)

#=
println("\nrandom_number_tester()")
for the_test in 1:2
    val = random_number_tester(the_test, 10_000_000)
    stats = @timed random_number_tester(the_test, 10_000_000)
    println("$(the_test) avg=$(val) time=$(stats.time)")
end

quartile_test()

test_generate_random_integer(3, 12)

println(crossover("asdfghij", "qwertyui", 2))
println(crossover("asdfghij", "qwertyui", 5))
println(crossover("asdfghij", "qwertyui", 1))
println(crossover("asdfghij", "qwertyui", 20))
println(crossover("asdfghij", "qwertyui", -3))

genelen = 20
a1 = create_gene(genelen)
println(a1)
a2 = create_gene(genelen)
println(a2)
println(crossover(a1, a2, 2))
println(crossover(a1, a2, 5))
println(crossover(a1, a2, 1))
try
    println(crossover(a1, a2[1:end-2], 3))
catch ex
    if ex isa DomainError
        println("DomainError (code = $(ex.val) $(ex.msg)")
    else
        println("Unknown error")
    end
end

try
    println(crossover(a1, a2, 25))
catch ex
    if ex isa DomainError
        println("DomainError (code = $(ex.val) $(ex.msg)")
    else
        println("Unknown error")
    end
end

try
    println(crossover(a1, a2, -3))
catch ex
    if ex isa DomainError
        println("DomainError (code = $(ex.val) $(ex.msg)")
    else
        println("Unknown error")
    end
end
test_mutation_count()

@time test_mutation_count()

crossover_exercise(50)
@time crossover_exercise(100)
@time crossover_exercise(1000)
@time crossover_exercise(10000)
@time crossover_exercise(100_000)
=#

end
