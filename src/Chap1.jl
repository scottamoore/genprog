module Chap1
using Printf, Random

export crossover, mutation, mutation_count, create_gene, test_mutation,
    test_mutation_count, crossover_exercise,
    test_generate_random_integer, random_number_tester, quartile_test

function create_gene(genelen)
    Array{Bool, 1}([rand(MersenneTwister(), Bool) for i in 1:genelen])
end
function create_gene(x::Array{Bool, 1})
    thelen = length(x)
    Array{Bool, 1}([x[i] for i in 1:thelen])
end

function crossover_exercise(population_size=200)
    half_population = div(population_size, 2)
    genelen = 5
    generation = 50
    thepop = Vector{Vector{Bool}}(undef, population_size)
    for i in 1:population_size
        thepop[i] = create_gene([false, false, false, false, false])
    end
    for i in 1:half_population
        thepop[i] = create_gene([true, true, true, false, false])
        thepop[i+half_population] = create_gene([false, false, false, true, true])
    end
    for g in 1:generation
        nextgen = Vector{Vector{Bool}}(undef, population_size)
        selection_order = randperm(MersenneTwister(), population_size)
        for c in 1:half_population
            crossover_site = rand(1:genelen-1)
            (nextgen[c], nextgen[c+half_population]) =
                crossover(thepop[selection_order[c]],
                    thepop[selection_order[c+half_population]],
                    crossover_site)
        end
        thepop = deepcopy(nextgen)
    end
    counts = zeros(Int16, 32)
    for i in 1:population_size
        #println("$i : $(thepop[i])")
        the_value = 0
        for p in 1:genelen
            the_value += thepop[i][p] ? 2^(p-1) : 0
        end
        counts[the_value+1] += 1
    end
    for i in 1:32
        println("$(i-1) : $(counts[i])")
    end
end

function mutation(g, prob)
    println(g)
    g = [rand() < prob ? !g[i] : g[i] for i in 1:length(g)]
    println("$(g)\n")
    return g
end

function mutation_count(g, prob)
    count::Int64 = 0
    for i in 1:length(g)
        if rand() < prob
            g[i] = !g[i]
            count += 1
        end
    end
    return (g, count)
end

function test_mutation()
    thelen = 10
    total_reps = 10
    odds_of_mutation = 0.1
    for reps in 1:total_reps
        gene = create_gene(thelen)
        #print("gene = $(gene)")
        newgene = mutation(gene, odds_of_mutation)
    end
end

function test_mutation_count()
    thelen = 100
    total_count = 0
    total_reps = 100000
    odds_of_mutation = 0.0001
    expected_mutations = thelen * total_reps * odds_of_mutation
    for reps in 1:total_reps
        gene = create_gene(thelen)
        #print("gene = $(gene)")
        (newgene, c) = mutation_count(gene, odds_of_mutation)
        total_count += c
        #println(" | $(newgene) | $(c)/$(thelen * odds_of_mutation)")
    end
    println("len=$(thelen) | reps=$(total_reps) | odds=$(odds_of_mutation) | exp=$(expected_mutations) | total=$(total_count)")
end

function crossover(s1::String, s2::String, site)
    if length(s1) == length(s2)
        if site >= 1
            if site <= length(s1)
                (s1[1:site] * s2[site+1:end], s2[1:site] * s1[site+1:end])
            else
                throw(DomainError(site, "'site' must be less than $(length(s2))"))
            end
        else
            throw(DomainError(site, "'site' must be >= 1"))
        end
    else
        throw(DomainError("The strings must be of the same length"))
    end
end
function crossover(s1::Array{<:Bool, 1}, s2::Array{<:Bool, 1}, site)
    if length(s1) == length(s2)
        if site >= 1
            if site < length(s1)
                c1 = deepcopy(s1)
                retval1 = c1[1:site]
                c2 = deepcopy(s2)
                retval2 = c2[1:site]
                append!(retval1, s2[site+1:end])
                append!(retval2, s1[site+1:end])
            else
                throw(DomainError(site, "'site' must be less than $(length(s2))"))
            end
        else
            throw(DomainError(site, "'site' must be >= 1"))
        end
    else
        throw(DomainError("The vectors must be of the same length"))
    end
    return (retval1, retval2)
end


function generate_random_integer(lower::Int64, upper::Int64)
    rand(lower:upper)
end

function test_generate_random_integer(lower::Int64, upper::Int64)
    num_tests = 10_000_000
    values_in_range = upper - lower + 1
    counter = zeros(Int64, values_in_range)
    num_to_subtract = lower - 1
    expected = num_tests/values_in_range
    for i in 1:num_tests
        counter[rand(lower:upper) - num_to_subtract] += 1
    end
    println("\ntest_generate_random_integer()")
    for i in lower:upper
        @printf("%d : %d : %.6f\n",
            i,
            counter[i-num_to_subtract],
            (counter[i-num_to_subtract] - expected)/expected)
    end
end

function random_number_tester(whichtest, num_tests::Int64)
    res::Float64 = 0.0
    if whichtest == 1
        thesum::Float64 = 0.0
        for i in (rand() for x=1:num_tests)
            thesum += i
        end
        res = thesum/num_tests
    elseif whichtest == 2
        res = sum(rand() for x=1:num_tests)/num_tests
    end
end

function quartile_test()
    counter = zeros(Int64, 4)
    num_tests = 10_000_000
    println("\nquartile_test()")
    for i in (rand() for x=1:num_tests)
        if i <= 0.25
            counter[1] += 1
        elseif i <= 0.5
            counter[2] += 1
        elseif i <= 0.75
            counter[3] += 1
        else
            counter[4] += 1
        end
    end
    expected = div(num_tests, 4)
    for i in 1:4
        @printf("%d : %d : %.6f \n",
            i, counter[i], (counter[i] - expected)/expected)
    end
end

end
