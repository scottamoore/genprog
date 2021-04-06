module Chap3
using Printf, Random

export sga

mutable struct Individual
    chrom::Vector{Bool}
    length::Int16
    x::Float64
    fitness::Float64
    parent1::Int16
    parent2::Int16
    xsite::Int16
end

"""
    create_individual(thelen[, is_empty=false])

Create an individual for the population.

If 'is_empty' is true, then do not create values for the chromosome.

# Examples
```julia-repl
julia> create_individual(5)

julia> create_individual(20, true)

julia> create_individual(30, )
```
"""
function create_individual(thelen::Signed, is_empty::Bool = false)
    Individual(create_gene(thelen, is_empty),
        thelen, 0.0, 0.0, 0, 0, 0)
end
function create_individual(thelen::Signed, chrom::Vector{Bool})
    theval = decode(chrom, thelen)
    Individual(create_gene(chrom),
        thelen, theval,
        objfunc(theval),
        0, 0, 0)
end

function create_gene(genelen::Signed, is_empty::Bool = false)
    if is_empty
        Vector{Bool}([false for i in 1:genelen])
    else
        Vector{Bool}([rand(MersenneTwister(), Bool) for i in 1:genelen])
    end
end
function create_gene(x::Vector{Bool})
    thelen = length(x)
    Vector{Bool}([x[i] for i in 1:thelen])
end

function mutation_count(g, prob)
    count::Int16 = 0
    for i in 1:length(g)
        if rand() < prob
            g[i] = !g[i]
            count += 1
        end
    end
    return (g, count)
end

rand_crossover_site(thelen::Signed) = Int16(rand(1:thelen-1))

function do_crossover(s1::Vector{Bool}, s2::Vector{Bool}, site::Signed)
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

function crossover(s1::Vector{Bool}, s2::Vector{Bool}, site::Signed,
            len_of_gene::Signed,
            prob_of_crossover::Float64,
            prob_of_mutation::Float64)
    if rand() < prob_of_crossover
        num_crossovers = 1
        (child_chrom1, child_chrom2) = do_crossover(s1, s2, site)
    else
        num_crossovers = 0
        (child_chrom1, child_chrom2) = (s1, s2)
    end
    (child_chrom1, num_mutations) = mutation_count(child_chrom1, prob_of_mutation)
    return (child_chrom1, child_chrom2, num_crossovers, num_mutations)
end

function select(popsize::Signed,
            sumfitness::Float64,
            population::Vector{Individual})
    partsum::Float64 = 0.0
    indiv_choice::Int16 = 0
    selector_value::Float64 = rand() * sumfitness
    while partsum < selector_value
        indiv_choice += 1
        partsum += population[indiv_choice].fitness
    end
    return indiv_choice
end

empty_population(popsize::Signed) = Vector{Individual}(undef, popsize)

function initialize_population(population_size::T,
            thelen::T) where {T <: Signed}
    thepop = empty_population(population_size)
    for i in 1:population_size
        thepop[i] = create_individual(thelen, false)
        thepop[i].length = thelen
        thepop[i].x = decode(thepop[i].chrom, thelen)
        thepop[i].fitness = objfunc(thepop[i].x)
    end
    return thepop
end

function decode(chrom::Vector{Bool}, thelen)
    sum((chrom[i] ? 2^(i-1) : 0 for i in 1:thelen))
end

objfunc(x) = (x/1_073_741_823.0)^10

function population_fitness(p::Vector{Individual}, popsize::Signed)
    sum((p[i].fitness for i in 1:popsize))
end

function min_population_fitness(p::Vector{Individual}, popsize::Signed)
    min((p[i].fitness for i in 1:popsize)...)
end

function max_population_fitness(p::Vector{Individual}, popsize::Signed)
    max((p[i].fitness for i in 1:popsize)...)
end

function generation(oldpop::Vector{Individual},
            maxpopsize::Signed,
            maxgenelength::Signed,
            prob_of_crossover::Float64,
            prob_of_mutation::Float64)
    i::Int16 = 0
    total_crossovers::Int16 = 0
    total_mutations::Int16 = 0
    newpop = empty_population(maxpopsize)
    sumfitness = population_fitness(oldpop, maxpopsize)
    while i < maxpopsize
        mate1 = select(maxpopsize, sumfitness, oldpop)
        mate2 = select(maxpopsize, sumfitness, oldpop)
        site = rand_crossover_site(maxgenelength)
        (child1_chrom, child2_chrom, new_crossovers, new_mutations) =
            crossover(oldpop[mate1].chrom, oldpop[mate2].chrom,
                site, maxgenelength,
                prob_of_crossover, prob_of_mutation)
        newpop[i+1] = create_individual(maxgenelength, child1_chrom)
        newpop[i+1].parent1 = mate1
        newpop[i+1].parent2 = mate2
        newpop[i+1].xsite = site
        newpop[i+2] = create_individual(maxgenelength, child2_chrom)
        newpop[i+2].parent1 = mate1
        newpop[i+2].parent2 = mate2
        newpop[i+2].xsite = site
        total_crossovers += new_crossovers
        total_mutations += new_mutations
        i += 2  # the new population now has 2 new members
    end
    return (newpop, total_crossovers, total_mutations)
end

function statistics(thepop::Vector{Individual}, popsize::Signed)
    popfitness = population_fitness(thepop, popsize)
    minfitness = min_population_fitness(thepop, popsize)
    maxfitness = max_population_fitness(thepop, popsize)
    avgfitness = popfitness/popsize
    return (minfitness, maxfitness, avgfitness)
end

function print_results(id::String,
            newpop::Vector{Individual},
            current_gen::Signed,
            popsize::Signed,
            minfitness, maxfitness, avgfitness)
    @printf("%s | %.10f | %.10f | %.10f\n",
        id, minfitness, avgfitness, maxfitness)
end

function sga(id::Signed,
            maxpopsize::T = 100,
            maxgenelength::T = 30,
            maxgen::T = 50,
            prob_of_crossover::Float64 = 0.8,
            prob_of_mutation::Float64 = 0.0001) where {T <: Signed}
    total_crossovers::Int16 = 0
    total_mutations::Int16 = 0
    oldpop = initialize_population(maxpopsize, maxgenelength)
    newpop = Vector{Individual}(undef, maxpopsize)
    println("\nID  |     Min      |      Avg     |      Max ")
    for current_gen in 1:maxgen
        unique_generation_id = @sprintf("%d-%d", id, current_gen)
        (newpop, new_crossovers, new_mutations) =
            generation(oldpop, maxpopsize, maxgenelength,
                prob_of_crossover, prob_of_mutation)
        (minfitness, maxfitness, avgfitness) = statistics(newpop, maxpopsize)
        print_results(unique_generation_id, newpop, current_gen, maxpopsize,
            minfitness, maxfitness, avgfitness)
        oldpop = deepcopy(newpop)
    end
    println("Settings: id=$(id)  popsize=$(maxpopsize)  gens=$(maxgen)  xover=$(prob_of_crossover)  mut=$(prob_of_mutation)")
end

end
