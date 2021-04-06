using LoopVectorization

const SIZE = 100_000
const a = rand(SIZE)
const b = rand(SIZE)

```
The base case.
```
test1(a, b) = a .> b

```
A vector of Booleans in which we do not
check the validity of the array index
and also pre-allocate the memory for
the result array.
```
function test2(a, b, the_size)
    c = Vector{Bool}(undef, the_size)
    @inbounds for i in 1:the_size
        c[i] = a[i] > b[i]
    end
    c
end

'''
This is the fastest one.
A Vector of Booleans takes up 8x the memory compared with a BitVector.
Times in μs:
test1: 44.9
test2: 35.9
test3: 24.1
test3a: 220.7
test4: 37.5
eltype(c) == Bool
'''
function test3(a, b, the_size)
    c = BitVector(undef, the_size)
    @avx for i in 1:the_size
        c[i] = a[i] > b[i]
    end
    c
end

function test3a(a, b, the_size)
    c = BitVector(undef, the_size)
    for i in 1:the_size
        c[i] = a[i] > b[i]
    end
    c
end

function test4(a, b, the_size)
    c = Vector{Bool}(undef, the_size)
    @avx for i in 1:the_size
        c[i] = a[i] > b[i]
    end
    c
end
