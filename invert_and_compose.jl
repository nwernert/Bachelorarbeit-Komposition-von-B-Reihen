using BSeries
using RootedTrees
using Latexify
using Test
import SymPyPythonCall; sp = SymPyPythonCall;
import BSeries.compose
import BSeries.SplittingIterator
import RootedTrees.binary_digits!
import RootedTrees._subtree_last_index

#Implementation der Inversen fuer B-Reihen basierend auf zwei- und einfarbigen Wurzelbäumen.
function invert(a :: TruncatedBSeries)
    series_keys = keys(a)
    a_inv = empty(a, keytype(a), valtype(a))

    for t in series_keys
        coefficient = invert(a, a_inv, t)
        a_inv[t] = coefficient
    end

    return a_inv
end

function invert(a :: TruncatedBSeries, a_inv :: TruncatedBSeries, t::RootedTree)
    if(t == rootedtree(Int[]))
        return 1
    end

    result = zero(first(values(a)))

    for (forest, subtree) in SplittingIterator(t)
        if subtree != t
            update = a_inv[subtree]
            update isa Union{Rational, Bool} && iszero(update) && continue
            for tree in forest
                update *= a[tree]
            end
            result -= update
        end
    end

    return result
end

function invert(a :: TruncatedBSeries, a_inv :: TruncatedBSeries, t::BicoloredRootedTree)
    if(t == rootedtree(Int[], Bool[]))
        return 1
    end

    result = zero(first(values(a)))

    for (forest, subtree) in ColoredSplittingIterator(t)
        if subtree != t
            update = a_inv[subtree]
            update isa Union{Rational, Bool} && iszero(update) && continue
            for tree in forest
                update *= a[tree]
            end
            result -= update
        end
    end

    return result
end

#Implementation einer Verallgemeinerung von compose(b, a; normalize_stepsize = false) für möglicherweise unterschiedliche Schrittweiten.
function compose(b, a, factor_b, factor_a)
    series_keys = keys(b)
    @assert keytype(b) == keytype(a)
    V = promote_type(valtype(b), valtype(a))
    series = empty(b, keytype(b), V)

    b = rescale_stepsize(b, factor_b)
    a = rescale_stepsize(a, factor_a) 

    for t in series_keys
        coefficient = compose(b, a, t)
        series[t] = coefficient
    end

    return series
end

function rescale_stepsize(a :: TruncatedBSeries, factor)
    series_keys = keys(a)
    a_new = empty(a, keytype(a), valtype(a))

    for t in series_keys
        a_new[t] = a[t] * factor^order(t)
    end

    return a_new
end

#Implementation von compose fuer B-Reihen fuer partitionierte ODEs der Form u'(t) = f(u(t))+g(u(t)).
function compose(b, a, t::BicoloredRootedTree)
    result = zero(first(values(a)) * first(values(b)))

    for (forest, subtree) in ColoredSplittingIterator(t)
        update = a[subtree]
        update isa Union{Rational, Bool} && iszero(update) && continue
        for tree in forest
            update *= b[tree]
        end
        result += update
    end

    return result
end
#der zugehoerige SplittingIterator fuer farbige Baeume
struct ColoredSplittingIterator{T <: ColoredRootedTree}
    t::T
    node_set::Vector{Bool}
    max_node_set_value::Int

    function ColoredSplittingIterator(t::T) where {T <: ColoredRootedTree}
        node_set = zeros(Bool, order(t))
        return new{T}(t, node_set, 2^order(t) - 1)
    end
end

Base.IteratorSize(::Type{<:ColoredSplittingIterator}) = Base.SizeUnknown()
Base.eltype(::Type{ColoredSplittingIterator{T}}) where {T} = Tuple{Vector{T}, T}

@inline function Base.iterate(splittings::ColoredSplittingIterator)
    node_set_value = 0
    return iterate(splittings, node_set_value)
end

@inline function Base.iterate(splittings::ColoredSplittingIterator, node_set_value)
    node_set_value > splittings.max_node_set_value && return nothing

    node_set = splittings.node_set
    t = splittings.t
    ls = t.level_sequence
    cs = t.color_sequence
    T = eltype(ls)
    forest = Vector{typeof(t)}()

    while node_set_value <= splittings.max_node_set_value
        binary_digits!(node_set, node_set_value)

        # Check that if a node is removed then all of its descendants are removed
        subtree_root_index = 1
        empty!(forest)
        while subtree_root_index <= order(t)
            if node_set[subtree_root_index] == false # This node is removed
                subtree_last_index = _subtree_last_index(subtree_root_index, ls)

                # Check that subtree is all removed
                if !any(@view node_set[subtree_root_index:subtree_last_index])
                    # If `iscanonical(t)`, the subtree starting at the root of `t`
                    # is also in canonical representation. Thus, we don't need to
                    # use the more expensive version
                    #   push!(forest, rootedtree!(level_sequence))
                    # but can use the cheaper version below.
                    level_sequence = ls[subtree_root_index:subtree_last_index]
                    color_sequence = cs[subtree_root_index:subtree_last_index]
                    push!(forest, ColoredRootedTree(level_sequence, color_sequence, RootedTrees.iscanonical(t)))
                    subtree_root_index = subtree_last_index + 1
                else
                    break
                end
            else
                subtree_root_index += 1
            end
        end

        if subtree_root_index == order(t) + 1
            # This is a valid ordered subtree.
            # The `level_sequence` will not automatically be a canonical representation.
            # TODO: splittings;
            #       Decide whether canonical representations should be used. Disabling
            #       them will increase the performance.
            level_sequence = ls[node_set]
            color_sequence = cs[node_set]
            subtree = rootedtree!(level_sequence, color_sequence)
            return ((forest, subtree), node_set_value + 1)
        else
            node_set_value = node_set_value + 1
        end
    end

    return nothing
end


#helper for tests to make the identity b-series for additive bseries
function turnIntoIdentity(a :: TruncatedBSeries)
    series_keys = keys(a)
    id = empty(a, keytype(a), valtype(a))
    emptyT = rootedtree(Int[],Bool[])

    for t in series_keys
        coefficient = 0
        if t == emptyT
            coefficient = 1
        end
        id[t] = coefficient
    end

    return id
end
#helper for tests that randomize the B-series while keeping a(emptyset) = 1
#randomizer for B-series for additive B-series
function randomizeAdd(a:: TruncatedBSeries)
    series_keys = keys(a)
    rando = empty(a, keytype(a), valtype(a))
    emptyT = rootedtree(Int[],Bool[])

    for t in series_keys
        coefficient = rand(-100:100)
        if t == emptyT
            coefficient = 1
        end
        rando[t] = coefficient
    end

    return rando
end
#randomizer for ordinary B-series
function randomize(a:: TruncatedBSeries)
    series_keys = keys(a)
    rando = empty(a, keytype(a), valtype(a))
    emptyT = rootedtree(Int[])

    for t in series_keys
        coefficient = rand(-100:100)
        if t == emptyT
            coefficient = 1
        end
        rando[t] = coefficient
    end

    return rando
end


#tests fuer invert auf gewoehnlichen rooted trees
α = sp.symbols("α", real = true)
A = [0 0; 1/(2*α) 0]; b = [1-α, α]; c = [0, 1/(2*α)]
coeffs2 = bseries(A, b, c, 9)
inverse = invert(coeffs2)
identity = bseries(IdentityMap{Rational{Int}}(), 9)

@test invert(inverse) == coeffs2
@test compose(inverse, coeffs2) == identity
@test compose(coeffs2, inverse) == identity

rando = randomize(coeffs2)
inverse = invert(rando)

@test invert(inverse) == rando
@test compose(inverse, rando) == identity
@test compose(rando, inverse) == identity

#tests fuer invert/compose auf bicolored rooted trees
As = [[0 0; 1//2 1//2],[1//2 0; 1//2 0]]; bs = [[1//2, 1//2],[1//2, 1//2]]
ark = AdditiveRungeKuttaMethod(As, bs)
arkS = bseries(ark, 9)
inverse = invert(arkS)
identity = turnIntoIdentity(arkS)

@test invert(inverse) == arkS
@test compose(inverse, arkS) == identity
@test compose(arkS, inverse) == identity

rando = randomizeAdd(arkS)
inverse = invert(rando)

@test invert(inverse) == rando
@test compose(inverse, rando) == identity
@test compose(rando, inverse) == identity


#tests fuer compose auf bicolored rooted trees und tests fuer compose(b, a, factor_b, factor_a)
#strang splitting tests
#explizite mittelpunktsregel, ordnung p=2
F = [0 0; 1//2 0]
b_f = [0, 1]
#explizite Trapezregel, ordnung p=2
G = [0 0; 1 0]
b_g = [1//2, 1//2]

As_f = [F, [0 0; 0 0]]
bs_f = [b_f, [0, 0]]
rk_f = AdditiveRungeKuttaMethod(As_f, bs_f)
f = bseries(rk_f, 9)
As_g = [[0 0; 0 0], G]
bs_g = [[0, 0], b_g]
rk_g = AdditiveRungeKuttaMethod(As_g, bs_g)
g = bseries(rk_g, 9)

@test compose(f,g,1//2,1//2) == compose(f,g,normalize_stepsize = true)
@test compose(g,f,1//2,1//2) == compose(g,f,normalize_stepsize = true)

fg = compose(f,g, 1//2, 1)
fgf = compose(fg, f, 1, 1//2)
@test order_of_accuracy(fgf) == 2

gf = compose(g,f,1,1//2)
fgf = compose(f,gf,1//2,1)
@test order_of_accuracy(fgf) == 2


#test fuer compose(b, a, factor_b, factor_a) fuer einfarbige Baeume
rando1 = randomize(coeffs2)
rando2 = randomize(coeffs2)
@test compose(rando1, rando2, normalize_stepsize = true) == compose(rando1, rando2, 1//2, 1//2)
