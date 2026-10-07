# Verbatim copy of the spider implementation from jeffwack/Mandelbrot.jl@619df42
# (src/spidermap.jl + src/spiderfuncs.jl), wrapped in a module so it can be
# benchmarked side by side with the upgraded implementation.
module LegacySpider

using Mandelbrot
using Mandelbrot: Sequence, KneadingSymbol, orbit, thetaitinerary, goestoidx
const partition = Mandelbrot.IterTools.partition

# Sequence indexing as it behaved at 619df42 (ignores the preperiod offset in the
# periodic part). Set FIX_GETINDEX[] = true to get the corrected behaviour.
const FIX_GETINDEX = Ref(false)
function legacy_getindex(S::Sequence, ii::Int)
    FIX_GETINDEX[] && return S[ii]
    ii <= S.preperiod && return S.items[ii]
    k = length(S.items) - S.preperiod
    return S.items[mod1(ii - S.preperiod, k)]
end

struct Spider
    angle::Rational
    orbit::Sequence
    kneading_sequence::Sequence{KneadingSymbol{4}}
    legs::Vector{Vector{ComplexF64}}
end

function Base.show(io::IO, S::Spider)
    npoints = sum([length(leg) for leg in S.legs])
    return println("A "*repr(S.angle)*"-spider with kneading sequence "*repr(S.kneading_sequence)*" and "*repr(npoints)*" points.")
end

function standardspider(theta)
    orb = orbit(theta)
    K = thetaitinerary(theta,orb) #the kneading sequence is calculated like so to avoid recalculating the orbit of theta
    legs = standardlegs(orb)
    return Spider(theta,orb,K,legs)
end

function standardlegs(orbit::Sequence)
    r = collect(LinRange(100,1,10))  #NOTE - it may be better to keep this as a 'linrange' but I don't understand what that means
    legs  = Vector{ComplexF64}[]
    for theta in orbit.items
        push!(legs,(cos(theta*2*pi)+1.0im*sin(theta*2*pi)) .* r)
    end
    #legs[1] .+= -legs[1][end]
    return legs
end

function grow!(legs::Vector{Vector{ComplexF64}},scale::Real,num::Int)
    #Get the arrow pointing from the second to last point to the last point
    for leg in legs
        arrow = scale*(leg[1]-leg[2])
        #add new elements which are the last element plus the arrow
        for i in 1:num
            pushfirst!(leg,leg[1]+arrow)
        end
    end
end


function prune!()
end

function stats(legs::Vector{Vector{ComplexF64}})
    println("~~~~~~~~~~~~~~~~~~~~")
    λ = legs[2][end]
    println(λ/2)
    points = sum(length.(legs))
    println(points)
    radius = sqrt(abs2(legs[1][1]))
    println(radius)
end

function sector(boundary1, boundary2, theta, phi)                                                                                                                                # Shift everything so boundary1 sits at 0
    b = mod(boundary2 - boundary1, 1)   # for a diameter, this is 0.5                                                                                                      
    t = mod(theta - boundary1, 1)
    p = mod(phi - boundary1, 1)                                                                                                                                            

    if p == b
        if t<b
            return KneadingSymbol{4}("*₂")
        else
            return KneadingSymbol{4}("*₁")
        end
    elseif p == 0
        if t<b
            return KneadingSymbol{4}("*₁")
        else
            return KneadingSymbol{4}("*₂")
        end
    end                                                                                                                                                                    
   
    # Same region iff both fall on the same side of b                                                                                                                      
    if (t < b) == (p < b)
        return KneadingSymbol{4}("A")
    else
        return KneadingSymbol{4}("B")
    end                                                                                                           
  end    
    
function spider_map(S::Spider)
    λ = S.legs[2][end] #the parameter of our polynomial map z -> λ(1+z/2)^2

    star1 = (angle(sqrt(S.legs[1][1]/λ))/(2*pi)+1)%1
    star2 = (star1 + 0.5)%1
    
    n = length(S.orbit.items)

    newLegs = Vector{Vector{ComplexF64}}(undef,n)

    newLegs[1] = reverse(path_sqrt(reverse(S.legs[2])./λ)) #reverse so that path_sqrt goes from foot to shoulder

    #println(abs(newLegs[1][end]-1.0)) #This should be zero
    
    thetaA = (angle(newLegs[1][1])/(2*pi)+1)%1 #gives the angle of the shoulder in full turns
    #this angle lays in region A by definition.
   
    sourceidxs = goestoidx(S.orbit)

    for ii in 2:n
        #first find the right half plane preimage of the shoulder
        sourceleg = S.legs[sourceidxs[ii]]
        u = sqrt(sourceleg[1]/λ)
        thetau = (angle(u)/(2*pi)+1)%1
        if legacy_getindex(S.kneading_sequence, ii) == sector(star1,star2,thetaA,thetau) #then this is the correct preimage
            newLegs[ii] = path_sqrt(sourceleg./λ)
        else
            newLegs[ii] = -1 .* path_sqrt(sourceleg./λ)
        end
    end
    
    for leg in newLegs
        leg .+= (-1.0+0.0im)
        leg .*= (2.0+0.0im)
    end
    
    grow!(newLegs,10,10)

    return Spider(S.angle,S.orbit,S.kneading_sequence,newLegs)
end
        
function spideriterates(S0::Spider,n_iter::Int)
    list = [S0]

    for i in 1:n_iter
        push!(list, spider_map(list[end]))
    end
    return list
end

"""
    parameter(S0::Spider, max_iter::Int)
    parameter(theta::Rational, max_iter::Int)

Computes the parameter (center) of a hyperbolic component using the spider algorithm 
described in [Hubbard_Schleicher_1995](@cite).

For a given external angle θ, this iteratively applies the spider map until convergence
to find the parameter c at the center of the corresponding hyperbolic component.

# References
- [Hubbard_Schleicher_1995](@cite): The spider algorithm
"""
#TODO modify below to have tolerance-based convergence behavior
function parameter(S::Spider,max_iter::Int)
    c_last = S.legs[2][end]/2 
    for ii in 1:max_iter
        S = spider_map(S)
        c = S.legs[2][end]/2
        #println(repr(c)*" delta "*repr(abs(c-c_last)))
        if abs(c-c_last)<(1e-15)
            return c
        end
        c_last = c
    end
    return S.legs[2][end]/2 
end

function parameter(theta::Rational,max_iter::Int)
    return parameter(standardspider(theta),max_iter)
end

#This function starts at the 'shoulder' (path[1]) and works in towards the 'foot' (path[end])
function path_sqrt(path::Vector{ComplexF64})
    branch = 1
    newpath = [sqrt(path[1])]
    cut = (0.0 + 0.0im, -1000 + 0.0im)
    for segment in partition(path,2,1)
        if test_intersection(cut...,segment...)
            branch *= -1
        end
        append!(newpath,branch*sqrt(segment[2]))
    end
    return newpath
end

function test_intersection(z1::Complex,z2::Complex,w1::Complex,w2::Complex)
    
    x1 = real(z1)
    y1 = imag(z1)

    x2 = real(z2)
    y2 = imag(z2)

    x3 = real(w1)
    y3 = imag(w1)

    x4 = real(w2)
    y4 = imag(w2)

    tn = (x1 - x3)*(y3 - y4) - (y1 - y3)*(x3 - x4)
    un = (x1 - x3)*(y1 - y2) - (y1 - y3)*(x1 - x2)
    d = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)

    if d == 0
        return false
    elseif sign(tn)*sign(d) == -1
        return false
    elseif sign(un)*sign(d) == -1
        return false
    elseif abs2(tn) > abs2(d)
        return false
    elseif abs2(un) > abs2(d)
        return false
    else
        return true
    end

end
end # module
