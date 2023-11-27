"""
    Properties

Photosynthetic, growth, and death properties.

# Fields
- `maint::Float64`: assimilates for maintenance (kg CH2O/palm/day)
- `frac::Float64`: DM (dry matter) partitioning (fraction)
- `growth::Float64`: growth rate (kg DM/palm/day)
- `death::Float64`: death rate (kg DM/palm/day)
- `weight::Float64`: weight of the tree part (kg DM/palm)
"""
@with_kw mutable struct Properties
    maint::Float64 = 0.0
    frac::Float64 = 0.0
    growth::Float64 = 0.0
    death::Float64 = 0.0
    weight::Float64 = 0.0
end


"""
    Parts

Tree parts and their respective `Properties`.

# Fields
- `pinnae::Properties`: pinnae properties (vegetative)
- `rachis::Properties`: rachis properties (vegetative)
- `trunk::Properties`: trunk properties (vegetative)
- `roots::Properties`: roots properties (vegetative)
- `maleflo::Properties`: male flowers properties (flower)
- `femaflo::Properties`: female flowers properties (flower)
- `bunches::Properties`: bunches properties (bunches)
"""
@with_kw struct Parts
    pinnae::Properties = Properties()    # 1, vegetative
    rachis::Properties = Properties()    # 2, vegetative
    trunk::Properties = Properties()     # 3, vegetative
    roots::Properties = Properties()     # 4, vegetative
    maleflo::Properties = Properties()   # 5, flowers
    femaflo::Properties = Properties()   # 6, flowers
    bunches::Properties = Properties()   # 7, bunches
end


"""
    getfields(allparts::Bool)

Return the `Parts` fields for either vegetative or all tree parts.

# Arguments
- `allparts::Bool`: `true` to return all tree parts, else `false` for only
                    vegetative tree parts

# Return
- `Vector{Symbol}`: list of `Parts` fields
"""
function getfields(allparts::Bool)
    n = allparts ? 7 : 4
    [:pinnae, :rachis, :trunk, :roots, :maleflo, :femaflo, :bunches][1:n]
end


"""
    vegparts(parts::Parts, prop::Symbol)

Return a given property for vegetative tree parts.

# Arguments
- `parts::Parts`: all tree parts
- `prop::Symbol`: property of tree parts to retrieve

# Returns
- `Vector{Float64}`: property values from every vegetative tree part
"""
function vegparts(parts::Parts, prop::Symbol)
    getproperty.(getproperty.(Ref(parts), getfields(false)), prop)
end


"""
    allparts(parts::Parts, prop::Symbol)

Return a given property for all tree parts.

# Arguments
- `parts::Parts`: all tree parts
- `prop::Symbol`: property of tree parts to retrieve

# Returns
- `Vector{Float64}`: property values from every tree part
"""
function allparts(parts::Parts, prop::Symbol)
    getproperty.(getproperty.(Ref(parts), getfields(true)), prop)
end


"""
    vegparts!(parts::Parts, prop::Symbol, val)

Set a given property for vegetative tree parts.

# Arguments
- `parts::Parts`: all tree parts
- `prop::Symbol`: property of tree parts to retrieve
- `val`: value to set

# Returns
- nothing
"""
function vegparts!(parts::Parts, prop::Symbol, val)
    setproperty!.(getproperty.(Ref(parts), getfields(false)), prop, val)
end


"""
    allparts!(parts::Parts, prop::Symbol, val)

Set a given property for all tree parts.

# Arguments
- `parts::Parts`: all tree parts
- `prop::Symbol`: property of tree parts to retrieve
- `val`: value to set

# Returns
- nothing
"""
function allparts!(parts::Parts, prop::Symbol, val)
    setproperty!.(getproperty.(Ref(parts), getfields(true)), prop, val)
end
