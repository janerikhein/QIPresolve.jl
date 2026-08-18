import MathOptInterface as MOI
import MathOptInterface.FileFormats as FF

"""
    load_moi_model(file_path; format=:mof)

Load an MOI model from disk.

# Arguments
- `file_path`: Path to the serialized optimization model.

# Keywords
- `format=:mof`: File format identifier. Only `:mof` is currently supported.

# Returns
- An `MOI.ModelLike` instance populated from the file contents.

# Throws
- `ArgumentError` if `format` is unsupported.
"""
function load_moi_model(file_path::AbstractString; format::Symbol = :mof)
    if format == :mof
        ff = FF.FORMAT_MOF
    else
        throw(ArgumentError("invalid format identifier: $format"))
    end

    moi = FF.Model(format = ff)
    MOI.read_from_file(moi, file_path)

    return moi
end

"""
    save_moi(model, filename)

Write `model` to `filename` using MathOptInterface file I/O.

When writing LP files, split ranged scalar-function constraints into two
one-sided constraints so the generated LP can be read back by MOI's LP reader.
Variable interval bounds are left unchanged because they belong in the LP
`Bounds` section.

# Returns
- `nothing`.
"""
function save_moi(model::MOI.ModelLike, filename::AbstractString)
    output_model = lowercase(splitext(filename)[2]) == ".lp" ?
        _lp_reader_compatible_model(model) :
        model
    MOI.write_to_file(output_model, filename)
    return nothing
end

function _lp_reader_compatible_model(model::MOI.ModelLike)
    output_model = MOI.Utilities.Model{Float64}()
    MOI.copy_to(output_model, model)
    _split_lp_interval_function_constraints!(output_model, MOI.ScalarAffineFunction{Float64}, "affine")
    _split_lp_interval_function_constraints!(output_model, MOI.ScalarQuadraticFunction{Float64}, "quadratic")
    return output_model
end

function _split_lp_interval_function_constraints!(
        model::MOI.ModelLike,
        ::Type{F},
        fallback_prefix::AbstractString,
    ) where {F}
    indices = copy(MOI.get(model, MOI.ListOfConstraintIndices{F, MOI.Interval{Float64}}()))
    for (idx, ci) in enumerate(indices)
        func = MOI.get(model, MOI.ConstraintFunction(), ci)
        set = MOI.get(model, MOI.ConstraintSet(), ci)
        name = MOI.get(model, MOI.ConstraintName(), ci)
        isempty(name) && (name = "$(fallback_prefix)_ranged_$(idx)")

        MOI.delete(model, ci)

        lower_ci = MOI.add_constraint(model, deepcopy(func), MOI.GreaterThan(set.lower))
        upper_ci = MOI.add_constraint(model, deepcopy(func), MOI.LessThan(set.upper))
        MOI.set(model, MOI.ConstraintName(), lower_ci, "$(name)_lb")
        MOI.set(model, MOI.ConstraintName(), upper_ci, "$(name)_ub")
    end
    return model
end
