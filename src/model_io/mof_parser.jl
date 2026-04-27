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

# Returns
- `nothing`.
"""
function save_moi(model::MOI.ModelLike, filename::AbstractString)
    MOI.write_to_file(model, filename)
    return nothing
end
