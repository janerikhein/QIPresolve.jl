"""
    PresolveConfig

Define package-level configuration hooks for presolve features.

The constants in this module define the package-level defaults for combined
presolve behavior.
"""
module PresolveConfig

export DEFAULT_PRESOLVE_RESIDUE_STRATEGY,
    DEFAULT_PRESOLVE_RESIDUE_THRESHOLD,
    DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD

const DEFAULT_PRESOLVE_RESIDUE_STRATEGY = :divisor_free
const DEFAULT_PRESOLVE_RESIDUE_THRESHOLD = 64
const DEFAULT_PRESOLVE_TREEWIDTH_THRESHOLD = 2

end # module
