const _VALID_PARITY_STRATEGIES = (:mod2_basic, :mod4_basic, :full)

function _parity_strategy_error(parity_strategy)
    return ArgumentError(
        "Invalid parity_strategy: $(repr(parity_strategy)). " *
        "Expected one of :mod2_basic, :mod4_basic, or :full.",
    )
end

function _normalize_parity_strategy(parity_strategy)::Symbol
    normalized = if parity_strategy isa Symbol
        Symbol(replace(lowercase(String(parity_strategy)), "-" => "_"))
    elseif parity_strategy isa AbstractString
        Symbol(replace(lowercase(strip(parity_strategy)), "-" => "_"))
    else
        throw(_parity_strategy_error(parity_strategy))
    end

    normalized in _VALID_PARITY_STRATEGIES || throw(_parity_strategy_error(parity_strategy))
    return normalized
end

_parity_includes_mod4(parity_strategy::Symbol) = parity_strategy != :mod2_basic
_parity_uses_full_propagation(parity_strategy::Symbol) = parity_strategy == :full
