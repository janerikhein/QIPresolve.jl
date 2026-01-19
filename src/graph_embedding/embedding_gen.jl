"""
    build_embedding_model(V, E, d2, R; r, k)

Build a JuMP model for integer graph embedding with squared edge lengths.
Returns `(model, x, y)` with integer coordinate variables.
"""
function build_embedding_model(V, E, d2, R; r, k)
    #TTODO: once XPRESS licence is there, use XPRESS to save as XPRESS LP file
    model = Model()

    # Integer coordinates with box bounds
    @variable(model, x[i in V], Int, lower_bound = -R[i], upper_bound =  R[i])
    @variable(model, y[i in V], Int, lower_bound = -R[i], upper_bound =  R[i])

    # objective: minimize 0
    @objective(model, Min, 0)

    # Edge distance equalities: (xi-xj)^2 + (yi-yj)^2 == d^2({i,j})
    # Assumes d2[(i,j)] is the squared distance (an integer/float).
    for (i, j) in E
        @constraint(model, (x[i] - x[j])^2 + (y[i] - y[j])^2 == d2[(i, j)])
    end

    # Radius constraint: xi^2 + yi^2 <= dist(r,i)   (here: R[i])
    for i in V
        @constraint(model, x[i]^2 + y[i]^2 <= R[i])
    end

    # Symmetry breaking / anchoring
    @constraint(model, x[r] == 0)
    @constraint(model, y[r] == 0)

    @constraint(model, y[k] - x[k] <= 0)
    @constraint(model, x[k] >= 0)
    @constraint(model, y[k] >= 0)

    return model, x, y
end

function generate_laman_instance()
    #TODO
end

function generate_random_instance()
    #TODO
end

