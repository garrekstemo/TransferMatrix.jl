function add_layer!(s::Structure, i::Int64, layer::Layer)
    insert!(s.layers, i, layer)
    initialize!(s.λ, s)
end

function add_layer!(s::Structure, layer::Layer)
    push!(s.layers, layer)
    initialize!(s.λ, s)
end

function delete_layer!(s::Structure, i::Int64)
    if i == 1
        popfirst!(s.layers)

    elseif i == lastindex(s.layers)
        pop!(s.layers)

    else
        popat!(s.layers, i)

    end
    return s
end

function delete_layer!(s::Structure)
    pop!(s.layers)
    return s
end