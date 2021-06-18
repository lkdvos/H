function expand_params(params::Dict)
    # Credits to DrWatson.jl -- _dict_list
    iterable_fields = filter(x -> typeof(params[x]) <: Vector, keys(params))
    non_iterables = setdiff(keys(params), iterable_fields)
    
    iterable_dict = Dict(iterable_fields .=> getindex.(Ref(params), iterable_fields))
    non_iterable_dict = Dict(non_iterables .=> getindex.(Ref(params), non_iterables))

    vec(
        map(Iterators.product(values(iterable_dict)...)) do vals
            dd = [k=>convert(eltype(params[k]),v) for (k,v) in zip(keys(iterable_dict),vals)]
            if isempty(non_iterable_dict)
                Dict(dd)
            elseif isempty(iterable_dict)
                non_iterable_dict
            else
                # We can't use merge here because it promotes types.
                # The uniqueness of the dictionary keys is guaranteed.
                Dict(dd..., collect(non_iterable_dict)...)
            end
        end
    )
    
end
