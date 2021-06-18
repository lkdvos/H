function su2_initial_state(spin, spt)
    physical_spaces = [Rep[SU₂](spin => 1)]
    if spt == 0
        return InfiniteMPS(
            physical_spaces, 
            [Rep[SU₂](0 => 20, 1 => 20, 2 => 10, 3 => 5)]
        )
    elseif spt == 1
        return InfiniteMPS(
            physical_spaces,
            [Rep[SU₂](1//2 => 20, 3//2 => 20, 5//2 => 10, 7//2 => 3)]
        )
    end
end