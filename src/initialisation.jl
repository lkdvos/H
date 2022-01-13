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

function su3_initial_state(spt)
    physical_spaces = [RepresentationSpace(SUNIrrep((3,0,0))=>1)]
    if spt == 0
        return InfiniteMPS(
            physical_spaces,
            [RepresentationSpace(SUNIrrep((0,0,0))=>5, SUNIrrep((2,1,0)) =>25, SUNIrrep((3,0,0))=>15, SUNIrrep((3,3,0))=>15, SUNIrrep((4,2,0))=>30, SUNIrrep((5,1,0))=>15, SUNIrrep((5,4,0))=>15, SUNIrrep((6,3,0))=>10, SUNIrrep((6,0,0))=>2, SUNIrrep((6,6,0))=>2, SUNIrrep((7,2,0))=>2, SUNIrrep((7,5,0))=>2, SUNIrrep((8,4,0))=>1)]
        )
    elseif spt == 1
        return InfiniteMPS(
            physical_spaces, 
            [RepresentationSpace(SUNIrrep((1,0,0))=> 25, SUNIrrep((2,2,0))=>30, SUNIrrep((3,1,0))=>50, SUNIrrep((4,0,0))=>20, SUNIrrep((4,3,0))=>41, SUNIrrep((5,2,0))=>40, SUNIrrep((5,5,0))=>10, SUNIrrep((6,1,0))=>10, SUNIrrep((6,4,0))=>20, SUNIrrep((7,3,0))=>10)]
        )
    elseif spt == 2
        return InfiniteMPS(
            physical_spaces,
            [RepresentationSpace(SUNIrrep((1,1,0))=>25, SUNIrrep((2,0,0))=>30, SUNIrrep((3,2,0))=>50, SUNIrrep((4,1,0))=>40, SUNIrrep((4,4,0))=>20, SUNIrrep((5,0,0))=>10, SUNIrrep((5,3,0))=>40, SUNIrrep((6,2,0))=>20, SUNIrrep((6,5,0))=>10, SUNIrrep((7,4,0))=>10)]
        )
    else
        @error "invalid spt"
    end
end
