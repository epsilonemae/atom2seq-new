loop over atoms {
    if atom is already in a group {skip this atom}
    else {
        detect Ind {
            check if atom is C, N, or H
            if C, check if it is bonded to three things, one is C {
                others are both C
                others are C and H
                others are N and H
                others are C and N
            }
            if N, check if it is bonded to two Cs and an H
            if H, check if it is bonded to a C or an N
        }
        detect Imd {
            check if atom is C, N, or H
            if C, check if it is bonded to three things, one is N {
                others are N and H
                others are C and H
                others contain a C
            }
            if N, check if it is bonded to two Cs {
                if bonded to an H, it is in the group
            }
            if H, check if it is bonded to N or C
            verify other atoms in this group
        }
        detect PhOH {
            check if atom is C, O, or H
            if C, check if it is connected to three things, two are Cs, and the remaining is O or H
            if O, check if it is connected to a C and an H
            if H, check if it is connected to an O or a C
            verify other atoms in group
        }
        detect Ph {
            check if atom is C or H
            if C, check if it is connected to three things and two are Cs {
                if it is connected to a H, that H is in the group
            }
            if H, check if it is connected to a C
            verify other atoms in group
        }
        detect Amd {
            check if atom is C, O, N, or H
            if C, check if it is connected to three things: one is a O and one is a N
            if O, check if it is connected to only a C
            if N, check if it is connected to three things: one is a C and one is an H
            if H, check if it is connected to an N
            verify other atoms in group
        }
        detect COOH {
            check if atom is C, O, or H
            if C, check if it is connected to three things and two are Os
            if O, check if it is connected to either {
                only a C
                a C and an H
            }
            if H, check it is connected to an O
            verify other atoms in group

        }
        detect C=O {
            check if atom is a C or O
            if C, check it is connected to three things and one is a O
            if O, check it is connected to only a C
            verify other atoms in group
        }
        detect XHn cluster {
            if H, find the thing it's connected to {
                run XHn cluster on that
            }
            if not H, bond it to all Hs connected to it (can be none)
        }
    }
}