func carbonyl_finder(atom) {
    if atom is O {
        find C bonded to O
        group(O, C)
    } elif atom is C {
        find O bonded to only this C
        group(O, C)
    }
}

func cooh_finder(atom, is_initial) {
    if is_initial {
        to_group = []
        to_group.append(atom)
        find C bonded to atom
        to_group.append(C)
        find (O bonded to H) bonded to C
        to_group.append(O)
        to_group.append(H)
        group(to_group)
    }
    else {
        if atom is H {
            find bonded
            cooh_finder(bonded, False)
        } elif atom is O {
            if bonded to only a C {
                cooh_finder(atom, True)
            } elif bonded to H {
                find C bonded to O
                cooh_finder(C, False)
            }
        } elif atom is C {
            find O bonded to only C
            cooh_finder(O, True)
        }
    }
}

func amd_finder(atom, is_initial) {
    if is_initial{
        to_group = []
        to_group.append(atom)
        find C bonded to atom
        to_group.append(C)
        find N bonded to C
        to_group.append(N)
        find Hs bonded to N
        for H in Hs:
            to_group.append(H)
    }
    else {
        if atom is H {
            find bonded
            amd_finder(bonded, False)
        } elif atom is O {
            amd_finder(atom, True)
        } elif atom is C {
            find O bonded to C
            amd_finder(O, True)
        } elif atom is N {
            find C bonded to N
            amd_finder(C, False)
        }
    }
}