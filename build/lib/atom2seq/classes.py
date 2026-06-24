import math


class Atom:
    """A class representing an atom."""

    def __init__(self, symbol: str, coords: tuple[float]):
        self.symbol = symbol
        self.coords = coords

    def __repr__(self):
        return f"Atom('{self.symbol}', {self.coords})"

    def __eq__(self, other):
        # Returns True if and only if both the symbol and the coordinates are
        # the same
        return (self.symbol == other.symbol) and (self.coords == other.coords)


# Defining this class before putting anything in it so that I can use it for
# type annotations.
class Mol:
    pass


class Mol:
    """A class representing a molecule. Supports checking equality."""

    def __init__(self, atoms: list[Atom], bonds: list[list[int]]):
        self._atoms = atoms
        self._bonds = bonds
        self._backbone = [False for _ in self._atoms]
        self._aas = [False for _ in self._atoms]
        self._n_term = False

    def __repr__(self):
        return f"Mol({self._atoms}, {self._bonds})"

    def __eq__(self, other):
        # Fast tracks a False return if the lengths are not the same.
        if (len(self._atoms) != len(other.get_atoms())) or (
            len(self._bonds) != len(other.get_bonds())
        ):
            return False
        # If the lengths are the same, checks entrywise.
        else:
            for i in range(len(self._atoms)):
                if self._atoms[i] != other.get_atoms()[i]:
                    return False
            for j in range(len(self._bonds)):
                if self._bonds[j] != other.get_bonds()[j]:
                    return False
        return True

    def is_bond(self, n: int, m: int) -> bool:
        """Detects whether the atoms at indices n and m are bonded."""
        # Loops over every bond and checks whether that bond is between the two
        # given atoms.
        i = 0
        while i < len(self._bonds):
            bond = self._bonds[i]
            if (bond[0] == n) and (bond[1] == m):
                return True
            elif (bond[1] == n) and (bond[0] == m):
                return True
            i += 1
        return False

    def del_bond(self, n: int, m: int) -> None:
        """Deletes a bond between the two given atoms."""
        # Goes over each bond and checks if it is between the two given atoms.
        # Removes it if it is.
        for bond in self._bonds:
            if (bond[0] == n) and (bond[1] == m):
                self._bonds.remove(bond)
            elif (bond[1] == n) and (bond[0] == m):
                self._bonds.remove(bond)

    def add_bond(self, n: int, m: int) -> None:
        """Adds a bond between the two given atoms."""
        if n == m:
            raise ValueError("You cannot bond an atom with itself.")
        if self.is_bond(n, m):
            # Raises an AttributeError if there is already a bond between the
            # given atoms.
            raise AttributeError("You cannot bond two bonded atoms.")
        else:
            self._bonds.append([min(n, m), max(n, m)])

    def dist(self, n: int, m: int) -> float:
        """Calculates the Euclidean distance between the given atoms."""
        n_coords = self._atoms[n].coords
        m_coords = self._atoms[m].coords
        return math.sqrt(
            (n_coords[0] - m_coords[0]) ** 2
            + (n_coords[1] - m_coords[1]) ** 2
            + (n_coords[2] - m_coords[2]) ** 2
        )

    def get_bonds(self) -> list[list[int]]:
        """Returns the list of bonds."""
        return self._bonds

    def get_atoms(self) -> list[Atom]:
        """Returns the list of atoms."""
        return self._atoms

    def del_atom(self, idx: int) -> None:
        """Deletes the given atom from the list of atoms, also updating the
        lists of bonds, sidechain numbering, and backbone list."""
        print("deleting ", self._atoms[idx])
        # Removes the given index from the three lists that use this indexing
        # system.
        self._atoms.pop(idx)
        self._aas.pop(idx)
        self._backbone.pop(idx)
        # Creates a new bonds list, then adds each bond to it if and only if it
        # does not involve the given atom.
        new_bonds = []
        for i in range(len(self._bonds)):
            bond = self._bonds[i]
            if (bond[0] != idx) and (bond[1] != idx):
                new_bonds.append(bond)
        # Updates the indices of the remaining bonds to reflect the fact that
        # an index was removed.
        for bond in new_bonds:
            if bond[0] > idx:
                bond[0] -= 1
            if bond[1] > idx:
                bond[1] -= 1
        self._bonds = new_bonds

    def get_bonded(self, idx: int) -> list[int]:
        """Returns a list of the indices of atoms bonded to the atom at index
        idx."""
        out = []
        for bond in self._bonds:
            if bond[0] == idx:
                out.append(bond[1])
            elif bond[1] == idx:
                out.append(bond[0])
        return out

    def find_submol(self, idx: int) -> list[int]:
        """Returns a list of all indices connected to the atom at index idx."""
        # Initializes the current indices and the list to be returned.
        current_idcs = {idx}
        print(f"{current_idcs=}")
        out = {idx}
        print(f"{out=}")
        done = False
        while not done:
            for i in current_idcs:
                for elt in self.get_bonded(i):
                    out.add(elt)
                    print(f"{out=}")
                if out == current_idcs:
                    done = True
            for elt in out:
                current_idcs.add(elt)
                print(f"{current_idcs=}")
        return list(out)

    def del_submol(self, idx: int) -> None:
        """Deletes all atoms connected to the atom at index idx."""
        submol_idcs = self.find_submol(idx)
        submol_idcs = sorted(submol_idcs, reverse=True)
        for idx in submol_idcs:
            self.del_atom(idx)

    def split_submol(self, idx: int) -> Mol:
        """Returns a Mol containing all atoms connected to the atom at index
        idx."""
        submol_idcs = self.find_submol(idx)
        atoms = []
        bonds = []
        for idx in submol_idcs:
            atoms.append(self._atoms[idx])
        for bond in self._bonds:
            if bond[0] in submol_idcs:
                new_bond = [
                    submol_idcs.index(bond[0]),
                    submol_idcs.index(bond[1]),
                ]  # noqa
                bonds.append(new_bond)
        self.del_submol(idx)
        return Mol(atoms, bonds)

    def set_n_term(self, idx: int) -> None:
        """Sets the N-terminus of the molecule to the atom at index idx."""
        self._n_term = idx

    def get_n_term(self) -> int | bool:
        """Returns the N-terminus of the molecule. Returns False if it is not
        set."""
        return self._n_term

    def set_backbone(self, idx_list: list[int]) -> None:
        """Sets the backbone to the given indices."""
        self._backbone = [False for _ in self._atoms]
        for idx in idx_list:
            self._backbone[idx] = True

    def number_aas(self, idx_list: list[int], num_list: list[int]):
        """Numbers the atoms at the indices in idx_list with the numbers in
        num_list."""
        for i in range(len(idx_list)):
            idx, num = idx_list[i], num_list[i]
            self._aas[idx] = num

    def get_backbone(self) -> list[bool]:
        """Returns a list of booleans corresponding to the indices of atoms in
        the molecule, where it is True if that atom is in the backbone."""
        return self._backbone

    def get_aas(self) -> list[int | bool]:
        """Returns a list of ints and booleans corresponding to the indices of
        atoms in the molecule, where it is False if it is not the beginning of
        a sidechain and numbered if it is."""
        return self._aas
