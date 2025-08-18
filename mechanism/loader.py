import cantera as ct
from typing import List, Dict

class Mechanism:
    """Wrapper around Cantera Solution for mechanism editing."""

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.solution = ct.Solution(file_path)

    @property
    def species_names(self) -> List[str]:
        return [s.name for s in self.solution.species()]

    def molecular_weights(self) -> Dict[str, float]:
        return {s.name: s.molecular_weight for s in self.solution.species()}

    def reactions(self) -> List[ct.Reaction]:
        return list(self.solution.reactions())

    def remove_species(self, remove: list[str]) -> None:
        """Remove species and all reactions involving them."""
        # Get allowed species set
        remaining_species = [s for s in self.solution.species() if s.name not in remove]

        # Filter reactions that only involve retained species
        allowed_names = {s.name for s in remaining_species}
        valid_reactions = []
        for r in self.solution.reactions():
            reactants = set(r.reactants.keys())
            products = set(r.products.keys())
            if reactants.issubset(allowed_names) and products.issubset(allowed_names):
                valid_reactions.append(r)

        # Rebuild solution
        self.solution = ct.Solution(
            thermo='IdealGas',
            kinetics='GasKinetics',
            species=remaining_species,
            reactions=valid_reactions
        )
        assert set(self.species_names) == allowed_names

    def save(self, out_path: str):
        """Save the current reduced mechanism to a YAML file."""
        from cantera import Species
        Species.save_yaml(
            self.solution.species(),
            out_path,
            reactions=self.solution.reactions(),
            transport=False,
            thermo="IdealGas",
            kinetics="GasKinetics"
        )
