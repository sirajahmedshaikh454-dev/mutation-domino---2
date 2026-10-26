from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

import numpy as np

from app.domain.properties import residue_properties
from app.domain.structure import COMMON_ELEMENTS, ParsedStructure, ResidueRecord


@dataclass(frozen=True)
class ContactRow:
    structure_index: int
    chain_id: str
    sequence_id: int
    insertion_code: str
    residue_id: str
    resname: str
    one_letter: str
    canonical_position: int | None
    distance_angstrom: float
    category: str
    target_atom: str | None = None
    partner_atom: str | None = None
    properties: dict[str, Any] = field(default_factory=dict)

    def to_dict(self) -> dict[str, Any]:
        payload: dict[str, Any] = {
            "structure_index": self.structure_index,
            "chain_id": self.chain_id,
            "sequence_id": self.sequence_id,
            "insertion_code": self.insertion_code,
            "residue_id": self.residue_id,
            "resname": self.resname,
            "one_letter": self.one_letter,
            "canonical_position": self.canonical_position,
            "distance_angstrom": round(self.distance_angstrom, 2),
            "category": self.category,
            "properties": self.properties,
        }
        if self.target_atom is not None:
            payload["target_atom"] = self.target_atom
        if self.partner_atom is not None:
            payload["partner_atom"] = self.partner_atom
        return payload


def _contact_atoms(residue: ResidueRecord) -> list[Any]:
    return [
        atom for atom in residue.atoms if not atom.is_hydrogen and atom.element.upper() in COMMON_ELEMENTS
    ]


def _coords(residue: ResidueRecord) -> np.ndarray | None:
    atoms = _contact_atoms(residue)
    if not atoms:
        return None
    return np.array([[a.x, a.y, a.z] for a in atoms], dtype=np.float64)


def _row(
    residue: ResidueRecord,
    distance: float,
    canonical_position: int | None,
    category: str,
    *,
    target_atom: str | None = None,
    partner_atom: str | None = None,
) -> ContactRow:
    props = residue_properties(residue.one_letter)
    return ContactRow(
        structure_index=residue.index,
        chain_id=residue.chain_id,
        sequence_id=residue.sequence_id,
        insertion_code=residue.insertion_code,
        residue_id=residue.residue_id,
        resname=residue.resname,
        one_letter=residue.one_letter,
        canonical_position=canonical_position,
        distance_angstrom=distance,
        category=category,
        target_atom=target_atom,
        partner_atom=partner_atom,
        properties={
            "category": props.category,
            "charge_class": props.charge_class,
            "aromatic": props.aromatic,
            "hydropathy": props.hydropathy,
        },
    )


def ca_neighbors(
    structure: ParsedStructure,
    target_index: int,
    radius_angstrom: float,
    *,
    structure_to_canonical: dict[int, int] | None = None,
    max_rows: int = 500,
) -> tuple[list[ContactRow], list[str]]:
    target = structure.residue_by_index(target_index)
    warnings: list[str] = []
    if target.ca is None:
        return [], ["The target residue has no CA atom; CA-neighbor analysis is unavailable."]
    center = np.array([target.ca.x, target.ca.y, target.ca.z], dtype=np.float64)
    inverse = structure_to_canonical or {}
    rows: list[ContactRow] = []
    for residue in structure.residues:
        if residue.index == target_index or residue.ca is None:
            continue
        point = np.array([residue.ca.x, residue.ca.y, residue.ca.z], dtype=np.float64)
        distance = float(np.linalg.norm(point - center))
        if distance <= radius_angstrom:
            rows.append(
                _row(
                    residue,
                    distance,
                    _canonical(inverse, residue.index),
                    "ca_neighbor",
                    target_atom="CA",
                    partner_atom="CA",
                )
            )
    rows.sort(key=lambda row: (row.distance_angstrom, row.structure_index))
    if len(rows) > max_rows:
        warnings.append(f"CA-neighbor list truncated to {max_rows} of {len(rows)} rows.")
        rows = rows[:max_rows]
    return rows, warnings


def heavy_atom_contacts(
    structure: ParsedStructure,
    target_index: int,
    radius_angstrom: float,
    *,
    structure_to_canonical: dict[int, int] | None = None,
    max_rows: int = 500,
) -> tuple[list[ContactRow], list[str]]:
    target = structure.residue_by_index(target_index)
    target_atoms = _contact_atoms(target)
    warnings: list[str] = []
    if not target_atoms:
        return [], ["The target residue has no heavy atoms available for contact analysis."]
    center = np.array([[a.x, a.y, a.z] for a in target_atoms], dtype=np.float64)
    inverse = structure_to_canonical or {}
    rows: list[ContactRow] = []
    for residue in structure.residues:
        if residue.index == target_index:
            continue
        partner = _coords(residue)
        if partner is None:
            continue
        atoms = _contact_atoms(residue)
        deltas = center[:, None, :] - partner[None, :, :]
        distances = np.linalg.norm(deltas, axis=2)
        flat_index = int(np.argmin(distances))
        target_atom_index, partner_atom_index = np.unravel_index(flat_index, distances.shape)
        distance = float(distances[target_atom_index, partner_atom_index])
        if distance <= radius_angstrom:
            rows.append(
                _row(
                    residue,
                    distance,
                    _canonical(inverse, residue.index),
                    "heavy_atom_contact",
                    target_atom=target_atoms[target_atom_index].name,
                    partner_atom=atoms[partner_atom_index].name,
                )
            )
    rows.sort(key=lambda row: (row.distance_angstrom, row.structure_index))
    if len(rows) > max_rows:
        warnings.append(f"Heavy-atom contact list truncated to {max_rows} of {len(rows)} rows.")
        rows = rows[:max_rows]
    return rows, warnings


def _canonical(inverse: dict[int, int], structure_index: int) -> int | None:
    value = inverse.get(structure_index)
    return value + 1 if value is not None else None


def out_of_range_residues(structure: ParsedStructure, target_index: int, radius: float) -> list[int]:
    """Return structure indices farther than ``radius`` from the target CA (sanity/debug helper)."""
    target = structure.residue_by_index(target_index)
    if target.ca is None:
        return [r.index for r in structure.residues if r.index != target_index]
    center = np.array([target.ca.x, target.ca.y, target.ca.z])
    far: list[int] = []
    for residue in structure.residues:
        if residue.ca is None:
            continue
        distance = float(np.linalg.norm(np.array([residue.ca.x, residue.ca.y, residue.ca.z]) - center))
        if distance > radius:
            far.append(residue.index)
    return far


def _is_bonded_pair(
    target: ResidueRecord,
    partner: ResidueRecord,
    target_atom: Any,
    partner_atom: Any,
    distance: float,
    inverse: dict[int, int],
) -> bool:
    target_canonical = inverse.get(target.index)
    partner_canonical = inverse.get(partner.index)
    seq_sep = abs(partner.index - target.index)
    canonical_sep: int | None = None
    if target_canonical is not None and partner_canonical is not None:
        canonical_sep = abs(partner_canonical - target_canonical)
    sequence_adjacent = seq_sep == 1 or (canonical_sep is not None and canonical_sep == 1)

    if sequence_adjacent and {target_atom.name, partner_atom.name} == {"C", "N"} and distance <= 2.0:
        return True
    return (
        target.resname == "CYS"
        and partner.resname == "CYS"
        and target_atom.name == "SG"
        and partner_atom.name == "SG"
        and distance <= 2.5
    )


def nonbonded_contacts(
    structure: ParsedStructure,
    target_index: int,
    radius_angstrom: float,
    *,
    structure_to_canonical: dict[int, int] | None = None,
    max_rows: int = 500,
) -> tuple[list[ContactRow], list[str], dict[str, int]]:
    target = structure.residue_by_index(target_index)
    target_atoms = _contact_atoms(target)
    warnings: list[str] = []
    if not target_atoms:
        return (
            [],
            ["The target residue has no heavy atoms available for contact analysis."],
            {
                "eligible_atom_pairs": 0,
                "excluded_bonded_atom_pairs": 0,
            },
        )
    inverse = structure_to_canonical or {}
    target_canonical = inverse.get(target_index)
    rows: list[ContactRow] = []
    excluded_bonded = 0
    eligible_pairs = 0
    for residue in structure.residues:
        if residue.index == target_index:
            continue
        partner_atoms = _contact_atoms(residue)
        if not partner_atoms:
            continue

        seq_sep = abs(residue.index - target_index)
        partner_canonical = inverse.get(residue.index)
        canonical_sep: int | None = None
        if target_canonical is not None and partner_canonical is not None:
            canonical_sep = abs(partner_canonical - target_canonical)
        sequence_adjacent = seq_sep == 1 or (canonical_sep is not None and canonical_sep == 1)

        min_nonbonded_distance: float | None = None
        min_pair: tuple[Any, Any] | None = None
        for target_atom in target_atoms:
            for partner_atom in partner_atoms:
                eligible_pairs += 1
                delta = np.array(
                    [
                        target_atom.x - partner_atom.x,
                        target_atom.y - partner_atom.y,
                        target_atom.z - partner_atom.z,
                    ],
                    dtype=np.float64,
                )
                distance = float(np.linalg.norm(delta))
                if _is_bonded_pair(target, residue, target_atom, partner_atom, distance, inverse):
                    excluded_bonded += 1
                    continue
                if distance <= radius_angstrom and (
                    min_nonbonded_distance is None or distance < min_nonbonded_distance
                ):
                    min_nonbonded_distance = distance
                    min_pair = (target_atom, partner_atom)

        if min_nonbonded_distance is None or min_pair is None:
            continue
        target_atom, partner_atom = min_pair
        row = _row(
            residue,
            min_nonbonded_distance,
            _canonical(inverse, residue.index),
            "nonbonded_contact",
            target_atom=target_atom.name,
            partner_atom=partner_atom.name,
        )
        row.properties["sequence_adjacent"] = sequence_adjacent
        row.properties["sequence_separation"] = seq_sep
        if canonical_sep is not None:
            row.properties["canonical_sequence_separation"] = canonical_sep
        rows.append(row)

    rows.sort(key=lambda row: (row.distance_angstrom, row.structure_index))
    summary = {
        "eligible_atom_pairs": eligible_pairs,
        "excluded_bonded_atom_pairs": excluded_bonded,
    }
    if excluded_bonded:
        warnings.append(
            str(excluded_bonded) + " covalently bonded atom pair(s) were excluded from nonbonded contacts."
        )
    if len(rows) > max_rows:
        warnings.append(
            "Nonbonded contact list truncated to " + str(max_rows) + " of " + str(len(rows)) + " rows."
        )
        rows = rows[:max_rows]
    return rows, warnings, summary
