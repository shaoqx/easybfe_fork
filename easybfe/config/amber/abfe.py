from __future__ import annotations
from typing import Any, Optional
from pathlib import Path
from pydantic import BaseModel, Field, model_validator
from .simulation import AmberFepSimulationConfig


class BoreschRestraintGeneratorConfig(BaseModel):
    algorithm: str = 'rxrx'
    rst_wts: tuple[float, float, float, float, float, float] = Field(
        default=(10.0, 10.0, 10.0, 10.0, 10.0, 10.0))
    options: dict[str, Any] = Field(default_factory=dict)


class AmberAbfeConfig(BaseModel):
    protein: Optional[Path] = None
    ligand: Optional[Path] = None
    ligand_batch: Optional[list[Path]] = None
    output_dir: Optional[Path] = None
    ligand_base: Optional[Path] = None
    output_base: Optional[Path] = None

    # Keep default ABFE behavior; users can set boresch: null to disable.
    boresch: Optional[BoreschRestraintGeneratorConfig] = Field(
        default_factory=BoreschRestraintGeneratorConfig
    )

    # Legs are optional to support single-leg reorg workflows.
    complex: Optional[AmberFepSimulationConfig] = None
    solvent: Optional[AmberFepSimulationConfig] = None
    restraint: Optional[AmberFepSimulationConfig] = None

    @property
    def active_legs(self) -> list[str]:
        return [leg for leg in ('complex', 'solvent', 'restraint') if getattr(self, leg) is not None]

    @model_validator(mode='after')
    def validate_legs(self) -> 'AmberAbfeConfig':
        legs = self.active_legs
        if len(legs) == 0:
            raise ValueError('At least one simulation leg must be configured in AmberAbfeConfig')
        if self.restraint is not None and self.boresch is None:
            raise ValueError('restraint leg requires boresch restraints; set boresch config or disable restraint leg')
        return self
