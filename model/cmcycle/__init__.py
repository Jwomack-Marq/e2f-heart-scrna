"""Cardiomyocyte cell-cycle fate model (``cmcycle``).

Models the four-way outcome of a cardiomyocyte cell cycle -- quiescence,
productive division, binucleation, and nuclear polyploidization -- calibrated
against the two Bergmann-lab FUCCI papers:

  Baniol et al. 2021, Exp Cell Res 408:112880  (FUCCI mouse, in vivo + scRNA-seq)
  Murganti et al. 2022, Front Cardiovasc Med 9:840147 (TNNT2-FUCCI hiPSC-CM)

The two failure modes are kept mechanistically distinct because the data demands
it: binucleation is a completed mitosis with a failed furrow, whereas
polyploidization is G2 arrest in which the anaphase-promoting complex never fires
and the nuclear envelope never breaks down.

Modules
-------
``preflight``  closed-form consistency checks on the calibration set -- run these
               before fitting anything; they already falsified two assumptions.
``baniol``     re-analysis of the P0/P7 scRNA-seq behind every quoted number.
``svg``        dependency-free plotting.
``figures``    regenerates ``model/figures/``.

Stdlib only, deliberately: the whole point of the pre-flight is that it runs in a
bare interpreter, before any solver exists.
"""
from .preflight import Check, Target, load_targets, run_all

__all__ = ["Check", "Target", "load_targets", "run_all"]
