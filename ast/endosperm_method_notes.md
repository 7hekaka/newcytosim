# Endosperm Actin Modeling Notes

Context from the two methods sections:

- `PIIS1534580722002520.pdf` is the Cytosim paper. It uses a 3D stochastic agent-based endocytic actin model based on Akamatsu et al. 2020.
- Its useful physics for us is not vesicle geometry itself, but the separation between assembly sites, filament anchoring, and force readout. Filament anchoring points convert polymerization into pulling/squeezing forces; not every filament has to be directly bound to the object being moved.
- It also emphasizes robustness to network geometry. Multiple founding/mother filaments and local geometric constraints can still produce usable force even when exact filament organization varies.
- `s41598-024-69422-3.pdf` is not Cytosim. It is a continuum finite-element chemo-transport-mechanics model. It treats actin force generation as localized conversion of G-actin to an F-actin network, with cytosolic drag, network elasticity, and disassembly represented by an average filament lifetime.

Physics we may be missing or simplifying too aggressively:

- Local activation/assembly zones. Our free-filament tests should not be uniformly random everywhere if the biological claim is "asters in the nuclei's midst"; they should include perinuclear assembly zones without hard aster tethering.
- Active turnover during movement. Treadmilling is present, but we should explicitly test cutter/chewer-mediated remodeling in the division cases.
- Force-dependent polymerization. Current ACTIN8 uses `growing_force = inf`, so polymerization speed ignores load.
- Steric/crowding and mesh mechanics. Small smoke tests disable global sterics; production runs may need a steric-on comparison because crowding and mesh size affect force transmission.
- Crosslink remodeling. Static crosslinker counts may miss the gel-like elastic response that lets actin transmit force without every filament being anchored.
- Nuclear division should be modeled as a staged event: two nuclei become four, with actin assembly continuing through the event rather than ACTIN8 being killed and ACTIN10 restarted.

Immediate TODOs:

- Watch `05_division_bound_asters`, `06_division_free_filaments`, and `07_division_perinuclear_free_turnover` side by side.
- Check whether perinuclear unbound assembly can move or space nuclei without the clumping artifact from killing ACTIN8.
- Add a steric-on version only after the smoke tests show a promising qualitative behavior.
- For ACTIN8-to-ACTIN10 production sweeps, keep one fiber class and change polymerization speed instead of renucleating ACTIN10.
- Build analysis around observable links to coenocytic endosperm: nuclear spacing, pairwise distance change after division, actin density around nuclei, filament orientation relative to nuclei, and clump/bridge persistence.
