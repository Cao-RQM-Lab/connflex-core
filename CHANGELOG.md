# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project uses [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Packaging, code standards, CI, release workflow, and documentation.
- `Material`, `Layer`, `Stack`, `ExtractionParams`, and `SolveOptions` with
  validation, and `load_config` for TOML process files.
- `Layout` with shapely and GDS readers and an adapter registry.
- `extract` producing a `Structure` of rigid bodies and beams: layer healing,
  release-hole and oxide estimate, anchor selection, pruning, body assembly
  through the stack, beam qualification, a networkx graph, and a JSON tree.
- Mass, center of mass, and inertia tensor of every body by triangulated
  quadrature over the stack, with the 6 by 6 body mass matrix. Release openings are
  removed from the device-layer mass.
- Beam stiffness: torsion constant with the short-beam correction, the local
  12 by 12 Euler-Bernoulli matrix, and the transformation onto body degrees of
  freedom.
- Assembly of the 6N by 6N mass and stiffness matrices, modal analysis with
  dense and sparse solvers, `ModalResult` with `.npz` export, and a static solve.

[Unreleased]: https://github.com/Cao-RQM-Lab/connflex-core/commits/main
