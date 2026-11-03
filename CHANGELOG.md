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

[Unreleased]: https://github.com/Cao-RQM-Lab/connflex-core/commits/main
