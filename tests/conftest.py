from collections.abc import Callable

import pytest
from shapely import Polygon

from connflex_core import ExtractionParams, Layer, Layout, Material, Stack

UM = 1e-6


@pytest.fixture
def stack() -> Stack:
    return Stack(
        device=Layer(thickness=25 * UM, material=Material.silicon()),
        box=Layer(thickness=2 * UM, material=Material.oxide()),
        handle=Layer(thickness=400 * UM, material=Material.silicon()),
    )


@pytest.fixture
def params() -> Callable[..., ExtractionParams]:
    def make(**overrides: object) -> ExtractionParams:
        values: dict[str, object] = {
            "release_hole_max_area": 100 * UM * UM,
            "release_distance": 5 * UM,
            "body_min_width": 20 * UM,
        }
        values.update(overrides)
        return ExtractionParams(**values)  # type: ignore[arg-type]

    return make


@pytest.fixture
def layout_um() -> Callable[..., Layout]:
    def make(
        device: list[Polygon],
        handle: list[Polygon],
        release: list[Polygon] | None = None,
        vias: list[Polygon] | None = None,
    ) -> Layout:
        return Layout.from_shapely(
            device=device, handle=handle, release=release, vias=vias, unit=UM
        )

    return make
