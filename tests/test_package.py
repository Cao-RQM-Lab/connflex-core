import connflex_core


def test_installed_package_reports_a_version() -> None:
    assert connflex_core.__version__
