"""Entry point of a worker process: ``python -m connflex_core._serve``."""

from connflex_core.workers import serve

if __name__ == "__main__":
    serve()
