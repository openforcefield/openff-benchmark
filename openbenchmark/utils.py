from pathlib import Path

HERE = Path(__file__).parent


def datapath(relative_path):
    return HERE / "data" / relative_path
