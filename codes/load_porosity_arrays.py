from pathlib import Path
import numpy as np

def load_porosity_arrays(
    input_dir,
    x: int = 0,                          # changeable
    y_range = range(1, 6),               # 1..5 inclusive
    kinds = ("impactor", "target"),
    timestep_range = range(0, 100),      # 0..99 inclusive
    has_header: bool = False             # set True if your CSVs have a header row
) -> dict[tuple[int, int, str, int], np.ndarray]:
    """
    Loads CSVs named upright_porosity_{x}_case_{y}_{kind}_{timestep}.csv from input_dir
    using NumPy (no pandas). Missing or empty files are skipped.

    Returns a dict keyed by (x, y, kind, timestep) -> np.ndarray
    """
    input_dir = Path(input_dir)
    results: dict[tuple[int, int, str, int], np.ndarray] = {}

    for y in y_range:
        for kind in kinds:
            for t in timestep_range:
                fname = f"upright_porosity_{x}_case_{y}_{kind}_{t}.csv"
                fpath = input_dir / fname

                # Skip if missing or empty file (“if empty skip”)
                if not fpath.is_file() or fpath.stat().st_size == 0:
                    continue

                try:
                    arr = np.genfromtxt(
                        fpath,
                        delimiter=",",
                        dtype=float,
                        autostrip=True,
                        skip_header=1 if has_header else 0
                    )
                    # Ensure 2D shape even for single-row files
                    arr = np.atleast_2d(arr)

                    # If load produced an empty array for some reason, skip
                    if arr.size == 0:
                        continue

                    results[(x, y, kind, t)] = arr

                except Exception as e:
                    # unreadable → skip (you can print or log if you want)
                    # print(f"[warn] Failed to read {fpath}: {e}")
                    pass

    return results

# ---------------------------
# Example usage:
# ---------------------------
# data = load_porosity_arrays("path/to/input_folder", x=0, has_header=False)
# print(f"Loaded {len(data)} files")
# arr = data.get((0, 1, "impactor", 0))  # Example access
