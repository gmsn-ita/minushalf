"""Observe the initial VBM projection without changing MinusHalf's return values."""
import json
from pathlib import Path


def projection_details(projection):
    groups = (("s", 0, 1), ("p", 1, 4), ("d", 4, 9), ("f", 9, 16))
    total = sum(sum(values) for values in projection.values())
    if total <= 0:
        raise ValueError("Initial VBM projection has zero norm")
    return {f"{element}-{orbital}": 100 * sum(values[start:end]) / total
            for element, values in projection.items()
            for orbital, start, end in groups if start < len(values)}


def main():
    from minushalf.utils.band_structure import BandStructure
    from minushalf.minushalf import minushalf
    destination = Path.cwd() / "initial_vbm.json"
    original = BandStructure.vbm_projection
    captured = False

    def observe(self, *args, **kwargs):
        nonlocal captured
        projection = original(self, *args, **kwargs)
        if not captured:
            captured = True
            try:
                record = {"source": "initial machine-readable VBM projection",
                          "raw_projection": projection,
                          "contributions_percent": projection_details(projection)}
                destination.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
            except Exception as exc:
                destination.write_text(json.dumps({"error": str(exc)}) + "\n")
        return projection

    BandStructure.vbm_projection = observe
    minushalf()


if __name__ == "__main__":
    main()
