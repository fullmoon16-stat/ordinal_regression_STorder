"""Reproduce Supplementary Table S2 for the original seven categories."""

from pathlib import Path

from white_wine_common import TableSpec, run_table


PROJECT_ROOT = Path(__file__).resolve().parents[1]


if __name__ == "__main__":
    run_table(
        PROJECT_ROOT,
        TableSpec(
            table_name="TableS2",
            n_classes=7,
            response_map={3: 1, 4: 2, 5: 3, 6: 4, 7: 5, 8: 6, 9: 7},
        ),
    )
