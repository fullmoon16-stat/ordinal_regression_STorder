"""Reproduce Supplementary Table S3 for the five-category response."""

from pathlib import Path

from white_wine_common import TableSpec, run_table


PROJECT_ROOT = Path(__file__).resolve().parents[1]


if __name__ == "__main__":
    run_table(
        PROJECT_ROOT,
        TableSpec(
            table_name="TableS3",
            n_classes=5,
            response_map={3: 1, 4: 1, 5: 2, 6: 3, 7: 4, 8: 5, 9: 5},
        ),
    )
