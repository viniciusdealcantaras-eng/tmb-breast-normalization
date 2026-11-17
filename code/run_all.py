"""Run the full analysis in order (after 01_extract_cbioportal.py)."""
import runpy
from pathlib import Path

HERE = Path(__file__).parent
for script in ["02_denominator_audit.py", "03_prevalence_cohorts.py", "04_repeated_samples.py",
               "05_gene_null.py", "06_survival.py", "07_continuous_ci.py", "08_subtype.py",
               "09_gene_length.py"]:
    print(f"\n{'=' * 70}\n{script}\n{'=' * 70}")
    runpy.run_path(str(HERE / script), run_name="__main__")
