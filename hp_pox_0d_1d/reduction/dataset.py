from pathlib import Path
import json
import pandas as pd
import numpy as np


def prepare_dataset(cases_root: str, cases: list[str], outdir: Path) -> int:
    outdir.mkdir(parents=True, exist_ok=True)
    manifest = {
        'cases_root': str(Path(cases_root).resolve()),
        'cases': cases,
        'files': [],
        'taps_pct': [5, 20, 50, 100]
    }

    for case in cases:
        csv_path = Path(cases_root) / case / 'axial_profiles.csv'
        kpi_path = Path(cases_root) / case / 'kpis.json'
        if not csv_path.exists():
            print(f"Missing axial CSV for {case}: {csv_path}")
            continue
        df = pd.read_csv(csv_path)
        # Minimal features: z, T_C, X_H2, X_CO, X_CO2, X_CH4
        cols = [c for c in df.columns if c in ['z_m','T_C','X_H2','X_CO','X_CO2','X_CH4']]
        feats = df[cols].copy()
        # Compute taps at 5/20/50/100% length
        if 'z_m' in feats.columns:
            z = feats['z_m'].values
            z_max = float(np.nanmax(z)) if len(z) else 0.0
            taps = {}
            for pct in manifest['taps_pct']:
                z_target = z_max * (pct/100.0)
                idx = int(np.argmin(np.abs(z - z_target))) if len(z) else 0
                row = feats.iloc[idx].to_dict()
                taps[str(pct)] = {k: (float(v) if pd.notna(v) else None) for k, v in row.items()}
        else:
            taps = {}

        f_out = outdir / f'{Path(case).name}_features.parquet'
        feats.to_parquet(f_out, index=False)
        entry = {'case': case, 'features': str(f_out.resolve()), 'taps': taps}
        manifest['files'].append(entry)
        # Copy KPIs path
        if kpi_path.exists():
            entry['kpis'] = str(kpi_path.resolve())

    with open(outdir / 'manifest.json', 'w') as f:
        json.dump(manifest, f, indent=2)
    print(f"Dataset prepared at: {outdir.resolve().as_posix()} | {len(manifest['files'])} cases")
    return 0


