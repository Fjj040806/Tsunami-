# Laboratory data provenance and decisions

Official benchmark page: https://nctr.pmel.noaa.gov/benchmark/Solitary_wave/

Archive: https://nctr.pmel.noaa.gov/benchmark/Solitary_wave/tsunami3_runup.zip

Pinned SHA-256: `b1eb73cbca75d7bb1430374e5a94aaf0abba0ce7a57b3e23fe9e67da943a85cb`

Fetch with `python scripts/fetch_noaa_data.py`. The fetcher rejects a different archive hash. The raw archive is stored in ignored `data/raw/` and is excluded from the delivered repository ZIP. Derived metrics, predictions and attributed comparison figures are included. The project MIT license does not confer rights over third-party data or publications.

## Audit and fixed choices

The archive Case A contains 600 rows at 0.05 s spacing (265.05–295.00 s), with time plus seven gauge columns. Water-surface values are in metres. The selected comparison ends at 283.00 s: 360 samples. First-20-sample means are subtracted independently from each gauge; there is no fitted trend, time shift or fitted model parameter.

Archive gauge positions are G4=12.64, G5=15.04, G6=17.22, G7=19.40, G8=20.86, G9=22.33 and G10=22.80 m from the wavemaker. The web page lists G4=12.22 m and nominal Case A H/d=0.038, whereas the archive describes a roughly 0.039 case. These discrepancies are unresolved. The delivered comparison prefers archive metadata because it describes the exact downloaded records; it does not silently reconcile the sources.

The official figure's shore-origin junction coordinates (8.181, 3.821, 0.891 m) are converted to wavemaker-origin coordinates using the wall position 23.22 m. The resulting bed has junctions 15.039, 19.399 and 22.329 m; slopes are 1:53, 1:150 and 1:13. Still-water depth is 0.218 m upstream and approximately 0.047664 m at the truncated reflecting wall. This is a fully wet truncated-domain reconstruction, not a wet/dry runup computation.

G4 is used as an incoming signal at the left boundary. Only G5–G10 are comparison gauges. A measured total surface trace is treated as incoming before substantial reflected waves return to G4. The chosen end time precedes the estimated return of the principal reflected peak, but reflected tails and the model's characteristic boundary approximation remain limitations. There is no independence claim for the G4 input itself.

The paper therefore labels this a **conditional laboratory comparison**, not an official benchmark certification or historical-tsunami validation. All numerical decisions are fixed in `configs/final_studies.json` and `scripts/run_final_studies.py`.
