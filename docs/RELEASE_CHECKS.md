# Final release checks

The extracted repository passed all 35 tests (0 failures, 0 errors; 4.125 s). All recorded numerical-source hashes for the final studies match the delivered sources. The ZIP passed CRC and per-file SHA-256 checks and excludes third-party raw laboratory observations.

The final 16-page PDF was rendered and visually inspected on every page; no overfull boxes or unresolved references remain. Five notebook code cells executed in order in a shared namespace, including the independent lake-at-rest assertion; five PNG outputs were captured. The animation contains 61 frames; selected frames were visually inspected.

Offline explorer JavaScript and data passed functional checks for initial draw, scenario selection, slider, play/pause and all finite frame arrays using a mocked DOM/canvas. A headless-browser rendering check was unavailable in the build environment; no full browser-rendering verification is claimed.

`results/release_checks.json` records these checks. `results/release_manifest.json` identifies delivered file bytes; it excludes itself.
