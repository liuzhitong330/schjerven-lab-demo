# Early IKAROS-response readouts

An assistant-supported, reproducible application demo prepared by Cathy Liu, inspired by the Schjerven Lab at UCSF. Not affiliated with or endorsed by the lab.

**Question:** Which early transcript responses remain consistent across two deposited B-ALL models, and which need a model/preparation-specific check before selecting a follow-up measurement?

The static explorer uses all 60,624 unique gene IDs from eight GSE310958 libraries. It shows within-model condition-mean effects, the four unpaired sample-comparison ratios, individual measurements, normalization sensitivity, and a downloadable proposed readout list. It also includes a provenance manifest and an explicitly unperformed follow-up recording template.

- [Live demo](https://liuzhitong330.github.io/schjerven-lab-demo/)
- [Reproduction instructions and interpretation limits](analysis/README.md)
- [Source hashes and provenance](analysis/provenance.json)
- [Numerical checks](analysis/checks.json)
- [Expected numerical results](analysis/summary.json)
- [Third-party attribution and licenses](analysis/THIRD_PARTY_NOTICES.md)

## Run locally

No build or external service is needed. Open `index.html`, or run `python -m http.server 8000` from this directory. The browser does not upload data or contact an API.

## Reproduce the data

```sh
python -m pip install -r analysis/requirements.txt
python analysis/build_data.py --inputs /path/to/public-input-cache --download
```

All six source inputs are hash-checked. The data builder tests median-ratio size factors and the attributed Pearson-correlation component ported from the Schjerven-coauthored ssvQC project. It does not run the entire R package or reproduce DESeq2 significance.

## Evidence boundaries

Two biological replicates per condition per model, one 24-hour time point. PDX2 is historical data reprocessed for this study, not an independent validation cohort. Model and library preparation are confounded. Compare within-model changes, not cross-model baseline expression. Four sample ratios are not four independent replicates or a confidence interval. Candidate counts depend on exploratory thresholds, not significance tests. No pseudocount is added to response ratios. RNA does not establish protein, direct regulation, therapeutic vulnerability or clinical utility.

ERG is an established source-study result. The added output is a transparent way to qualify readouts and controls for independent follow-up, not a new dependency claim.

Contact: Cathy Liu · cathyliu014@gmail.com · 415-216-3799
