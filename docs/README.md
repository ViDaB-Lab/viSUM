# viSUM user documentation

Start with [installation and usage](usage.md), then use the reference pages as
needed. The [README](../README.md) explains the workflow and what's new.

| Task | Guide |
| --- | --- |
| Install; run one sample or a batch; resume | [Usage](usage.md) |
| Download/build or reuse databases | [Database setup](database-setup.md) |
| Prepare viCAT viral/nonviral references | [viCAT databases](vicat-database-preparation.md) |
| Select or update a VITAP database | [VITAP databases](vitap-database-lifecycle.md) |
| Set CPU/RAM budgets or use SLURM | [Resources and execution](resources.md) |
| Understand viral/TE candidates and taxonomy | [Interpreting results](interpretation.md) |
| Find output files and column definitions | [Reporting bundle](reporting-bundle.md) |
| Understand TE flags | [TE evidence](te-candidate-routing.md) |
| Interpret undetermined CheckV quality | [CheckV context](checkv-evidence-context.md) |
| Understand extracted proviral regions | [Provirus boundaries](provirus-boundary-policy.md) |
| Diagnose failed or unexpected runs | [Troubleshooting](troubleshooting.md) |

Options not commonly needed are documented beside their defaults in
`visum.config`. Display command-line help with:

```bash
nextflow run visum_nextflow.nf --help
```

This branch contains runtime code and user documentation. Test suites and
development records are not required for installation and are not shipped here.
