# Tier calibration tables

Rendered by `run_tier_calibration.py` from `TIER_CALIBRATION.json`.

## All four stages

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | 21 / 21 | 1.00 (0.84-1.00) | 0 | 0 |
| `divergent` | 0 / 5 | — | 0 | 0 |
| `low-depth` | 0 / 10 | — | 0 | 9 |
| `truncated` | 0 / 11 | — | 0 | 11 |
| `multi-copy` | 4 / 5 | 1.00 (0.40-1.00) | 0 | 0 |
| `degraded` | 0 / 10 | — | 0 | 8 |
| **all** | 25 / 62 | 1.00 (0.86-1.00) | 0 | 28 |

## Stage A

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | 6 / 6 | 1.00 (0.54-1.00) | 0 | 0 |
| `truncated` | 0 / 11 | — | 0 | 11 |
| `degraded` | 0 / 10 | — | 0 | 8 |
| **all** | 6 / 27 | 1.00 (0.54-1.00) | 0 | 19 |

## Stage B

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | 5 / 5 | 1.00 (0.48-1.00) | 0 | 0 |
| `multi-copy` | 4 / 5 | 1.00 (0.40-1.00) | 0 | 0 |
| **all** | 9 / 10 | 1.00 (0.66-1.00) | 0 | 0 |

## Stage C

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `supported` | 10 / 10 | 1.00 (0.69-1.00) | 0 | 0 |
| `divergent` | 0 / 5 | — | 0 | 0 |
| **all** | 10 / 15 | 1.00 (0.69-1.00) | 0 | 0 |

## Stage D

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `low-depth` | 0 / 10 | — | 0 | 9 |
| **all** | 0 / 10 | — | 0 | 9 |

## Stage D, --isolate control

| case class | reached the top tier | top-tier precision (95% CI) | false accepts | withheld but exact |
|---|---|---|---|---|
| `low-depth` | 0 / 10 | — | 0 | 9 |
| **all** | 0 / 10 | — | 0 | 9 |
