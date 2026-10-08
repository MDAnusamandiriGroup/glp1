# Data dictionary and eligibility mapping

## NHANES files and variables

| Construct | File | Variable(s) | Operational rule |
|---|---|---|---|
| Age | `DEMO_L` | `RIDAGEYR` | >=45 years |
| Survey design | `DEMO_L` | `SDMVPSU`, `SDMVSTRA` | NHANES complex design |
| Analysis weight | `GHB_L` / `DEMO_L` | `WTPH2YR` | Phlebotomy two-year weight because HbA1c is required |
| BMI | `BMX_L` | `BMXBMI` | >=27 kg/m2 |
| Diagnosed diabetes | `DIQ_L` | `DIQ010` | Exclude value 1 (yes) |
| Biochemical diabetes | `GHB_L` | `LBXGH` | Exclude HbA1c >=6.5% |
| Prior myocardial infarction | `MCQ_L` | `MCQ160E` | value 1 (yes) |
| Prior stroke | `MCQ_L` | `MCQ160F` | value 1 (yes) |

## Mapping to SELECT

The public-data operationalization requires age >=45, BMI >=27, no diagnosed
diabetes, HbA1c <6.5%, and self-reported prior MI or stroke. SELECT also allowed
symptomatic peripheral arterial disease and applied additional exclusions
(recent cardiovascular/neurologic event, NYHA class IV heart failure, end-stage
kidney disease/dialysis, recent glucose-lowering therapy, and planned
revascularization). Those features cannot all be reconstructed reliably in the
selected public NHANES files.

Therefore:

- use NHANES for the national scale and age distribution only;
- standardize the main dynamic results per 100,000 eligible initiators;
- do not call the NHANES denominator an exact SELECT-eligible count;
- include a sensitivity multiplier around the NHANES count if projecting an
  absolute US total.

## Outcome definitions

- `first_mace`: first cardiovascular death/undetermined death, nonfatal MI, or
  nonfatal stroke.
- `fatal_or_undetermined_first_mace`: the fatal/undetermined component when it
  is the first primary event.
- `D_ncv`: non-cardiovascular death before a first MACE, included as a competing
  event.
- All-cause death is deliberately not claimed as a model output in version 1.

