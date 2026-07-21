## Test environments

* local macOS install, R 4.6.1
* win-builder (devel)

## R CMD check results

0 errors | 0 warnings | 1 note

The one NOTE is the expected archival notice:

```
* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Nicholas Tierney <nicholas.tierney@gmail.com>'
New submission
Package was archived on CRAN
CRAN repository db overrides:
  X-CRAN-Comment: Archived on 2025-09-20 as requires archived package
    'greta'.
```

## Submission notes

This is a resubmission of a previously archived package.

`greta.gp` was archived from CRAN on 2025-09-20 as a consequence of its
dependency, `greta`, being archived. `greta` has since been fixed and returned
to CRAN as version 0.6.0. This release of `greta.gp` depends on
`greta` (>= 0.6.0) and has been checked against it, so the reason for the
archival no longer applies.

The tests and examples require Python, TensorFlow and TensorFlow Probability,
which are not available on CRAN's check machines. As with the `greta` package
itself, examples are wrapped in `\dontrun{}` and the tests skip themselves when
those dependencies are absent, so the checks remain fast and self-contained on
CRAN. The full test suite and all examples have been run locally against
`greta` 0.6.0 and pass.

## revdepcheck results

As the package was archived, there are no reverse dependencies to check.

## Method references

The methods implemented in this package are described in the reference already
given in the `Description` field, Golding (2019) <doi:10.21105/joss.01601>,
which describes the `greta` software that this package extends.
