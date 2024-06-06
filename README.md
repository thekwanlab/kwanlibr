# kwanlibr
R library for commonly used functions in the Kwan Lab.

# Development

See more details in the Kwan Lab Bioinformatics google drive

1. Obtain the package repository with git and start a new branch
1. Make changes
1. Write documentation under the `roxygen2` standard
1. Test loading the package with `devtools::load_all()`
1. Run the test suite corresponding with your changes
    1. This will require obtaining a `informal_test/test_suite_config.yaml` file containing paths to personal files used for testing
1. Update package metadata (such as imports, description, version number) with `usethis`
1. Commit, push to non-main branch, and make pull request