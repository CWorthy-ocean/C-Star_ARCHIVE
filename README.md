[![codecov](https://codecov.io/github/dafyddstephenson/C-Star/branch/automated_testing/graph/badge.svg?token=Z0L4U76WSG)](https://codecov.io/github/dafyddstephenson/C-Star)
# C-Star
Computational Systems for Tracking Ocean Carbon

## Setup (shell-script-based version)
### First time setup
- Clone this repository and run (from the command line) the script `setup_cstar` in the `setup_cstar` directory, providing your system as an argument (e.g. `./setup_cstar osx_arm64_gnu`). For a list of supported systems, run `setup_cstar --help`.
- The setup script will obtain and compile any external code, and also make an environment on your machine to use when running C-Star in future. To activate this environment, run the command `cstar_env` after completing the setup (note you will have to restart your shell).

### Obtaining C-Star blueprints
- With the C-Star environment active (see above), use the command `cstar_get_blueprint` to obtain and compile a C-Star blueprint. For a list of available blueprint, run `cstar_get_blueprint --help`.
- C-Star blueprints are saved to `${CSTAR_ROOT}/blueprints`
- For help with a specific blueprint, see its README file (e.g. `${CSTAR_ROOT}/blueprints/roms_marbl_example/README.md`)

## Setup (python package)
An early version of the python package implementation of C-Star with specific documentation is available in the `cstar_ocean` directory.

## See Also
- [ROMS-tools](https://github.com/CWorthy-ocean/roms-tools)
