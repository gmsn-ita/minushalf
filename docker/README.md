# Quantum ESPRESSO test environment

This container lets us test MinusHalf with Quantum ESPRESSO in GitHub Actions.
It includes Python 3.14, QE 7.5 and Open MPI, with the four executables used by
the workflow: `pw.x`, `projwfc.x`, `ld1.x` and `virtual_v2.x`.

MinusHalf is installed from the source checked out for each CI run, so the
tests use the code being reviewed rather than a copy stored in the image.

## How the tests run

The **QE CI** workflow builds the container from `docker/Dockerfile` and uses
cached build layers when available. The first build takes longer because it
needs to compile QE. The tests do not depend on downloading our published
image from GHCR.

Inside the container, the workflow:

1. Installs MinusHalf and its test dependencies.
2. Runs the Python tests with `pytest`.
3. Runs a silicon SCF calculation with two MPI processes and checks convergence.
4. Runs `minushalf execute` and checks that it produces a nonempty results file.

The silicon inputs and pseudopotential are in `tests/integration/qe/si/`.
A passing run confirms that these steps completed successfully. Checking the
physical gap and whether the correction was applied correctly is separate
work; those checks are not part of this CI yet.

## Running and checking a test

Tests start when you push to `qe_container` or open or update a pull request
targeting `qe_execute` or `main`, provided the workflow is included. The workflow
also supports manual dispatch. For pull requests, GitHub tests the proposed
merge with the base branch.

Open **Actions > QE CI** to follow a run. Its downloadable artifact contains
the available logs, test reports, silicon results and version information.
The upload step also runs after test failures to help with debugging.

Each job has a 45-minute limit, including the image build and installation.
This limit applies only to GitHub Actions, not to cluster calculations.

## Publishing the base image

The **Build QE Docker Image** workflow publishes the QE environment to
`ghcr.io/gmsn-ita/minushalf-qe-base`. It runs when the Dockerfile or the
build workflow changes on `qe_container`, and also supports manual dispatch.

The image is tagged with `qe-7.5` and `qe-7.5-<commit>`. Publishing uses
the token supplied by GitHub Actions and requires package write access.

The image contains QE and the build tools. MinusHalf is installed from
the checked-out source when the tests run.
