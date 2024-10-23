Depmap Omics Long Read RNA
---

This repository contains WDL workflows and a Python module to continuously run workflows on new RNA long read samples and sync BAM files between Terra and Gumbo.

# Installation

1. Install the required system dependencies:

    - [pyenv](https://github.com/pyenv/pyenv)
    - [Poetry](https://python-poetry.org/)
    - [pre-commit](https://pre-commit.com/)

2. Install the required Python version (3.12.3):

   ```bash
   pyenv install "$(cat .python-version)"
   ```

3. Confirm that `python` maps to the correct version:

   ```
   python --version
   ```

4. Set the Poetry interpreter and install the Python dependencies:

   ```bash
   poetry env use "$(pyenv which python)"
   poetry install
   ```

5. Copy `env.dist` into a new file `.env` and fill it out:

    - `WOMTOOL_JAR` is the path to the .jar file downloaded from [https://github.com/broadinstitute/cromwell/releases](broadinstitute/cromwell)
    - `HASURA_ADMIN_SECRET` is a list of Terra users (email addresses or group names) that should be considered owners of workflow method configs

## Credentials

Your GCP `DEFAULT_APPLICATION_CREDENTIALS` should already by configured in order to run commands.

# Development

Run `pre-commit run --all-files` to automatically format your code with [Ruff](https://docs.astral.sh/ruff/) and check static types with [Pyright](https://microsoft.github.io/pyright).

Whenever possible, function/method arguments and return values should be validated with Pydantic or Pandera (if a data frame).

## GraphQL code generation

This repo uses [ariadne-codegen](https://github.com/mirumee/ariadne-codegen) to generate the `gumbo_gql_client` module. It uses the folder of GraphQL queries (`./gql`) and the current GraphQL schema for a particular Gumbo environment to generate all of the Python classes, Pydantic models, and query/mutation methods for interacting with the [Gumbo GraphQL Service](https://github.com/broadinstitute/gumbo_client/tree/main/gumbo-gql-service). To regenerate the module using the current production schema:

```shell
HASURA_ADMIN_SECRET=... poetry run ariadne-codegen --config ariadne-prod.toml
```

# Scratch files

Some Python files in `./scratch` are available to seed existing data (i.e. legacy uBAM files).
