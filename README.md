Depmap Omics Long Read RNA
---

This repository contains WDL workflows and a Python module to continuously run workflows on new RNA long read samples and sync BAM files between Terra and Gumbo.

# Installation

1. Install [uv](https://docs.astral.sh/uv).

2. Create a new virtualenv and install the project dependencies:

   ```sh
   uv venv
   uv sync
   ```

3. Copy `env.dist` into a new file `.env` and fill it out:

    - `WOMTOOL_JAR` is the path to the .jar file downloaded from [https://github.com/broadinstitute/cromwell/releases](broadinstitute/cromwell)
    - `FIRECLOUD_OWNERS` is a list of Terra users (email addresses or group names) that should be considered owners of workflow method configs

## Credentials

Your GCP `DEFAULT_APPLICATION_CREDENTIALS` must already be configured in order to run commands.

# Development

Configure your editor or IDE to automatically format your code with [Ruff](https://docs.astral.sh/ruff/) and check static types with [Pyright](https://microsoft.github.io/pyright) by running `uv run pyright`.

Whenever possible, function/method arguments and return values should be validated with Pydantic or Pandera (if a data frame).

## GraphQL code generation

This repo uses [ariadne-codegen](https://github.com/mirumee/ariadne-codegen) to generate the `gumbo_gql_client` module. It uses the folder of GraphQL queries (`./gql`) and the current GraphQL schema for a particular Gumbo environment to generate all of the Python classes, Pydantic models, and query/mutation methods for interacting with the [Gumbo GraphQL Service](https://github.com/broadinstitute/gumbo_client/tree/main/gumbo-gql-service). To regenerate the module using the current production schema:

```shell
HASURA_ADMIN_SECRET=... uv run ariadne-codegen --config ariadne-prod.toml
```

# Scratch files

Some Python files in `./scratch` are available to seed existing data (i.e. legacy uBAM files).
