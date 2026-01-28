FROM python:3.11.14

ENV PYTHONUNBUFFERED=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=on \
    PIP_DEFAULT_TIMEOUT=100 \
    PYTHONPATH=/app \
    UV_VENV=/app/.venv

ENV PATH="$UV_VENV/bin:$PATH"

WORKDIR /app

RUN pip install uv

COPY pyproject.toml uv.lock README.md ./
COPY combine_requantify_tools/__init__.py ./combine_requantify_tools/__init__.py
RUN uv sync --no-cache --no-group=dev

COPY combine_requantify_tools ./combine_requantify_tools
