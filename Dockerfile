# Fluxsheet: simulator, API and web editor in one container.
# Build:  docker build -t fluxsheet .
# Run:    docker run -p 8000:8000 fluxsheet      (open http://localhost:8000)
# The image runs with the demo limits on (FLUXSHEET_DEMO=1); set it to 0 for private use.

FROM node:20.18.1-alpine3.20 AS web
WORKDIR /src/web
COPY web/package.json web/package-lock.json ./
RUN npm ci --no-audit --no-fund
COPY web/ ./
RUN npm run build

FROM python:3.12.7-slim-bookworm
ENV PYTHONDONTWRITEBYTECODE=1 \
    PYTHONUNBUFFERED=1 \
    PIP_NO_CACHE_DIR=1 \
    PIP_DISABLE_PIP_VERSION_CHECK=1 \
    FLUXSHEET_DEMO=1
WORKDIR /src
COPY pyproject.toml hatch_build.py README.md LICENSE CHANGELOG.md ./
COPY src ./src
COPY examples ./examples
COPY --from=web /src/web/dist ./web/dist
RUN FLUXSHEET_REQUIRE_WEB=1 pip install ".[server]" \
 && rm -rf /src \
 && useradd --system --uid 10001 --no-create-home fluxsheet
WORKDIR /
USER fluxsheet
EXPOSE 8000
HEALTHCHECK --interval=60s --timeout=5s --start-period=60s \
  CMD python -c "import urllib.request; urllib.request.urlopen('http://127.0.0.1:8000/api/health', timeout=4)"
CMD ["fluxsheet", "serve", "--host", "0.0.0.0", "--port", "8000"]
