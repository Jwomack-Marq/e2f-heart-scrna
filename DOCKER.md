# Running the app with Docker

This packages the **live server-side R Shiny app** (`shiny_app/app.R` + `app_data.rds`)
into a container that runs on a real R/Shiny server. It is independent of the shinylive
static export used for GitHub Pages.

## Prerequisites

- [Docker](https://docs.docker.com/get-docker/) installed (Docker Desktop on macOS/Windows,
  or Docker Engine on Linux). Docker Compose ships with recent Docker.

## Quick start (Docker Compose)

```bash
docker compose up --build
```

Then open **http://localhost:3838**. Stop with `Ctrl+C`, or `docker compose down`.

## Quick start (plain Docker)

```bash
docker build -t e2f-heart-scrna .
docker run --rm -p 3838:3838 e2f-heart-scrna
```

Open **http://localhost:3838**.

## Changing the port

Map any host port to the container's `3838`. For example, to serve on `8080`:

```bash
docker run --rm -p 8080:3838 e2f-heart-scrna     # http://localhost:8080
```

With Compose, edit the `ports` line in `docker-compose.yml` to `"8080:3838"`.

## Notes

- **First build takes a few minutes** (installs R packages); subsequent builds are cached.
  Startup then takes a few seconds while R loads packages and reads the ~23 MB data file.
- **Data:** the image bundles `shiny_app/app_data.rds` (the current slim data build).
- **Concurrency:** one R process serves all sessions — fine for a demo or a modest launch.
  For heavier traffic, run multiple replicas behind a load balancer (or move to
  shiny-server / ShinyProxy).
- **HTTPS / public launch:** this setup serves plain HTTP for local use. Before exposing it
  publicly (e.g. linked from the splash page), put an HTTPS reverse proxy in front
  (Caddy / Nginx / Traefik) or deploy to a platform that terminates TLS (Cloud Run,
  Fly.io, Render) and forwards to container port `3838`.
