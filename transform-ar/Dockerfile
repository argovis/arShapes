FROM rust:1.70.0

RUN apt-get update -y ; apt-get install -y nano netcdf-bin libhdf5-serial-dev libnetcdff-dev

WORKDIR /app
COPY Cargo.toml Cargo.toml
COPY src/main.rs src/main.rs
RUN chown -R 1000660000 /app
