default:
    @just --list

[private]
warn := "\\033[33m"
error := "\\033[31m"
info := "\\033[34m"
success := "\\033[32m"
reset := "\\033[0m"
bold := "\\033[1m"

# Print formatted headers without shell scripts
[private]
header msg:
    @printf "{{info}}{{bold}}==> {{msg}}{{reset}}\n"

# Get native architecture
[private]
native_arch := if `uname -m` == "arm64" { "aarch64" } else { `uname -m` }

# Install cargo tools
install-cargo-tools:
    @just header "Installing Cargo tools"
    # cargo-udeps
    if ! command -v cargo-udeps > /dev/null; then \
        printf "{{info}}Installing cargo-udeps...{{reset}}\n" && \
        cargo install cargo-udeps --locked; \
    else \
        printf "{{success}}✓ cargo-udeps already installed{{reset}}\n"; \
    fi
    # cargo-semver-checks
    if ! command -v cargo-semver-checks > /dev/null; then \
        printf "{{info}}Installing cargo-semver-checks...{{reset}}\n" && \
        cargo install cargo-semver-checks; \
    else \
        printf "{{success}}✓ cargo-semver-checks already installed{{reset}}\n"; \
    fi
    # taplo
    if ! command -v taplo > /dev/null; then \
        printf "{{info}}Installing taplo...{{reset}}\n" && \
        cargo install taplo-cli; \
    else \
        printf "{{success}}✓ taplo already installed{{reset}}\n"; \
    fi
    # wasm-pack
    if ! command -v wasm-pack > /dev/null; then \
        printf "{{info}}Installing wasm-pack...{{reset}}\n" && \
        cargo install wasm-pack; \
    else \
        printf "{{success}}✓ wasm-pack already installed{{reset}}\n"; \
    fi

# Install nightly rust
install-rust-nightly:
    @just header "Installing Rust nightly"
    rustup install nightly
    rustup component add clippy --toolchain nightly
    rustup component add rustfmt --toolchain nightly

# Install required Rust targets
install-targets:
    @just header "Installing Rust targets"
    rustup target add \
        aarch64-unknown-linux-gnu \
        x86_64-unknown-linux-gnu \
        wasm32-unknown-unknown

# Setup complete development environment
setup: install-cargo-tools install-rust-nightly install-targets
    @printf "{{success}}{{bold}}Development environment setup complete!{{reset}}\n"

# Build with local OS target
build:
    @just header "Building workspace"
    cargo build --workspace --all-targets

# Build all architecture targets
build-all: build-aarch64 build-x86 build-wasm
    @printf "{{success}}{{bold}}All arch builds completed!{{reset}}\n"

# Build aarch64 target
build-aarch64:
    @just header "Building aarch64"
    cargo build --workspace --target aarch64-unknown-linux-gnu

# Build x86_64 target
build-x86:
    @just header "Building x86_64"
    cargo build --workspace --target x86_64-unknown-linux-gnu

# Build wasm target
build-wasm:
    @just header "Building wasm32"
    cargo build --workspace --target wasm32-unknown-unknown

# Run tests for native architecture and wasm
test:
    @just header "Running native architecture tests"
    cargo test --workspace --all-targets --all-features
    @just header "Running wasm tests"
    wasm-pack test --node

# Run clippy for the workspace
lint:
    @just header "Running clippy"
    cargo clippy --workspace --all-targets --all-features -- --deny warnings

# Run clippy for all targets
lint-all: lint-aarch64 lint-x86 lint-wasm
    @printf "{{success}}{{bold}}All arch lint completed!{{reset}}\n"

# Run clippy for the aarch64 target
lint-aarch64:
    @just header "Checking lint on aarch64"
    cargo clippy --workspace --all-targets --target aarch64-unknown-linux-gnu -- --deny warnings

# Run clippy for the x86_64 target
lint-x86:
    @just header "Checking lint on x86_64"
    cargo clippy --workspace --all-targets --target x86_64-unknown-linux-gnu -- --deny warnings

# Run clippy for the wasm target
lint-wasm:
    @just header "Checking lint on wasm32"
    cargo clippy --workspace --all-targets --target wasm32-unknown-unknown -- --deny warnings

# Check for semantic versioning
semver:
    @just header "Checking semver compatibility"
    rustup override set stable && cargo semver-checks check-release --workspace && rustup override unset

# Run format for the workspace
fmt:
    @just header "Formatting code"
    cargo fmt --all
    taplo fmt

# Check for unused dependencies
udeps:
    @just header "Checking unused dependencies"
    cargo +nightly udeps --workspace

# Clean build artifacts
clean:
    @just header "Cleaning build artifacts"
    cargo clean

# Show environment information
info:
    @just header "Environment Information"
    @printf "{{info}}OS:{{reset}} %s\n" "$(uname -s)"
    @printf "{{info}}Architecture:{{reset}} %s\n" "$(uname -m)"
    @printf "{{info}}Rust:{{reset}} %s\n" "$(rustc --version)"
    @printf "{{info}}Cargo:{{reset}} %s\n" "$(cargo --version)"
    @printf "{{info}}Installed targets:{{reset}}\n"
    @rustup target list --installed | sed 's/^/  /'

# Run all CI checks
ci:
    @printf "{{bold}}Starting CI checks{{reset}}\n\n"
    @ERROR=0; \
    just run-single-check "Rust formatting" "cargo fmt --all -- --check" || ERROR=1; \
    just run-single-check "TOML formatting" "taplo fmt --check" || ERROR=1; \
    just run-single-check "aarch64 build" "cargo build --target aarch64-unknown-linux-gnu --workspace" || ERROR=1; \
    just run-single-check "x86_64 build" "cargo build --target x86_64-unknown-linux-gnu --workspace" || ERROR=1; \
    just run-single-check "wasm32 build" "cargo build --target wasm32-unknown-unknown --workspace" || ERROR=1; \
    just run-single-check "aarch64 clippy" "cargo clippy --target aarch64-unknown-linux-gnu --all-targets --all-features -- --deny warnings" || ERROR=1; \
    just run-single-check "x86_64 clippy" "cargo clippy --target x86_64-unknown-linux-gnu --all-targets --all-features -- --deny warnings" || ERROR=1; \
    just run-single-check "wasm32 clippy" "cargo clippy --target wasm32-unknown-unknown --all-targets --all-features -- --deny warnings" || ERROR=1; \
    just run-single-check "Native tests" "cargo test --verbose --workspace" || ERROR=1; \
    just run-single-check "WASM tests" "wasm-pack test --node" || ERROR=1; \
    just run-single-check "Unused dependencies" "cargo +nightly udeps --workspace" || ERROR=1; \
    just run-single-check "Semver compatibility" "(rustup override set stable && cargo semver-checks check-release --workspace && rustup override unset)" || ERROR=1; \
    printf "\n{{bold}}CI Summary:{{reset}}\n"; \
    if [ $ERROR -eq 0 ]; then \
        printf "{{success}}{{bold}}All checks passed successfully!{{reset}}\n"; \
    else \
        printf "{{error}}{{bold}}Some checks failed. See output above for details.{{reset}}\n"; \
        exit 1; \
    fi

# Run a single check and return status (0 = pass, 1 = fail)
[private]
run-single-check name command:
    #!/usr/bin/env sh
    printf "{{info}}{{bold}}Running{{reset}} {{info}}%s{{reset}}...\n" "{{name}}"
    if {{command}} > /tmp/check-output 2>&1; then
        printf "  {{success}}{{bold}}PASSED{{reset}}\n"
        exit 0
    else
        printf "  {{error}}{{bold}}FAILED{{reset}}\n"
        printf "{{error}}----------------------------------------\n"
        while IFS= read -r line; do
            printf "{{error}}%s{{reset}}\n" "$line"
        done < /tmp/check-output
        printf "{{error}}----------------------------------------{{reset}}\n"
        exit 1
    fi