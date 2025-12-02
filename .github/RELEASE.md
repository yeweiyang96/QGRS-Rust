# GitHub Actions Release Guide

This repository uses GitHub Actions to automatically build and release binaries for multiple platforms.

## How It Works

### Automated CI (`ci.yml`)
- **Triggers**: Push/PR to `master`, `main`, or `develop` branches
- **Actions**:
  - Runs tests on Linux, macOS, and Windows
  - Runs clippy linter
  - Checks code formatting
- **Purpose**: Ensures code quality before merging

### Automated Release (`release.yml`)
- **Triggers**: When you push a version tag (e.g., `v1.0.0`, `v0.2.1`)
- **Actions**:
  1. Builds binaries for 5 platforms:
     - `x86_64-unknown-linux-gnu` (Linux x64)
     - `aarch64-unknown-linux-gnu` (Linux ARM64)
     - `x86_64-apple-darwin` (macOS Intel)
     - `aarch64-apple-darwin` (macOS Apple Silicon)
     - `x86_64-pc-windows-msvc` (Windows x64)
  2. Packages binaries with README and LICENSE
  3. Generates SHA256 checksums
  4. Creates a GitHub Release with all archives

## How to Trigger a Release

### 1. Update version in Cargo.toml (optional but recommended)
```toml
[package]
name = "QGRS-Rust"
version = "1.0.0"  # Update this
```

### 2. Commit your changes
```bash
git add .
git commit -m "Release v1.0.0"
git push origin master
```

### 3. Create and push a version tag
```bash
# Create a tag
git tag v1.0.0

# Push the tag to GitHub
git push origin v1.0.0
```

### 4. Watch the action run
- Go to your repository on GitHub
- Click on "Actions" tab
- You'll see the "Release" workflow running
- Wait for all builds to complete (~10-15 minutes)

### 5. Check your release
- Go to "Releases" tab
- Your new release will be published with downloadable binaries

## Release Archives

Each platform gets its own archive:

**Unix/macOS** (`.tar.gz`):
```
qgrs-linux-x64.tar.gz
qgrs-linux-arm64.tar.gz
qgrs-macos-x64.tar.gz
qgrs-macos-arm64.tar.gz
```

**Windows** (`.zip`):
```
qgrs-windows-x64.zip
```

Each archive contains:
- `qgrs` (or `qgrs.exe` on Windows)
- `compare_modes` (or `.exe`)
- `compare_csv_outputs` (or `.exe`)
- `README.md`
- `LICENSE`

Each archive also has a `.sha256` checksum file for verification.

## Version Numbering

Follow [Semantic Versioning](https://semver.org/):
- `vMAJOR.MINOR.PATCH` (e.g., `v1.2.3`)
- Increment MAJOR for breaking changes
- Increment MINOR for new features
- Increment PATCH for bug fixes

Examples:
```bash
git tag v0.1.0   # Initial release
git tag v0.1.1   # Bug fix
git tag v0.2.0   # New feature
git tag v1.0.0   # Stable release
```

## Troubleshooting

### Build fails for a specific platform
- Check the Actions log for that platform
- Common issues:
  - Missing dependencies (especially on ARM64 Linux)
  - Compilation errors (fix code and re-tag)
  - Network issues (re-run the workflow)

### How to re-run a failed build
1. Go to Actions tab
2. Click on the failed workflow
3. Click "Re-run all jobs" button

### How to delete a release and re-create it
```bash
# Delete the tag locally
git tag -d v1.0.0

# Delete the tag remotely
git push origin :refs/tags/v1.0.0

# Delete the release on GitHub (via web UI)
# Then recreate the tag and push again
git tag v1.0.0
git push origin v1.0.0
```

## Manual Release (if needed)

If you need to build locally instead:

```bash
# Build for your current platform
cargo build --release --all-targets

# Binaries will be in target/release/
ls target/release/qgrs
ls target/release/compare_modes
ls target/release/compare_csv_outputs
```

For cross-compilation, use [cross](https://github.com/cross-rs/cross):
```bash
cargo install cross
cross build --release --target x86_64-unknown-linux-gnu
```

## Notes

- The workflow requires GitHub Actions to be enabled on your repository
- The `GITHUB_TOKEN` is automatically provided by GitHub Actions
- Release drafts are disabled by default (change `draft: false` to `draft: true` in release.yml if needed)
- Pre-releases can be marked by changing `prerelease: false` to `true` or using tags like `v1.0.0-beta.1`
