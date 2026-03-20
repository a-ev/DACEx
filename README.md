# DACEx

DACEx extends DACE with additional polynomial-family support and example/test programs.

## Requirements

- CMake >= 3.15
- A C++17 compiler
- DACE installed or built locally

## DACE Setup

DACEx depends on DACE. A common workflow is:

```bash
git clone https://github.com/dacelib/dace.git dace
cmake -S dace -B dace-build
cmake --build dace-build
cmake --install dace-build   # optional
```

## First-Time Build (Recommended)

### Windows (Visual Studio 2022, x64)

Always use `x64` unless you intentionally built a 32-bit DACE.

```powershell
if (Test-Path _build) { Remove-Item -Recurse -Force _build }
cmake -S . -B _build -G "Visual Studio 17 2022" -A x64 -T host=x64
cmake --build _build --config Debug -j
```

Run examples:

```powershell
.\_build\bin\Debug\DACEx_Gaussian_example.exe
.\_build\bin\Debug\DACEx_Gaussian_2D_example.exe
```

### Linux / macOS (Ninja or Make)

```bash
rm -rf _build
cmake -S . -B _build -DCMAKE_BUILD_TYPE=Debug
cmake --build _build -j
```

Run examples:

```bash
./_build/bin/DACEx_Gaussian_example
./_build/bin/DACEx_Gaussian_2D_example
```

## How DACE Discovery Works

DACEx tries, in order:

1. `find_package(dace CONFIG)`
2. Common local roots such as `~/dace-install` and `~/dace-build`
3. Manual overrides via CMake cache variables

Useful overrides:

- `-Ddace_DIR=<path containing daceConfig.cmake>`
- `-DCMAKE_PREFIX_PATH=<dace install prefix>`
- `-DDACEX_DACE_ROOT=<dace root with include/ and lib/>`
- `-DDACEX_DACE_INCLUDE_DIR=<dir containing dace/dace.h>`
- `-DDACEX_DACE_LIBRARY=<full path to dace library file>`

## Notes

- `DACEX_FETCH_DACE` is `OFF` by default. DACEx will not auto-install DACE unless explicitly enabled.
- On Windows, DACEx blocks Win32 by default (`DACEX_ALLOW_WIN32=OFF`) because most users have x64 DACE binaries.

## License

DACEx is open source and licensed under the Apache License 2.0. See `LICENSE`.

Copyright ownership is retained by the project author:

- Copyright (c) 2026 Adam Evans (@a-ev)

Attribution and ownership details are also captured in `NOTICE` and `AUTHORS.md`.


