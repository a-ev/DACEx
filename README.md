# DACEx

DACEx extends the [DACE](https://github.com/dacelib/dace) library to additional polynomial families.

## Requirements

- CMake >= 3.15
- A C++17 compiler
- DACE installed or built locally

## DACE Setup

DACEx depends on [DACE](https://github.com/dacelib/dace). You can build DACE as per the instructions on the Wiki:

```bash
git clone https://github.com/dacelib/dace.git dace
cmake -S dace -B dace-build
cmake --build dace-build
cmake --install dace-build
```

## First-Time Build (Recommended)

Build is currently only configured for windows.

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

## DACE Discovery Notes

DACEx tries, in order:

1. `find_package(dace CONFIG)`
2. Common local roots such as `~/dace-install` and `~/dace-build`
3. Manual overrides via CMake cache variables


## Notes

- `DACEX_FETCH_DACE` is `OFF` by default. DACEx will not auto-install DACE unless explicitly enabled.
- On Windows, DACEx blocks Win32 by default (`DACEX_ALLOW_WIN32=OFF`) because most users have x64 DACE binaries.

## License

DACEx is open source and licensed under the Apache License 2.0. See `LICENSE`. Attribution and ownership details are also captured in `NOTICE` and `AUTHORS.md`.


