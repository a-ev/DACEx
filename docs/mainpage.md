# DACEx API Documentation

DACEx extends DACE-based differential algebra workflows with Taylor and Chebyshev polynomial families.

## Contents

- `PolynomialBase`: shared setup, coefficient access, evaluation, and utility wrappers.
- `Polynomial<Taylor>`: Taylor-family algebra and transcendental operations.
- `Polynomial<Chebyshev>`: Chebyshev-family algebra and series-oriented operations.

## Build Docs Locally

```bash
doxygen docs/Doxyfile
```

Generated HTML output is written to `_docs/html/`.
