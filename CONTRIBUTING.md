# Contributing

Open an issue before changing an inferential default or public return schema.
Bug fixes should include a minimal reproducible example and a regression test.
New pathway statistics should include size-stratified type-I error simulations,
power scenarios, and a plain-language interpretation before they are promoted
from experimental status.

Run the package checks before submitting a pull request:

```bash
R CMD build .
R CMD check --no-manual PASSAGE_*.tar.gz
```

Do not commit controlled-access data, identifiable sample metadata, credentials,
cluster logs, or generated analysis outputs. Public example data and gene sets
must have verified redistribution terms, source attribution, checksums and a
description of preprocessing. Do not add bulk MSigDB collections or gene sets
whose licenses do not permit redistribution.
