#!/bin/bash
set -euo pipefail

# Download the cellranger-8.0.1.tar.gz file from the specified URL
wget -O cellranger-8.0.1.tar.gz \
  "https://cf.10xgenomics.com/releases/cell-exp/cellranger-8.0.1.tar.gz?Expires=1721449718&Key-Pair-Id=APKAI7S6A5RYOXBWRPDA&Signature=iD8KiXRVmnegl~HthVrAk~tUZEDYTGa~bhldJVI08dZ6dXNeoCqy-2DjuFQDHTrw2mU~GR9LWfecyKyv9hIHg1giU5wmWEgv~UGycoTGQ5ntZb~ISQ3fjiLAMLLs2KVEDnS4iX13M9CAcUvJeULTq0s-KJ97GK9-mPBDfzH18Ky3uGJK~Pr3WAZ-YdJ1rN0rE315dl8M9EVacW89v5L5gYCf34vUX0Wjki9SLSJYD1jcNULTCVjZytllY8zgyviEYG4pv31X3sPQDKnl7Kofi3zn395eybwnw0qQWOgPZENTEFur4DdNfTdUqu6e0jF3Bnk0tZe486mCmpNFjcljQg__"
