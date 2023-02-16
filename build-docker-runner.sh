#!/bin/bash

set -e

VERSION=$(git describe --tags `git rev-list --tags --max-count=1` --always)

# dynamically pull more interesting stuff from latest git commit
HASH=$(git show-ref --head --hash=8 head)  # first 8 letters of hash should be enough; that's what GitHub uses

# Build runner docker image
docker build -f Dockerfile.runner -t quartet-dseqc-report:runner-${VERSION}-${HASH} . && \
docker tag quartet-dseqc-report:runner-${VERSION}-${HASH} ghcr.io/chinese-quartet/quartet-dseqc-report:runner-${VERSION}-${HASH} && \
docker push ghcr.io/chinese-quartet/quartet-dseqc-report:runner-${VERSION}-${HASH}