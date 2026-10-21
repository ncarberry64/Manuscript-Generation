"""Keep archived repair evidence separate from executable repository tests."""

# These directories preserve delivered/superseded test sources byte-for-byte.
# Their current executable counterparts live in tests/; collecting both copies
# creates import-name collisions and would execute historical evidence as tests.
collect_ignore = [
    "data/strong_lensing_response_repair_20261021/delivery",
    "artifacts/strong_lensing_response_repair_20261021",
]
