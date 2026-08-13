# Engineering record

Development history for the Rust implementation: what was measured, what was
tried, and why each decision went the way it did. None of this is user
documentation — for that, see [`../index.md`](../index.md) — and none of it is
shipped: `cargo package` includes only each crate's own `README.md`, so nothing
here reaches a published artifact.

It is kept because most of its value is *negative*. A large share of the
optimizations attempted here were measured and reverted, and without the record
the obvious ones invite being tried again. Anyone reaching for an idea below
should read its entry first.

## Contents

| file | what it holds |
| --- | --- |
| [performance-record.md](performance-record.md) | Every measurement, attempt and reversal, in chronological order |
| [review-findings.md](review-findings.md) | The audit that opened the optimization campaign, and the K > 31 defects |
| [cpp-nondeterminism.md](cpp-nondeterminism.md) | The C++ implementation dropping unitigs, and why it is not a reliable oracle |
| [color-query-index.md](color-query-index.md) | Design sketch for a k-mer → unitig locator (not built) |
| [module-boundaries.md](module-boundaries.md) | Safe decomposition boundaries for the large modules |
| [compat-harness.md](compat-harness.md) | How compatibility fixtures and expectations are produced |
| [edge-matrix-containers.md](edge-matrix-containers.md) | The blocked edge matrix's container layout |
| [roadmap.md](roadmap.md) | Original rewrite plan |

## Measured and reverted — do not re-attempt without new evidence

Each of these was implemented, measured, and removed. The performance record
carries the numbers and the reasoning.

| attempt | outcome |
| --- | --- |
| Inline `(key, state)` entries in the vertex map | +16% local contraction; the probed array grew and evicted |
| Flatten the wanted-colour map the same way | +17.9% on the colored walk, same reason |
| Membership filter in front of the wanted-colour map | Neutral; the table is L1-resident, so a miss was already cheap |
| Software-pipelined prefetch in `add_packed_parts` | Neutral at distance 8 and 16; the phase is bandwidth-bound, not latency-bound |
| Overlapped flush crew for the colored window | Worse at every crew size; the flush is CPU-bound, so the barrier was never wasting the machine |
| Read-side sort for colour exactness | +29%; destroyed the streaming locality of bucket consumption |
| Interleaved colored compression | 2.1% faster but 11.8% more intermediate bytes, paid again downstream |
| LTO, mimalloc, jemalloc decay tuning, bucket-compression defaults | See the rejected-experiments section of the performance record |

## What worked, and what the wins had in common

| change | effect |
| --- | --- |
| Flat prefetched vertex map (R1) | Local contraction −39% uncolored, −27% colored |
| Colour slot: uncolored builds stop carrying a colour hash | Uncolored build −3.0%; matches C++'s `State_Config` specialization |
| Parallel wide partition table | k = 55 contraction 15.0 s → 1.9 s |
| Streamed wide expansion | k = 55 row load/decode 2.46 s → 0.15 s |
| Skipping the colored window sort | Colored partition −2.1%; it was repairing payloads already correct |
| Split colored staging | Colored partition −2.7% at identical on-disk volume |

The pattern is worth stating: the largest wins came from finding work *this*
implementation was doing that the C++ one was not — a redundant sort, a
redundant deinterleave, a colour hash uncolored builds never read. The attempts
that came from profiling our own hot loops in isolation were mostly neutral.
