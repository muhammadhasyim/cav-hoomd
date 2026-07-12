# Packed Aging SLURM Submission Design

## Goal

Top up each aging-campaign coupling to a target of 500 scientifically usable
replicas while reducing the number of SLURM array tasks. No new work is
scheduled beyond the exact top-up requirement. Each array task will run two
replicas concurrently on one A100 GPU. Existing valid results count toward the
target and are never deleted merely because a coupling already has more than
500 replicas.

## Scope

The design applies to the aging campaign with coupling values
`0.0`, `0.01`, `0.016667`, `0.023333`, and `0.03`, a 2500 ps runtime, and
13 structural-correlation reference times. It covers:

- scientifically robust completion checks;
- selection of the replicas required to reach 500 valid results per coupling;
- pairing selected replicas into manifest-backed array tasks;
- concurrent execution of two replicas on one A100;
- resource and walltime selection;
- dry-run output, failure recovery, and automated tests; and
- submission and queue verification without affecting the `wrap` job.

It does not delete valid replicas above the new target, alter simulation physics,
change analysis definitions, or modify the unrelated `wrap` workflow.

## Current Evidence

Single-replica aging jobs normally take 56 to 58 minutes. Across completed jobs,
the measured 99th-percentile runtime is approximately 60 to 66 minutes, excluding
one anomalous control run. The existing 24-hour request is therefore much longer
than required.

The old HDF5 byte-size completion rule is insufficient. It classified three
short `lambda=0.016667` trajectories as complete and did not detect four
interleaved or malformed `lambda=0` correlation sets. Submission planning must
use the data products required by the relaxation analysis, not HDF5 size alone.

## Target Selection

The target is 500 total valid replicas per coupling, including every valid
existing replica in the configured ID domain `0..999`.

For each coupling:

1. Scan replica IDs `0..999`.
2. Count scientifically valid existing replicas.
3. Compute `needed = max(0, 500 - valid_count)`.
4. Select the lowest missing or invalid replica IDs until `needed` IDs have been
   selected.
5. Submit nothing when the valid count is already at least 500.

This means no new `lambda=0` replicas are expected. `lambda=0.01` is topped up
only if robust validation finds fewer than 500 valid outputs. The exact counts
for the remaining coupling values are computed immediately before submission.

## Scientific Completion Policy

A replica is valid only when all of the following hold:

1. `observables_replica_<id>.h5` exists and meets the established minimum-size
   guard.
2. All 13 files named `prod-<id>_fkt_ref_<ref>.txt` exist.
3. Every non-comment correlation row contains exactly two finite numeric values.
4. Lag times are strictly increasing and contain no duplicates.
5. The reference-time header is present.
6. Each reference reaches its expected final lag,
   `runtime_ps - reference_time_ps`, within 2.0 ps.

The robust policy is used for target selection and task-time skip checks.
The lightweight HDF5-only policy may remain available for quick inventory, but
it cannot authorize skipping a production replica in this campaign.

Invalid artifacts are removed only for replica IDs selected for rerun. Cleanup
includes the replica HDF5 and GSD files and its `prod-<id>_fkt_ref_*.txt` files.
Valid siblings and replicas above the target are preserved.

## Manifest and Array Layout

Selected replica IDs are grouped into pairs independently within each coupling.
A group may contain one replica when a coupling has an odd number remaining.

Groups from underfilled coupling values are interleaved into one combined
manifest so all coupling strengths begin accumulating data. Each manifest row
contains:

```text
lambda_tag<TAB>lambda_value<TAB>replica_1<TAB>replica_2_or_empty
```

The combined array has one index per manifest row and a global `%4` throttle.
Consequently, at most four A100 jobs and eight simulations run concurrently.
Using one combined array avoids accidentally allowing four active jobs for each
of several independently submitted arrays.

## Array-Task Execution

Each array task:

1. Reads its coupling and one or two replica IDs from the manifest.
2. Revalidates each replica in case it completed after manifest generation.
3. Cleans only invalid artifacts for replicas that still need execution.
4. Starts both simulation commands concurrently, each with a unique replica ID
   and the same allocated GPU 0.
5. Records each replica's standard output and standard error separately.
6. Waits for both child processes, preserving a successful result if its sibling
   fails.
7. Exits nonzero if either replica fails.

Signal traps forward termination to both child processes and wait for them,
preventing orphaned simulations and overlapping writers.

## Resources and Walltime

Each packed array task requests:

- one A100 GPU;
- eight CPU cores;
- 32 GB host memory;
- a three-hour walltime; and
- the existing chemistry account and A100 partition.

The three-hour limit is based on twice the measured single-replica p99 runtime,
plus operational margin. This covers the conservative case where two GPU
processes serialize most of their work while avoiding the scheduling cost of a
24-hour request.

## Failure Recovery

If one child succeeds and one fails, the valid child remains complete. A later
resume scan selects only the invalid replica. Manifest generation is
deterministic for a fixed filesystem snapshot, so dry-run and production plans
can be compared exactly.

The production submission command refuses to submit an empty plan. It also
prints per-coupling valid, needed, selected, and grouped counts before invoking
`sbatch`.

## Dry Run and Operational Flow

Dry-run mode performs the full robust scan and manifest planning without
deleting artifacts or calling `sbatch`. It reports:

- valid and invalid counts per coupling;
- selected replica IDs;
- pair count and singleton count;
- total array size and `%4` throttle;
- requested CPUs, memory, GPU, and walltime; and
- the exact intended `sbatch` command.

After tests and dry-run verification, production execution will:

1. confirm that only the unrelated `wrap` job is running;
2. generate the final manifest and packed SBATCH file;
3. submit the combined array;
4. inspect `squeue` to verify job name, task count, throttle, and states; and
5. leave `wrap` untouched.

## Testing Strategy

Tests are written before implementation and cover:

- selecting only enough missing IDs to reach a total target;
- submitting nothing when valid results already meet or exceed 500;
- pairing even and odd replica counts;
- preventing pairs from mixing coupling values;
- deterministic interleaving and manifest serialization;
- rejection of missing, malformed, duplicate, out-of-order, and short FKT data;
- preservation of valid replicas during cleanup;
- dry-run behavior with no `sbatch` side effect; and
- packed-task exit status when zero, one, or both child commands fail.

The existing aging status tests remain green. Shell syntax checks and a complete
dry run are required before any production submission.

## Acceptance Criteria

The work is accepted when:

1. all automated tests and shell syntax checks pass;
2. dry-run reports no new `lambda=0` work, selects no coupling already at 500,
   and never selects more replicas than its exact top-up requirement;
3. every manifest row contains at most two replicas from one coupling;
4. the generated array uses one A100, eight CPUs, 32 GB, three hours, and `%4`;
5. production submission leaves `wrap` running;
6. the queued array count matches the generated manifest; and
7. a failed child can be resumed without rerunning its successful sibling.
