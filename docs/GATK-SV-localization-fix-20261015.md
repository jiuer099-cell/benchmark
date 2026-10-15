# GATK-SV input-localization repair

Attempt 003 failed in DetermineGermlineContigPloidyCaseMode: the official
task invokes `gunzip` on localized counts. Cromwell's soft-link localization
created a link to a valid upstream gzip file. The pinned GATK image's gunzip
rejects that link with ELOOP ("Too many levels of symbolic links"); this was
not an actual recursive link on disk.

The backend now uses independent copies and disallows Docker soft links.
No hard-link fallback is allowed because tasks may modify localized inputs.
The official WDL, images, biological inputs and filters are unchanged.
The main adapter also receives the already-used gzip evidence reader and
read-only checkpoint mounting, which previously existed only in recovery code.
Stage plans bind the adapter, backend and container helper hashes.

Validation on the server: 26 unit tests passed. A real Cromwell 91/Apptainer
two-task regression used the original 9,822,124-byte condensed counts file
and the pinned GATK image. The old soft-link operation reproduced ELOOP;
the copy backend completed both tasks, preserved the input and upstream
producer output, and returned the exact decompressed bytes.
Source SHA256: `64bd75a7adc9f92b60e4ac7dd8ca5292fa5bdc23567e6e6b2ab8d2a396aa819e`.

Recovery uses a new `HG002-attempt-004-copy-localization` directory and a
frozen helper snapshot. Authenticated attempt-002 RG/Manta/GatherSampleEvidence
outputs are reused through official case inputs; downstream stages rerun.
Attempt 003 remains a failed historical run. Successful localization testing
does not imply the complete HG002 workflow has finished or produced a score.
