# Use the SRA v4.3.0 reads and the standard Clair3 v4.3.0 calling models

We reuse the v4.3.0 simplex hac and sup reads deposited in SRA rather than rebasecalling the figshare pod5s. Dorado's `--bacteria` model supports reads basecalled with v4.2.0 or later, and keeping the same reads makes the results directly comparable with the eLife paper. Ryan Wick's 2026 polishing post also found the bacterial model matched or beat the newer basecaller-specific polishing models, so rebasecalling was unlikely to change the conclusion.

Clair3 now ships a bacteria fine-tuned model (`r1041_e82_400bps_sup_v430_bacteria_finetuned`) that is the obvious "latest" choice. We deliberately don't use it: it was trained on 12 of our 14 samples (only ATCC_25922 and ATCC_10708 were held out), so scoring it here would leak training data into the test. The post must say this explicitly.

## Consequences

- `dorado smallvar` (hac v5.2/v6.0 only) and the newer polishing models can't be tested without rebasecalling. They are listed as caveats in the post, not benchmarked.
- Duplex reads are out of scope because `dorado polish` rejects duplex models.
