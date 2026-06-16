# ⟨7,7,7⟩ recombination multisets — base ⟨2,2,2⟩

Every recombination multiset the ⟨2,2,2⟩ base can realise (exact GL-orbit via `RecombinationMultisetOrbit`), and the number of multiplications each yields when recombined to ⟨7,7,7⟩.

**The allocation is independent per axis (n, m, p).** Each axis is split into 2/2/2 parts; the `best` column is the minimum over the *joint* per-axis split space and `best alloc` is the winning `(n / m / p)` split. The `a+b` columns are the *symmetric* splits (same on all three axes) for reference — the per-axis optimum can beat them.

- **40 multisets** total over ℚ, **6 ternary** (constructable, cited by base hash; the other 34 need a rational change-of-basis — predictions only).
- `multiset` is allocation-independent (B = larger block ⌈⌉, S = smaller block ⌊⌋; the per-product triple is axis-sorted); **bold** = column-best.

| #  | multiset                      | 4+3     | 5+2     | 6+1     | best    | best alloc (n / m / p) | ternary | base (hash)                                                    |
| --- | ----------------------------- | ------- | ------- | ------- | ------- | ---------------------- | ------- | -------------------------------------------------------------- |
| 1  | 1×BBB + 3×BBS + 3×BSS         | **249** | **267** | **279** | **249** | 4+3 / 4+3 / 4+3        | ✓       | 2x2x2-r7-strassen-db11bcc.json (@db11bcc)                      |
| 2  | 2×BBB + 1×BBS + 4×BSS         | 250     | 298     | 366     | 250     | 4+3 / 4+3 / 4+3        | ✓       | 2x2x2-r7-winograd_1971-511df05.json (@511df05)                 |
| 3  | 1×BBB + 4×BBS + 1×BSS + 1×SSS | 252     | 278     | 304     | 252     | 4+3 / 4+3 / 4+3        | ✓       | 2x2x2-r7-perminov_cr15_cn24_ZT_reduced-e498eb7.json (@e498eb7) |
| 4  | 2×BBB + 2×BBS + 2×BSS + 1×SSS | 253     | 309     | 391     | 253     | 4+3 / 4+3 / 4+3        | ✓       | 2x2x2-r7-alphatensor_Z-18678ba.json (@18678ba)                 |
| 5  | 3×BBB + 3×BSS + 1×SSS         | 254     | 340     | 478     | 254     | 4+3 / 4+3 / 4+3        | ✓       | 2x2x2-r7-solven_orbit_c9x3_c8x1-37b6062.json (@37b6062)        |
| 6  | 2×BBB + 3×BBS + 2×SSS         | 256     | 320     | 416     | 256     | 4+3 / 4+3 / 4+3        | ✓       | 2x2x2-r7-solven_orbit_c9x2_c8x2-c9abd45.json (@c9abd45)        |
| 7  | 1×BBB + 4×BBS + 2×BSS         | 258     | 289     | 309     | 258     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 8  | 2×BBB + 2×BBS + 3×BSS         | 259     | 320     | 396     | 259     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 9  | 2×BBB + 2×BBS + 3×BSS         | 259     | 320     | 396     | 259     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 10 | 2×BBB + 3×BBS + 1×BSS + 1×SSS | 262     | 331     | 421     | 262     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 11 | 2×BBB + 3×BBS + 1×BSS + 1×SSS | 262     | 331     | 421     | 262     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 12 | 3×BBB + 1×BBS + 2×BSS + 1×SSS | 263     | 362     | 508     | 263     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 13 | 1×BBB + 5×BBS + 1×BSS         | 267     | 311     | 339     | 267     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 14 | 1×BBB + 5×BBS + 1×BSS         | 267     | 311     | 339     | 267     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 15 | 2×BBB + 3×BBS + 2×BSS         | 268     | 342     | 426     | 268     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 16 | 2×BBB + 3×BBS + 2×BSS         | 268     | 342     | 426     | 268     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 17 | 2×BBB + 3×BBS + 2×BSS         | 268     | 342     | 426     | 268     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 18 | 3×BBB + 1×BBS + 3×BSS         | 269     | 373     | 513     | 269     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 19 | 3×BBB + 1×BBS + 3×BSS         | 269     | 373     | 513     | 269     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 20 | 2×BBB + 4×BBS + 1×SSS         | 271     | 353     | 451     | 271     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 21 | 3×BBB + 2×BBS + 1×BSS + 1×SSS | 272     | 384     | 538     | 272     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 22 | 3×BBB + 2×BBS + 1×BSS + 1×SSS | 272     | 384     | 538     | 272     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 23 | 4×BBB + 2×BSS + 1×SSS         | 273     | 415     | 625     | 273     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 24 | 1×BBB + 6×BBS                 | 276     | 333     | 369     | 276     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 25 | 2×BBB + 4×BBS + 1×BSS         | 277     | 364     | 456     | 277     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 26 | 2×BBB + 4×BBS + 1×BSS         | 277     | 364     | 456     | 277     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 27 | 3×BBB + 2×BBS + 2×BSS         | 278     | 395     | 543     | 278     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 28 | 3×BBB + 2×BBS + 2×BSS         | 278     | 395     | 543     | 278     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 29 | 3×BBB + 2×BBS + 2×BSS         | 278     | 395     | 543     | 278     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 30 | 3×BBB + 3×BBS + 1×SSS         | 281     | 406     | 568     | 281     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 31 | 4×BBB + 1×BBS + 1×BSS + 1×SSS | 282     | 437     | 655     | 282     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 32 | 2×BBB + 5×BBS                 | 286     | 386     | 486     | 286     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 33 | 3×BBB + 3×BBS + 1×BSS         | 287     | 417     | 573     | 287     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 34 | 4×BBB + 1×BBS + 2×BSS         | 288     | 448     | 660     | 288     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 35 | 3×BBB + 4×BBS                 | 296     | 439     | 603     | 296     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 36 | 4×BBB + 2×BBS + 1×BSS         | 297     | 470     | 690     | 297     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 37 | 5×BBB + 2×BSS                 | 298     | 501     | 777     | 298     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 38 | 4×BBB + 3×BBS                 | 306     | 492     | 720     | 306     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 39 | 5×BBB + 2×BBS                 | 316     | 545     | 837     | 316     | 4+3 / 4+3 / 4+3        | —       | —                                                              |
| 40 | 7×BBB                         | 336     | 651     | 1071    | 336     | 4+3 / 4+3 / 4+3        | —       | —                                                              |

## Block structure (per ternary base)

The cost-multiset above uses only block *size* (B/S). But a `B` reaches the big size two ways: a **plain** big block (touches only the big block) or a **mixed** `B⊕S` sum (touches both). Same cost — yet the mixed block carries small-block support that blocks the further shaving a plain `S` would permit, which is *why* two bases can share a B/S cost-multiset and still be distinct. Legend: **B** plain big, **M** mixed big (B⊕S), **S** small.

| base (hash)                                                    | fine multiset (B / M / S)                     |
| -------------------------------------------------------------- | --------------------------------------------- |
| 2x2x2-r7-perminov_cr15_cn24_ZT_reduced-e498eb7.json (@e498eb7) | 2×BMS + 1×BSS + 1×MMM + 2×MMS + 1×SSS         |
| 2x2x2-r7-alphatensor_Z-18678ba.json (@18678ba)                 | 1×BBM + 1×BMS + 1×MMM + 1×MMS + 2×MSS + 1×SSS |
| 2x2x2-r7-winograd_1971-511df05.json (@511df05)                 | 1×BBS + 1×BMM + 1×MMM + 4×MSS                 |
| 2x2x2-r7-strassen-db11bcc.json (@db11bcc)                      | 3×BSS + 1×MMM + 3×MMS                         |
| 2x2x2-r7-solven_orbit_c9x3_c8x1-37b6062.json (@37b6062)        | 3×BSS + 3×MMM + 1×SSS                         |
| 2x2x2-r7-solven_orbit_c9x2_c8x2-c9abd45.json (@c9abd45)        | 1×BBB + 1×MMM + 3×MMS + 2×SSS                 |

_Generated by `RecombinationMultisetReport` (target ⟨7,7,7⟩, base ⟨2,2,2⟩)._
