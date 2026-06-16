# ⟨17,17,17⟩ recombination multisets — base ⟨2,2,2⟩

Every recombination multiset the ⟨2,2,2⟩ base can realise (exact GL-orbit via `RecombinationMultisetOrbit`), and the number of multiplications each yields when recombined to ⟨17,17,17⟩.

**The allocation is independent per axis (n, m, p).** Each axis is split into 2/2/2 parts; the `best` column is the minimum over the *joint* per-axis split space and `best alloc` is the winning `(n / m / p)` split. The `a+b` columns are the *symmetric* splits (same on all three axes) for reference — the per-axis optimum can beat them.

- **40 multisets** total over ℚ, **6 ternary** (constructable, cited by base hash; the other 34 need a rational change-of-basis — predictions only).
- `multiset` is allocation-independent (B = larger block ⌈⌉, S = smaller block ⌊⌋; the per-product triple is axis-sorted); **bold** = column-best.

| #  | multiset                      | 9+8      | 10+7     | 11+6     | 12+5     | 13+4     | 14+3     | 15+2     | 16+1     | best     | best alloc (n / m / p) | ternary | base (hash)                                                    |
| --- | ----------------------------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | -------- | ---------------------- | ------- | -------------------------------------------------------------- |
| 1  | 1×BBB + 4×BBS + 1×BSS + 1×SSS | **2930** | 3157     | 3254     | 3289     | 3485     | 3589     | 3498     | 3345     | **2930** | 9+8 / 9+8 / 9+8        | ✓       | 2x2x2-r7-perminov_cr15_cn24_ZT_reduced-e498eb7.json (@e498eb7) |
| 2  | 2×BBB + 3×BBS + 2×SSS         | 2934     | 3234     | 3522     | 3730     | 4336     | 4798     | 5165     | 5378     | 2934     | 9+8 / 9+8 / 9+8        | ✓       | 2x2x2-r7-solven_orbit_c9x2_c8x2-c9abd45.json (@c9abd45)        |
| 3  | 1×BBB + 3×BBS + 3×BSS         | 2940     | **3120** | **3147** | **3116** | **3275** | **3318** | **3252** | **3120** | 2940     | 9+8 / 9+8 / 9+8        | ✓       | 2x2x2-r7-strassen-db11bcc.json (@db11bcc)                      |
| 4  | 2×BBB + 2×BBS + 2×BSS + 1×SSS | 2944     | 3197     | 3415     | 3557     | 4126     | 4527     | 4919     | 5153     | 2944     | 9+8 / 9+8 / 9+8        | ✓       | 2x2x2-r7-alphatensor_Z-18678ba.json (@18678ba)                 |
| 5  | 2×BBB + 1×BBS + 4×BSS         | 2954     | 3160     | 3308     | 3384     | 3916     | 4256     | 4673     | 4928     | 2954     | 9+8 / 9+8 / 9+8        | ✓       | 2x2x2-r7-winograd_1971-511df05.json (@511df05)                 |
| 6  | 3×BBB + 3×BSS + 1×SSS         | 2958     | 3237     | 3576     | 3825     | 4767     | 5465     | 6340     | 6961     | 2958     | 9+8 / 9+8 / 9+8        | ✓       | 2x2x2-r7-solven_orbit_c9x3_c8x1-37b6062.json (@37b6062)        |
| 7  | 1×BBB + 4×BBS + 2×BSS         | 2982     | 3253     | 3369     | 3400     | 3589     | 3661     | 3544     | 3360     | 2982     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 8  | 2×BBB + 3×BBS + 1×BSS + 1×SSS | 2986     | 3330     | 3637     | 3841     | 4440     | 4870     | 5211     | 5393     | 2986     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 9  | 2×BBB + 3×BBS + 1×BSS + 1×SSS | 2986     | 3330     | 3637     | 3841     | 4440     | 4870     | 5211     | 5393     | 2986     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 10 | 2×BBB + 2×BBS + 3×BSS         | 2996     | 3293     | 3530     | 3668     | 4230     | 4599     | 4965     | 5168     | 2996     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 11 | 2×BBB + 2×BBS + 3×BSS         | 2996     | 3293     | 3530     | 3668     | 4230     | 4599     | 4965     | 5168     | 2996     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 12 | 3×BBB + 1×BBS + 2×BSS + 1×SSS | 3000     | 3370     | 3798     | 4109     | 5081     | 5808     | 6632     | 7201     | 3000     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 13 | 1×BBB + 5×BBS + 1×BSS         | 3024     | 3386     | 3591     | 3684     | 3903     | 4004     | 3836     | 3600     | 3024     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 14 | 1×BBB + 5×BBS + 1×BSS         | 3024     | 3386     | 3591     | 3684     | 3903     | 4004     | 3836     | 3600     | 3024     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 15 | 2×BBB + 4×BBS + 1×SSS         | 3028     | 3463     | 3859     | 4125     | 4754     | 5213     | 5503     | 5633     | 3028     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 16 | 2×BBB + 3×BBS + 2×BSS         | 3038     | 3426     | 3752     | 3952     | 4544     | 4942     | 5257     | 5408     | 3038     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 17 | 2×BBB + 3×BBS + 2×BSS         | 3038     | 3426     | 3752     | 3952     | 4544     | 4942     | 5257     | 5408     | 3038     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 18 | 2×BBB + 3×BBS + 2×BSS         | 3038     | 3426     | 3752     | 3952     | 4544     | 4942     | 5257     | 5408     | 3038     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 19 | 3×BBB + 2×BBS + 1×BSS + 1×SSS | 3042     | 3503     | 4020     | 4393     | 5395     | 6151     | 6924     | 7441     | 3042     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 20 | 3×BBB + 2×BBS + 1×BSS + 1×SSS | 3042     | 3503     | 4020     | 4393     | 5395     | 6151     | 6924     | 7441     | 3042     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 21 | 3×BBB + 1×BBS + 3×BSS         | 3052     | 3466     | 3913     | 4220     | 5185     | 5880     | 6678     | 7216     | 3052     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 22 | 3×BBB + 1×BBS + 3×BSS         | 3052     | 3466     | 3913     | 4220     | 5185     | 5880     | 6678     | 7216     | 3052     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 23 | 4×BBB + 2×BSS + 1×SSS         | 3056     | 3543     | 4181     | 4661     | 6036     | 7089     | 8345     | 9249     | 3056     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 24 | 1×BBB + 6×BBS                 | 3066     | 3519     | 3813     | 3968     | 4217     | 4347     | 4128     | 3840     | 3066     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 25 | 2×BBB + 4×BBS + 1×BSS         | 3080     | 3559     | 3974     | 4236     | 4858     | 5285     | 5549     | 5648     | 3080     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 26 | 2×BBB + 4×BBS + 1×BSS         | 3080     | 3559     | 3974     | 4236     | 4858     | 5285     | 5549     | 5648     | 3080     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 27 | 3×BBB + 3×BBS + 1×SSS         | 3084     | 3636     | 4242     | 4677     | 5709     | 6494     | 7216     | 7681     | 3084     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 28 | 3×BBB + 2×BBS + 2×BSS         | 3094     | 3599     | 4135     | 4504     | 5499     | 6223     | 6970     | 7456     | 3094     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 29 | 3×BBB + 2×BBS + 2×BSS         | 3094     | 3599     | 4135     | 4504     | 5499     | 6223     | 6970     | 7456     | 3094     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 30 | 3×BBB + 2×BBS + 2×BSS         | 3094     | 3599     | 4135     | 4504     | 5499     | 6223     | 6970     | 7456     | 3094     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 31 | 4×BBB + 1×BBS + 1×BSS + 1×SSS | 3098     | 3676     | 4403     | 4945     | 6350     | 7432     | 8637     | 9489     | 3098     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 32 | 2×BBB + 5×BBS                 | 3122     | 3692     | 4196     | 4520     | 5172     | 5628     | 5841     | 5888     | 3122     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 33 | 3×BBB + 3×BBS + 1×BSS         | 3136     | 3732     | 4357     | 4788     | 5813     | 6566     | 7262     | 7696     | 3136     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 34 | 4×BBB + 1×BBS + 2×BSS         | 3150     | 3772     | 4518     | 5056     | 6454     | 7504     | 8683     | 9504     | 3150     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 35 | 3×BBB + 4×BBS                 | 3178     | 3865     | 4579     | 5072     | 6127     | 6909     | 7554     | 7936     | 3178     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 36 | 4×BBB + 2×BBS + 1×BSS         | 3192     | 3905     | 4740     | 5340     | 6768     | 7847     | 8975     | 9744     | 3192     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 37 | 5×BBB + 2×BSS                 | 3206     | 3945     | 4901     | 5608     | 7409     | 8785     | 10396    | 11552    | 3206     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 38 | 4×BBB + 3×BBS                 | 3234     | 4038     | 4962     | 5624     | 7082     | 8190     | 9267     | 9984     | 3234     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 39 | 5×BBB + 2×BBS                 | 3290     | 4211     | 5345     | 6176     | 8037     | 9471     | 10980    | 12032    | 3290     | 9+8 / 9+8 / 9+8        | —       | —                                                              |
| 40 | 7×BBB                         | 3402     | 4557     | 6111     | 7280     | 9947     | 12033    | 14406    | 16128    | 3402     | 9+8 / 9+8 / 9+8        | —       | —                                                              |

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

_Generated by `RecombinationMultisetReport` (target ⟨17,17,17⟩, base ⟨2,2,2⟩)._
