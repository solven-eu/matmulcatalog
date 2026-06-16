# Candidate outer bases ranked by implied ω

Filter: leafOnly=true, ncOnly=true, maxDim≤12.  ω = 3·log(r)/log(n·m·p).

Allocation columns are evaluated at the reference cubic target ⟨17,17,17⟩.
  - **raw**:  product of per-axis compositions, no constraints.
  - **imb≤3**:  filtered by per-axis max−min ≤ 3.
  - **flip-canon**:  + per-axis allocation kept only if lex-smaller than its reverse (safe iff axis-flip orbit of base is covered by mask sweep or pool expansion).
  - **+perm-canon**:  + cubic axis-permutation orbit canonicalised on (aA,aB,aC); applies to cubic bases on the cubic target only. `--` means enumeration capped.

Source: regenerate via `eu.solven.matmul.docs.explore.RankBasesByOmega`.

| ω | Shape | Rank | Cubic? | Field | Source | raw | imb≤3 | flip-canon | +perm-canon |
|---|---|---|---|---|---|---:|---:|---:|---:|
| 2,7743 | ⟨3,3,6⟩ | 40 | no | R/Q/Z | Smirnov 2013 | 62,9M | 78,6k | 13,4k | 13,4k |
| 2,7773 | ⟨4,4,4⟩ | 47 | yes | F2 | AlphaTensor 2022 | 175,6M | 85,2k | 10,6k | 2,0k |
| 2,7856 | ⟨6,6,12⟩ | 280 | no | R/Q/Z | Fmm-lille | 83,3G | 1,3G | 157,4M | 157,4M |
| 2,7896 | ⟨9,9,12⟩ | 600 | no | R/Q/Z | Fmm-lille | 723,5G | 301,7G | 38,2G | 38,2G |
| 2,7925 | ⟨4,4,4⟩ | 48 | yes | C | Alphaevolve | 175,6M | 85,2k | 10,6k | 2,0k |
| 2,7925 | ⟨4,4,4⟩ | 48 | yes | R/Q/Z | Dumas-pernet-sedoglavic-2025 | 175,6M | 85,2k | 10,6k | 2,0k |
| 2,7937 | ⟨9,10,12⟩ | 668 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,7957 | ⟨12,12,12⟩ | 1040 | yes | R/Q/Z | Fmm-lille | 83,3G | 75,4G | 9,4G | -- |
| 2,7962 | ⟨9,10,12⟩ | 672 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,7972 | ⟨5,5,12⟩ | 204 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,7974 | ⟨8,8,8⟩ | 336 | yes | R/Q/Z | Solven-strassen-2026 | 1497,2G | 151,2G | 18,9G | -- |
| 2,7981 | ⟨9,12,12⟩ | 800 | no | R/Q/Z | Fmm-lille | 245,6G | 150,8G | 19,0G | 19,0G |
| 2,7982 | ⟨3,4,6⟩ | 54 | no | R/Q/Z | Alphaevolve | 293,5M | 288,3k | 42,0k | 42,0k |
| 2,7983 | ⟨9,11,12⟩ | 738 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 259,9G | 32,9G | 32,9G |
| 2,7988 | ⟨9,10,12⟩ | 676 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,7995 | ⟨11,12,12⟩ | 968 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 129,9G | 16,4G | 16,4G |
| 2,8000 | ⟨8,12,12⟩ | 720 | no | R/Q/Z | Fmm-lille | 218,3G | 95,1G | 11,9G | 11,9G |
| 2,8006 | ⟨9,11,12⟩ | 742 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 259,9G | 32,9G | 32,9G |
| 2,8012 | ⟨11,12,12⟩ | 972 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 129,9G | 16,4G | 16,4G |
| 2,8012 | ⟨6,8,12⟩ | 378 | no | R/Q/Z | Perminov (tensor decomposition) | 218,3G | 12,3G | 1,5G | 1,5G |
| 2,8023 | ⟨5,5,12⟩ | 206 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8028 | ⟨11,12,12⟩ | 976 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 129,9G | 16,4G | 16,4G |
| 2,8033 | ⟨9,12,12⟩ | 810 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 245,6G | 150,8G | 19,0G | 19,0G |
| 2,8036 | ⟨4,4,9⟩ | 104 | no | R/Q/Z | Fmm-lille | 4,0G | 16,4M | 2,1M | 2,1M |
| 2,8042 | ⟨6,6,11⟩ | 268 | no | R/Q/Z | Fmm-lille | 152,8G | 2,2G | 273,2M | 273,2M |
| 2,8045 | ⟨11,12,12⟩ | 980 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 129,9G | 16,4G | 16,4G |
| 2,8045 | ⟨6,7,12⟩ | 336 | no | R/Q/Z | Fmm-lille | 152,8G | 4,9G | 622,7M | 622,7M |
| 2,8048 | ⟨4,4,10⟩ | 115 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 3,6G | 17,9M | 2,2M | 2,2M |
| 2,8048 | ⟨5,5,12⟩ | 207 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8052 | ⟨3,4,7⟩ | 63 | no | C | Alphaevolve | 538,1M | 1,1M | 166,3k | 166,3k |
| 2,8052 | ⟨3,4,7⟩ | 63 | no | R/Q/Z | alphaevolve-2025 (over C); dumas-pernet-sedoglavic-2025 (this scheme, over Q without sqrt(-1)) | 538,1M | 1,1M | 166,3k | 166,3k |
| 2,8064 | ⟨12,12,12⟩ | 1068 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 83,3G | 75,4G | 9,4G | -- |
| 2,8070 | ⟨10,12,12⟩ | 902 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 164,9G | 20,6G | 20,6G |
| 2,8073 | ⟨9,9,11⟩ | 576 | no | R/Q/Z | AlphaTensor 2022 | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8074 | ⟨2,2,2⟩ | 7 | yes | F2 | Strassen 1969 | 4,1k | 64 | 8 | 4 |
| 2,8074 | ⟨2,2,2⟩ | 7 | yes | R/Q/Z | Strassen | 4,1k | 64 | 8 | 4 |
| 2,8074 | ⟨4,4,4⟩ | 49 | yes | R/Q/Z | AlphaTensor 2022 | 175,6M | 85,2k | 10,6k | 2,0k |
| 2,8074 | ⟨8,8,8⟩ | 343 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 1497,2G | 151,2G | 18,9G | -- |
| 2,8074 | ⟨9,11,12⟩ | 754 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 259,9G | 32,9G | 32,9G |
| 2,8074 | ⟨5,5,12⟩ | 208 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8075 | ⟨12,12,12⟩ | 1071 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 83,3G | 75,4G | 9,4G | -- |
| 2,8075 | ⟨6,6,6⟩ | 153 | yes | R/Q/Z | Moosbauer (symmetric flips) | 83,3G | 162,8M | 20,3M | -- |
| 2,8076 | ⟨8,9,12⟩ | 560 | no | R/Q/Z | AlphaTensor 2022 | 643,1G | 190,2G | 23,9G | 23,9G |
| 2,8076 | ⟨6,12,12⟩ | 560 | no | R/Q/Z | fmm-lille | 83,3G | 9,7G | 1,2G | 1,2G |
| 2,8077 | ⟨3,6,6⟩ | 80 | no | R/Q/Z | Fmm-lille | 2,3G | 3,6M | 521,7k | 521,7k |
| 2,8077 | ⟨3,3,12⟩ | 80 | no | R/Q/Z | Fmm-lille | 62,9M | 608,3k | 103,5k | 103,5k |
| 2,8078 | ⟨8,11,12⟩ | 676 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 163,9G | 20,6G | 20,6G |
| 2,8080 | ⟨6,6,10⟩ | 247 | no | R/Q/Z | Fmm-lille | 218,3G | 2,8G | 344,3M | 344,3M |
| 2,8081 | ⟨9,9,12⟩ | 626 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 301,7G | 38,2G | 38,2G |
| 2,8086 | ⟨11,12,12⟩ | 990 | no | R/Q/Z | AlphaTensor 2022 | 152,8G | 129,9G | 16,4G | 16,4G |
| 2,8088 | ⟨8,12,12⟩ | 735 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 95,1G | 11,9G | 11,9G |
| 2,8089 | ⟨10,12,12⟩ | 906 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 164,9G | 20,6G | 20,6G |
| 2,8089 | ⟨6,9,12⟩ | 429 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 245,6G | 19,5G | 2,5G | 2,5G |
| 2,8093 | ⟨4,6,6⟩ | 105 | no | R/Q/Z | Perminov (tensor decomposition) | 10,7G | 13,1M | 1,6M | 1,6M |
| 2,8094 | ⟨4,8,8⟩ | 180 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 73,3G | 1,2G | 156,1M | 156,1M |
| 2,8097 | ⟨6,6,12⟩ | 294 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 83,3G | 1,3G | 157,4M | 157,4M |
| 2,8098 | ⟨8,8,12⟩ | 504 | no | R/Q/Z | Fmm-lille | 571,7G | 119,9G | 15,0G | 15,0G |
| 2,8100 | ⟨8,9,9⟩ | 430 | no | R/Q/Z | Fmm-lille | 1894,9G | 380,5G | 48,1G | 48,1G |
| 2,8101 | ⟨10,10,12⟩ | 766 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 360,6G | 45,1G | 45,1G |
| 2,8107 | ⟨10,12,12⟩ | 910 | no | R/Q/Z | Perminov (tensor decomposition) | 218,3G | 164,9G | 20,6G | 20,6G |
| 2,8107 | ⟨9,11,12⟩ | 760 | no | R/Q/Z | AlphaTensor 2022 | 450,2G | 259,9G | 32,9G | 32,9G |
| 2,8108 | ⟨6,12,12⟩ | 564 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 83,3G | 9,7G | 1,2G | 1,2G |
| 2,8108 | ⟨8,9,12⟩ | 564 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 190,2G | 23,9G | 23,9G |
| 2,8108 | ⟨2,3,3⟩ | 15 | no | F2 | Hopcroft-Kerr 1971 | 230,4k | 576 | 98 | 98 |
| 2,8108 | ⟨2,3,3⟩ | 15 | no | R/Q/Z | AlphaTensor 2022 | 230,4k | 576 | 98 | 98 |
| 2,8108 | ⟨6,6,9⟩ | 225 | no | R/Q/Z | Perminov (tensor decomposition) | 245,6G | 2,5G | 316,8M | 316,8M |
| 2,8109 | ⟨9,9,12⟩ | 630 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 301,7G | 38,2G | 38,2G |
| 2,8110 | ⟨3,3,6⟩ | 42 | no | R/Q/Z | Moosbauer (symmetric flips) | 62,9M | 78,6k | 13,4k | 13,4k |
| 2,8111 | ⟨3,3,8⟩ | 55 | no | R/Q/Z | Solven-closure-2026 | 164,7M | 767,2k | 130,5k | 130,5k |
| 2,8112 | ⟨10,10,12⟩ | 768 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 360,6G | 45,1G | 45,1G |
| 2,8112 | ⟨5,6,7⟩ | 150 | no | R/Q/Z | meta-flip-graph search | 63,7G | 186,5M | 24,8M | 24,8M |
| 2,8112 | ⟨5,6,8⟩ | 170 | no | R/Q/Z | meta-flip-graph search | 90,9G | 465,5M | 61,1M | 61,1M |
| 2,8113 | ⟨9,10,12⟩ | 696 | no | R/Q/Z | AlphaTensor 2022 | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,8113 | ⟨6,10,12⟩ | 476 | no | R/Q/Z | Perminov (tensor decomposition) | 218,3G | 21,3G | 2,7G | 2,7G |
| 2,8118 | ⟨8,9,11⟩ | 521 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8118 | ⟨6,11,12⟩ | 521 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 16,8G | 2,1G | 2,1G |
| 2,8118 | ⟨8,10,12⟩ | 624 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 208,0G | 26,0G | 26,0G |
| 2,8118 | ⟨9,11,12⟩ | 762 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 259,9G | 32,9G | 32,9G |
| 2,8120 | ⟨5,6,6⟩ | 130 | no | R/Q/Z | meta-flip-graph search | 34,7G | 47,7M | 6,3M | 6,3M |
| 2,8121 | ⟨8,9,9⟩ | 432 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 380,5G | 48,1G | 48,1G |
| 2,8123 | ⟨10,10,12⟩ | 770 | no | R/Q/Z | Perminov (tensor decomposition) | 571,7G | 360,6G | 45,1G | 45,1G |
| 2,8127 | ⟨7,8,8⟩ | 306 | no | R/Q/Z | Fmm-lille | 1048,0G | 60,6G | 7,7G | 7,7G |
| 2,8127 | ⟨6,8,11⟩ | 357 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 21,2G | 2,7G | 2,7G |
| 2,8129 | ⟨11,11,12⟩ | 922 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 224,0G | 28,4G | 28,4G |
| 2,8131 | ⟨6,7,12⟩ | 342 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 4,9G | 622,7M | 622,7M |
| 2,8131 | ⟨6,8,9⟩ | 296 | no | R/Q/Z | Fmm-lille | 643,1G | 24,6G | 3,1G | 3,1G |
| 2,8131 | ⟨3,5,6⟩ | 68 | no | R/Q/Z | Alphaevolve | 954,0M | 1,0M | 160,5k | 160,5k |
| 2,8132 | ⟨6,9,12⟩ | 433 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 245,6G | 19,5G | 2,5G | 2,5G |
| 2,8134 | ⟨9,9,10⟩ | 534 | no | R/Q/Z | AlphaTensor 2022 | 1894,9G | 659,9G | 83,5G | 83,5G |
| 2,8134 | ⟨8,8,12⟩ | 508 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 119,9G | 15,0G | 15,0G |
| 2,8135 | ⟨6,8,10⟩ | 327 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 26,9G | 3,4G | 3,4G |
| 2,8135 | ⟨8,8,9⟩ | 388 | no | R/Q/Z | Fmm-lille | 1684,3G | 239,9G | 30,2G | 30,2G |
| 2,8135 | ⟨8,10,10⟩ | 528 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1497,2G | 454,9G | 56,9G | 56,9G |
| 2,8135 | ⟨11,12,12⟩ | 1002 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 129,9G | 16,4G | 16,4G |
| 2,8136 | ⟨10,10,10⟩ | 651 | yes | R/Q/Z | Perminov (tensor decomposition) | 1497,2G | 788,9G | 98,6G | -- |
| 2,8138 | ⟨5,5,7⟩ | 127 | no | R/Q/Z | meta-flip-graph search | 26,5G | 54,7M | 7,6M | 7,6M |
| 2,8140 | ⟨5,5,8⟩ | 144 | no | R/Q/Z | meta-flip-graph search | 37,9G | 136,4M | 18,8M | 18,8M |
| 2,8142 | ⟨4,9,10⟩ | 250 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 82,5G | 3,4G | 432,1M | 432,1M |
| 2,8142 | ⟨5,6,12⟩ | 250 | no | R/Q/Z | Fmm-lille | 34,7G | 369,0M | 48,4M | 48,4M |
| 2,8142 | ⟨6,9,12⟩ | 434 | no | R/Q/Z | Perminov (tensor decomposition) | 245,6G | 19,5G | 2,5G | 2,5G |
| 2,8143 | ⟨5,5,6⟩ | 110 | no | R/Q/Z | meta-flip-graph search | 14,5G | 14,0M | 1,9M | 1,9M |
| 2,8144 | ⟨4,4,5⟩ | 61 | no | R/Q/Z | Alphaevolve | 570,8M | 309,8k | 40,7k | 40,7k |
| 2,8144 | ⟨4,4,5⟩ | 61 | no | R/Q/Z | AlphaEvolve 2025 | 570,8M | 309,8k | 40,7k | 40,7k |
| 2,8146 | ⟨8,10,12⟩ | 628 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 208,0G | 26,0G | 26,0G |
| 2,8146 | ⟨5,5,9⟩ | 161 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 42,6G | 216,3M | 30,0M | 30,0M |
| 2,8147 | ⟨6,6,8⟩ | 203 | no | R/Q/Z | Perminov (tensor decomposition) | 218,3G | 1,6G | 198,5M | 198,5M |
| 2,8147 | ⟨4,12,12⟩ | 389 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 10,7G | 785,1M | 98,1M | 98,1M |
| 2,8149 | ⟨5,5,12⟩ | 211 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8149 | ⟨6,8,8⟩ | 266 | no | R/Q/Z | Perminov (tensor decomposition) | 571,7G | 15,5G | 1,9G | 1,9G |
| 2,8150 | ⟨9,10,12⟩ | 702 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,8152 | ⟨6,9,9⟩ | 332 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 39,0G | 4,9G | 4,9G |
| 2,8152 | ⟨8,10,10⟩ | 530 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1497,2G | 454,9G | 56,9G | 56,9G |
| 2,8153 | ⟨8,9,9⟩ | 435 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 380,5G | 48,1G | 48,1G |
| 2,8153 | ⟨6,9,12⟩ | 435 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 245,6G | 19,5G | 2,5G | 2,5G |
| 2,8154 | ⟨4,8,8⟩ | 182 | no | R/Q/Z | Perminov (tensor decomposition) | 73,3G | 1,2G | 156,1M | 156,1M |
| 2,8154 | ⟨5,5,10⟩ | 178 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 37,9G | 236,5M | 32,6M | 32,6M |
| 2,8155 | ⟨6,12,12⟩ | 570 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 83,3G | 9,7G | 1,2G | 1,2G |
| 2,8155 | ⟨8,9,12⟩ | 570 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 190,2G | 23,9G | 23,9G |
| 2,8155 | ⟨9,9,9⟩ | 486 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 2131,7G | 603,6G | 76,8G | -- |
| 2,8157 | ⟨10,11,12⟩ | 849 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 284,2G | 35,8G | 35,8G |
| 2,8158 | ⟨6,9,10⟩ | 367 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 42,6G | 5,4G | 5,4G |
| 2,8159 | ⟨9,9,10⟩ | 537 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 659,9G | 83,5G | 83,5G |
| 2,8159 | ⟨3,4,7⟩ | 64 | no | R/Q/Z | meta-flip-graph search | 538,1M | 1,1M | 166,3k | 166,3k |
| 2,8160 | ⟨8,10,12⟩ | 630 | no | R/Q/Z | Perminov (tensor decomposition) | 571,7G | 208,0G | 26,0G | 26,0G |
| 2,8160 | ⟨8,8,12⟩ | 511 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 119,9G | 15,0G | 15,0G |
| 2,8162 | ⟨4,9,10⟩ | 251 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 82,5G | 3,4G | 432,1M | 432,1M |
| 2,8162 | ⟨10,11,12⟩ | 850 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 284,2G | 35,8G | 35,8G |
| 2,8163 | ⟨5,5,5⟩ | 93 | yes | R/Q/Z | Alphaevolve | 6,0G | 4,1M | 592,7k | 102,3k |
| 2,8163 | ⟨7,12,12⟩ | 660 | no | R/Q/Z | Fmm-lille | 152,8G | 38,1G | 4,8G | 4,8G |
| 2,8164 | ⟨5,5,11⟩ | 195 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 26,5G | 186,4M | 25,9M | 25,9M |
| 2,8165 | ⟨6,8,10⟩ | 329 | no | R/Q/Z | AlphaTensor 2022 | 571,7G | 26,9G | 3,4G | 3,4G |
| 2,8165 | ⟨6,8,10⟩ | 329 | no | R/Q/Z | Perminov (tensor decomposition) | 571,7G | 26,9G | 3,4G | 3,4G |
| 2,8166 | ⟨3,5,7⟩ | 79 | no | R/Q/Z | meta-flip-graph search | 1,7G | 4,1M | 635,0k | 635,0k |
| 2,8166 | ⟨8,11,12⟩ | 690 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 163,9G | 20,6G | 20,6G |
| 2,8168 | ⟨9,10,12⟩ | 705 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,8169 | ⟨8,9,11⟩ | 527 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8169 | ⟨8,10,10⟩ | 532 | no | R/Q/Z | AlphaTensor 2022 | 1497,2G | 454,9G | 56,9G | 56,9G |
| 2,8169 | ⟨4,7,8⟩ | 161 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 500,5M | 63,3M | 63,3M |
| 2,8170 | ⟨7,10,12⟩ | 557 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 83,3G | 10,5G | 10,5G |
| 2,8170 | ⟨8,9,10⟩ | 482 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 416,0G | 52,3G | 52,3G |
| 2,8171 | ⟨4,9,11⟩ | 275 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 57,7G | 2,7G | 342,9M | 342,9M |
| 2,8171 | ⟨6,9,10⟩ | 368 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 42,6G | 5,4G | 5,4G |
| 2,8172 | ⟨8,8,9⟩ | 391 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 239,9G | 30,2G | 30,2G |
| 2,8173 | ⟨7,8,12⟩ | 452 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 48,0G | 6,1G | 6,1G |
| 2,8174 | ⟨6,7,11⟩ | 318 | no | R/Q/Z | Fmm-lille | 280,1G | 8,5G | 1,1G | 1,1G |
| 2,8175 | ⟨4,6,7⟩ | 123 | no | R/Q/Z | meta-flip-graph search | 19,6G | 51,3M | 6,5M | 6,5M |
| 2,8177 | ⟨9,10,11⟩ | 651 | no | R/Q/Z | Fmm-lille | 1179,0G | 568,6G | 72,0G | 72,0G |
| 2,8178 | ⟨10,10,11⟩ | 719 | no | R/Q/Z | Perminov (tensor decomposition) | 1048,0G | 621,7G | 78,2G | 78,2G |
| 2,8180 | ⟨3,3,7⟩ | 49 | no | R/Q/Z | Perminov (tensor decomposition) | 115,3M | 307,4k | 52,9k | 52,9k |
| 2,8182 | ⟨4,9,10⟩ | 252 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 82,5G | 3,4G | 432,1M | 432,1M |
| 2,8182 | ⟨6,6,10⟩ | 252 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 2,8G | 344,3M | 344,3M |
| 2,8183 | ⟨3,6,7⟩ | 94 | no | R/Q/Z | Fmm-lille | 4,2G | 14,0M | 2,1M | 2,1M |
| 2,8184 | ⟨9,9,10⟩ | 540 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 659,9G | 83,5G | 83,5G |
| 2,8185 | ⟨2,4,5⟩ | 32 | no | R/Q/Z | Alphaevolve | 16,3M | 28,2k | 3,7k | 3,7k |
| 2,8186 | ⟨6,7,9⟩ | 264 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 9,9G | 1,3G | 1,3G |
| 2,8188 | ⟨10,12,12⟩ | 928 | no | R/Q/Z | AlphaTensor 2022 | 218,3G | 164,9G | 20,6G | 20,6G |
| 2,8189 | ⟨6,10,11⟩ | 446 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 36,7G | 4,6G | 4,6G |
| 2,8189 | ⟨4,9,11⟩ | 276 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 57,7G | 2,7G | 342,9M | 342,9M |
| 2,8189 | ⟨6,6,11⟩ | 276 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 2,2G | 273,2M | 273,2M |
| 2,8189 | ⟨6,9,11⟩ | 404 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 33,6G | 4,3G | 4,3G |
| 2,8190 | ⟨3,4,4⟩ | 38 | no | F2 | Smirnov 2013 | 37,6M | 23,2k | 3,4k | 3,4k |
| 2,8190 | ⟨3,4,4⟩ | 38 | no | R/Q/Z | AlphaTensor 2022 | 37,6M | 23,2k | 3,4k | 3,4k |
| 2,8190 | ⟨9,10,10⟩ | 597 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 721,5G | 90,7G | 90,7G |
| 2,8190 | ⟨3,3,4⟩ | 29 | no | F2 | Smirnov 2013 | 8,1M | 6,3k | 1,1k | 1,1k |
| 2,8190 | ⟨3,3,4⟩ | 29 | no | R/Q/Z | AlphaTensor 2022 | 8,1M | 6,3k | 1,1k | 1,1k |
| 2,8190 | ⟨7,8,8⟩ | 310 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1048,0G | 60,6G | 7,7G | 7,7G |
| 2,8191 | ⟨11,11,12⟩ | 936 | no | R/Q/Z | Perminov (tensor decomposition) | 280,1G | 224,0G | 28,4G | 28,4G |
| 2,8191 | ⟨5,10,12⟩ | 408 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 6,2G | 819,6M | 819,6M |
| 2,8193 | ⟨7,8,12⟩ | 454 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 48,0G | 6,1G | 6,1G |
| 2,8193 | ⟨7,10,12⟩ | 560 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 83,3G | 10,5G | 10,5G |
| 2,8194 | ⟨4,11,12⟩ | 362 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 19,6G | 1,4G | 170,3M | 170,3M |
| 2,8195 | ⟨9,11,11⟩ | 715 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 448,1G | 57,1G | 57,1G |
| 2,8195 | ⟨4,5,7⟩ | 104 | no | R/Q/Z | meta-flip-graph search | 8,2G | 15,0M | 2,0M | 2,0M |
| 2,8196 | ⟨5,7,7⟩ | 176 | no | R/Q/Z | meta-flip-graph search | 116,7G | 729,3M | 98,0M | 98,0M |
| 2,8197 | ⟨9,10,10⟩ | 598 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 721,5G | 90,7G | 90,7G |
| 2,8197 | ⟨4,5,6⟩ | 90 | no | R/Q/Z | Alphaevolve | 4,5G | 3,8M | 504,5k | 504,5k |
| 2,8197 | ⟨3,5,8⟩ | 90 | no | R/Q/Z | meta-flip-graph search | 2,5G | 10,2M | 1,6M | 1,6M |
| 2,8197 | ⟨4,9,12⟩ | 300 | no | R/Q/Z | Hopcroft-Kerr 1971 | 31,5G | 1,6G | 197,5M | 197,5M |
| 2,8197 | ⟨4,4,11⟩ | 129 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 14,1M | 1,8M | 1,8M |
| 2,8198 | ⟨4,6,8⟩ | 140 | no | R/Q/Z | Perminov (tensor decomposition) | 28,0G | 128,0M | 16,0M | 16,0M |
| 2,8198 | ⟨9,10,12⟩ | 710 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,8199 | ⟨5,5,12⟩ | 213 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8200 | ⟨8,8,11⟩ | 475 | no | R/Q/Z | Perminov (tensor decomposition) | 1048,0G | 206,7G | 26,0G | 26,0G |
| 2,8200 | ⟨4,4,6⟩ | 73 | no | R/Q/Z | meta-flip-graph search | 1,4G | 1,1M | 132,1k | 132,1k |
| 2,8200 | ⟨3,4,8⟩ | 73 | no | R/Q/Z | meta-flip-graph search | 768,8M | 2,8M | 410,3k | 410,3k |
| 2,8200 | ⟨4,5,8⟩ | 118 | no | R/Q/Z | meta-flip-graph search | 11,7G | 37,5M | 4,9M | 4,9M |
| 2,8200 | ⟨4,4,10⟩ | 118 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 3,6G | 17,9M | 2,2M | 2,2M |
| 2,8200 | ⟨7,8,9⟩ | 347 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 96,1G | 12,2G | 12,2G |
| 2,8201 | ⟨7,9,12⟩ | 508 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8201 | ⟨5,6,9⟩ | 193 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,3G | 738,3M | 97,5M | 97,5M |
| 2,8201 | ⟨8,11,11⟩ | 641 | no | R/Q/Z | Fmm-lille | 733,6G | 282,5G | 35,8G | 35,8G |
| 2,8203 | ⟨2,4,4⟩ | 26 | no | F2 | Hopcroft-Kerr 1971 | 5,0M | 7,7k | 968 | 968 |
| 2,8203 | ⟨2,4,4⟩ | 26 | no | R/Q/Z | AlphaTensor 2022 | 5,0M | 7,7k | 968 | 968 |
| 2,8203 | ⟨2,4,4⟩ | 26 | no | R/Q/Z | Hopcroft-Kerr 1971 | 5,0M | 7,7k | 968 | 968 |
| 2,8203 | ⟨8,9,11⟩ | 531 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8204 | ⟨4,7,8⟩ | 162 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 500,5M | 63,3M | 63,3M |
| 2,8204 | ⟨9,10,10⟩ | 599 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 721,5G | 90,7G | 90,7G |
| 2,8207 | ⟨4,4,9⟩ | 107 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,0G | 16,4M | 2,1M | 2,1M |
| 2,8208 | ⟨4,5,9⟩ | 132 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 13,1G | 59,5M | 7,9M | 7,9M |
| 2,8209 | ⟨9,9,11⟩ | 594 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8210 | ⟨6,11,11⟩ | 490 | no | R/Q/Z | Fmm-lille | 280,1G | 29,0G | 3,7G | 3,7G |
| 2,8210 | ⟨6,9,10⟩ | 371 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 42,6G | 5,4G | 5,4G |
| 2,8211 | ⟨3,4,5⟩ | 47 | no | F2 | AlphaTensor 2022 | 122,3M | 84,5k | 12,9k | 12,9k |
| 2,8211 | ⟨3,4,5⟩ | 47 | no | R/Q/Z | AlphaTensor 2022 | 122,3M | 84,5k | 12,9k | 12,9k |
| 2,8211 | ⟨2,5,6⟩ | 47 | no | R/Q/Z | Alphaevolve | 127,2M | 349,4k | 45,9k | 45,9k |
| 2,8211 | ⟨7,7,12⟩ | 402 | no | R/Q/Z | Fmm-lille | 280,1G | 19,3G | 2,5G | 2,5G |
| 2,8211 | ⟨8,9,11⟩ | 532 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8212 | ⟨6,7,10⟩ | 293 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 10,8G | 1,4G | 1,4G |
| 2,8212 | ⟨9,10,10⟩ | 600 | no | R/Q/Z | Perminov (tensor decomposition) | 1684,3G | 721,5G | 90,7G | 90,7G |
| 2,8212 | ⟨4,5,5⟩ | 76 | no | F2 | AlphaTensor 2022 | 1,9G | 1,1M | 155,2k | 155,2k |
| 2,8212 | ⟨4,5,5⟩ | 76 | no | R/Q/Z | AlphaTensor 2022 | 1,9G | 1,1M | 155,2k | 155,2k |
| 2,8213 | ⟨11,11,12⟩ | 941 | no | R/Q/Z | AlphaTensor 2022 | 280,1G | 224,0G | 28,4G | 28,4G |
| 2,8214 | ⟨3,5,5⟩ | 58 | no | F2 | Sedoglavic-Smirnov 2021 | 397,5M | 307,2k | 49,4k | 49,4k |
| 2,8214 | ⟨3,5,5⟩ | 58 | no | R/Q/Z | AlphaTensor 2022 | 397,5M | 307,2k | 49,4k | 49,4k |
| 2,8214 | ⟨10,11,11⟩ | 793 | no | R/Q/Z | Perminov (tensor decomposition) | 733,6G | 490,0G | 62,1G | 62,1G |
| 2,8214 | ⟨7,8,9⟩ | 348 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 96,1G | 12,2G | 12,2G |
| 2,8214 | ⟨5,5,9⟩ | 163 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 42,6G | 216,3M | 30,0M | 30,0M |
| 2,8215 | ⟨5,5,10⟩ | 180 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 37,9G | 236,5M | 32,6M | 32,6M |
| 2,8216 | ⟨8,10,11⟩ | 588 | no | R/Q/Z | Perminov (tensor decomposition) | 1048,0G | 358,5G | 45,1G | 45,1G |
| 2,8217 | ⟨9,9,11⟩ | 595 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8217 | ⟨9,10,11⟩ | 657 | no | R/Q/Z | AlphaTensor 2022 | 1179,0G | 568,6G | 72,0G | 72,0G |
| 2,8218 | ⟨4,5,10⟩ | 146 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8218 | ⟨7,9,12⟩ | 510 | no | R/Q/Z | AlphaTensor 2022 | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8220 | ⟨8,9,11⟩ | 533 | no | R/Q/Z | AlphaTensor 2022 | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8220 | ⟨8,9,11⟩ | 533 | no | R/Q/Z | Perminov (tensor decomposition) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8221 | ⟨4,4,8⟩ | 96 | no | R/Q/Z | Perminov (tensor decomposition) | 3,6G | 10,3M | 1,3M | 1,3M |
| 2,8222 | ⟨6,7,8⟩ | 238 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 6,2G | 785,5M | 785,5M |
| 2,8222 | ⟨7,12,12⟩ | 669 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 38,1G | 4,8G | 4,8G |
| 2,8223 | ⟨5,8,12⟩ | 333 | no | R/Q/Z | Perminov (tensor decomposition) | 90,9G | 3,6G | 472,6M | 472,6M |
| 2,8223 | ⟨5,5,12⟩ | 214 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8224 | ⟨7,8,10⟩ | 385 | no | R/Q/Z | Perminov (tensor decomposition) | 1048,0G | 105,1G | 13,3G | 13,3G |
| 2,8224 | ⟨6,7,9⟩ | 266 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 9,9G | 1,3G | 1,3G |
| 2,8224 | ⟨9,9,11⟩ | 596 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8224 | ⟨8,9,12⟩ | 579 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 190,2G | 23,9G | 23,9G |
| 2,8224 | ⟨6,9,11⟩ | 407 | no | R/Q/Z | Perminov (tensor decomposition) | 450,2G | 33,6G | 4,3G | 4,3G |
| 2,8225 | ⟨7,10,12⟩ | 564 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 83,3G | 10,5G | 10,5G |
| 2,8227 | ⟨5,12,12⟩ | 488 | no | R/Q/Z | Fmm-lille | 34,7G | 2,9G | 374,7M | 374,7M |
| 2,8227 | ⟨6,9,12⟩ | 442 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 245,6G | 19,5G | 2,5G | 2,5G |
| 2,8228 | ⟨8,9,11⟩ | 534 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8228 | ⟨6,11,12⟩ | 534 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 16,8G | 2,1G | 2,1G |
| 2,8229 | ⟨3,3,10⟩ | 69 | no | R/Q/Z | Perminov (tensor decomposition) | 164,7M | 1,3M | 226,4k | 226,4k |
| 2,8229 | ⟨5,6,9⟩ | 194 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,3G | 738,3M | 97,5M | 97,5M |
| 2,8229 | ⟨4,5,11⟩ | 160 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 8,2G | 51,3M | 6,8M | 6,8M |
| 2,8230 | ⟨3,6,9⟩ | 120 | no | R/Q/Z | Fmm-lille | 6,7G | 55,4M | 8,1M | 8,1M |
| 2,8231 | ⟨9,11,11⟩ | 721 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 448,1G | 57,1G | 57,1G |
| 2,8233 | ⟨7,11,12⟩ | 618 | no | R/Q/Z | Fmm-lille | 280,1G | 65,7G | 8,4G | 8,4G |
| 2,8233 | ⟨6,8,11⟩ | 365 | no | R/Q/Z | AlphaTensor 2022 | 400,2G | 21,2G | 2,7G | 2,7G |
| 2,8233 | ⟨4,11,12⟩ | 365 | no | R/Q/Z | Perminov (tensor decomposition) | 19,6G | 1,4G | 170,3M | 170,3M |
| 2,8234 | ⟨7,7,12⟩ | 404 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 19,3G | 2,5G | 2,5G |
| 2,8235 | ⟨6,7,11⟩ | 322 | no | R/Q/Z | AlphaTensor 2022 | 280,1G | 8,5G | 1,1G | 1,1G |
| 2,8235 | ⟨6,7,11⟩ | 322 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 8,5G | 1,1G | 1,1G |
| 2,8235 | ⟨8,11,11⟩ | 646 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 282,5G | 35,8G | 35,8G |
| 2,8236 | ⟨6,9,10⟩ | 373 | no | R/Q/Z | AlphaTensor 2022 | 643,1G | 42,6G | 5,4G | 5,4G |
| 2,8236 | ⟨6,9,10⟩ | 373 | no | R/Q/Z | Perminov (tensor decomposition) | 643,1G | 42,6G | 5,4G | 5,4G |
| 2,8236 | ⟨8,9,10⟩ | 489 | no | R/Q/Z | AlphaTensor 2022 | 1684,3G | 416,0G | 52,3G | 52,3G |
| 2,8236 | ⟨6,10,12⟩ | 489 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 21,3G | 2,7G | 2,7G |
| 2,8237 | ⟨8,9,11⟩ | 535 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8237 | ⟨2,6,6⟩ | 56 | no | R/Q/Z | Moosbauer (symmetric flips) | 305,3M | 1,2M | 149,1k | 149,1k |
| 2,8237 | ⟨3,4,6⟩ | 56 | no | R/Q/Z | Moosbauer (symmetric flips) | 293,5M | 288,3k | 42,0k | 42,0k |
| 2,8237 | ⟨3,3,8⟩ | 56 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 164,7M | 767,2k | 130,5k | 130,5k |
| 2,8237 | ⟨5,10,12⟩ | 412 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 6,2G | 819,6M | 819,6M |
| 2,8238 | ⟨4,7,8⟩ | 163 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 500,5M | 63,3M | 63,3M |
| 2,8238 | ⟨4,4,12⟩ | 141 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,4G | 8,2M | 1,0M | 1,0M |
| 2,8240 | ⟨4,5,12⟩ | 174 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,5G | 29,7M | 3,9M | 3,9M |
| 2,8241 | ⟨11,11,11⟩ | 873 | yes | R/Q/Z | Perminov (tensor decomposition) | 513,5G | 386,1G | 49,3G | -- |
| 2,8241 | ⟨3,3,5⟩ | 36 | no | F2 | Smirnov 2013 | 26,2M | 23,0k | 4,1k | 4,1k |
| 2,8241 | ⟨3,3,5⟩ | 36 | no | R/Q/Z | AlphaTensor 2022 | 26,2M | 23,0k | 4,1k | 4,1k |
| 2,8242 | ⟨4,4,11⟩ | 130 | no | R/Q/Z | Perminov (tensor decomposition) | 2,5G | 14,1M | 1,8M | 1,8M |
| 2,8242 | ⟨4,9,10⟩ | 255 | no | R/Q/Z | AlphaTensor 2022 | 82,5G | 3,4G | 432,1M | 432,1M |
| 2,8243 | ⟨9,10,11⟩ | 661 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 568,6G | 72,0G | 72,0G |
| 2,8243 | ⟨6,7,8⟩ | 239 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 6,2G | 785,5M | 785,5M |
| 2,8244 | ⟨4,9,11⟩ | 279 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 57,7G | 2,7G | 342,9M | 342,9M |
| 2,8244 | ⟨7,8,11⟩ | 423 | no | R/Q/Z | Perminov (tensor decomposition) | 733,6G | 82,8G | 10,5G | 10,5G |
| 2,8245 | ⟨7,9,12⟩ | 513 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8245 | ⟨6,10,12⟩ | 490 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 21,3G | 2,7G | 2,7G |
| 2,8245 | ⟨5,5,10⟩ | 181 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 37,9G | 236,5M | 32,6M | 32,6M |
| 2,8245 | ⟨5,5,11⟩ | 198 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 26,5G | 186,4M | 25,9M | 25,9M |
| 2,8246 | ⟨4,4,7⟩ | 85 | no | R/Q/Z | Alphaevolve | 2,5G | 4,1M | 522,7k | 522,7k |
| 2,8246 | ⟨9,9,11⟩ | 599 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8246 | ⟨4,11,12⟩ | 366 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 19,6G | 1,4G | 170,3M | 170,3M |
| 2,8247 | ⟨3,5,7⟩ | 80 | no | R/Q/Z | Alphaevolve | 1,7G | 4,1M | 635,0k | 635,0k |
| 2,8248 | ⟨4,7,7⟩ | 144 | no | R/Q/Z | meta-flip-graph search | 35,9G | 200,6M | 25,7M | 25,7M |
| 2,8248 | ⟨5,10,12⟩ | 413 | no | R/Q/Z | Perminov (tensor decomposition) | 90,9G | 6,2G | 819,6M | 819,6M |
| 2,8250 | ⟨9,10,11⟩ | 662 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 568,6G | 72,0G | 72,0G |
| 2,8250 | ⟨4,4,10⟩ | 119 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 3,6G | 17,9M | 2,2M | 2,2M |
| 2,8253 | ⟨7,10,10⟩ | 478 | no | R/Q/Z | AlphaTensor 2022 | 1048,0G | 182,3G | 23,1G | 23,1G |
| 2,8254 | ⟨9,9,11⟩ | 600 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8255 | ⟨9,11,11⟩ | 725 | no | R/Q/Z | AlphaTensor 2022 | 825,3G | 448,1G | 57,1G | 57,1G |
| 2,8255 | ⟨7,9,10⟩ | 433 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 166,7G | 21,2G | 21,2G |
| 2,8255 | ⟨7,7,8⟩ | 277 | no | R/Q/Z | Fmm-lille | 733,6G | 24,3G | 3,1G | 3,1G |
| 2,8255 | ⟨8,11,11⟩ | 649 | no | R/Q/Z | AlphaTensor 2022 | 733,6G | 282,5G | 35,8G | 35,8G |
| 2,8256 | ⟨9,10,10⟩ | 606 | no | R/Q/Z | AlphaTensor 2022 | 1684,3G | 721,5G | 90,7G | 90,7G |
| 2,8256 | ⟨5,6,9⟩ | 195 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,3G | 738,3M | 97,5M | 97,5M |
| 2,8257 | ⟨4,5,10⟩ | 147 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8258 | ⟨8,9,9⟩ | 445 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 380,5G | 48,1G | 48,1G |
| 2,8259 | ⟨6,9,11⟩ | 410 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 33,6G | 4,3G | 4,3G |
| 2,8261 | ⟨4,9,11⟩ | 280 | no | R/Q/Z | AlphaTensor 2022 | 57,7G | 2,7G | 342,9M | 342,9M |
| 2,8261 | ⟨4,8,12⟩ | 272 | no | R/Q/Z | Fmm-lille | 28,0G | 990,2M | 123,8M | 123,8M |
| 2,8262 | ⟨6,7,9⟩ | 268 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 9,9G | 1,3G | 1,3G |
| 2,8262 | ⟨6,7,10⟩ | 296 | no | R/Q/Z | AlphaTensor 2022 | 400,2G | 10,8G | 1,4G | 1,4G |
| 2,8262 | ⟨6,7,10⟩ | 296 | no | R/Q/Z | Perminov (tensor decomposition) | 400,2G | 10,8G | 1,4G | 1,4G |
| 2,8262 | ⟨5,7,12⟩ | 296 | no | R/Q/Z | Fmm-lille | 63,7G | 1,4G | 191,6M | 191,6M |
| 2,8262 | ⟨7,9,12⟩ | 515 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8263 | ⟨4,5,11⟩ | 161 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 8,2G | 51,3M | 6,8M | 6,8M |
| 2,8263 | ⟨3,6,8⟩ | 108 | no | R/Q/Z | meta-flip-graph search | 6,0G | 34,9M | 5,1M | 5,1M |
| 2,8263 | ⟨4,4,9⟩ | 108 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,0G | 16,4M | 2,1M | 2,1M |
| 2,8263 | ⟨3,4,12⟩ | 108 | no | R/Q/Z | Perminov (tensor decomposition) | 293,5M | 2,2M | 325,2k | 325,2k |
| 2,8264 | ⟨6,6,7⟩ | 183 | no | R/Q/Z | meta-flip-graph search | 152,8G | 636,5M | 80,5M | 80,5M |
| 2,8265 | ⟨6,11,11⟩ | 496 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 29,0G | 3,7G | 3,7G |
| 2,8265 | ⟨9,9,12⟩ | 653 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 301,7G | 38,2G | 38,2G |
| 2,8265 | ⟨7,9,10⟩ | 434 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 166,7G | 21,2G | 21,2G |
| 2,8266 | ⟨5,6,11⟩ | 236 | no | R/Q/Z | Fmm-lille | 63,7G | 636,2M | 84,1M | 84,1M |
| 2,8266 | ⟨9,9,9⟩ | 498 | yes | R/Q/Z | AlphaTensor 2022 | 2131,7G | 603,6G | 76,8G | -- |
| 2,8269 | ⟨7,8,9⟩ | 352 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 96,1G | 12,2G | 12,2G |
| 2,8270 | ⟨6,9,11⟩ | 411 | no | R/Q/Z | AlphaTensor 2022 | 450,2G | 33,6G | 4,3G | 4,3G |
| 2,8270 | ⟨8,9,11⟩ | 539 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8271 | ⟨8,8,10⟩ | 441 | no | R/Q/Z | AlphaTensor 2022 | 1497,2G | 262,3G | 32,8G | 32,8G |
| 2,8271 | ⟨4,6,10⟩ | 175 | no | R/Q/Z | Perminov (tensor decomposition) | 28,0G | 222,0M | 27,7M | 27,7M |
| 2,8271 | ⟨4,5,12⟩ | 175 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,5G | 29,7M | 3,9M | 3,9M |
| 2,8271 | ⟨6,10,11⟩ | 454 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 36,7G | 4,6G | 4,6G |
| 2,8271 | ⟨5,11,12⟩ | 454 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 4,9G | 650,4M | 650,4M |
| 2,8271 | ⟨7,9,12⟩ | 516 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8272 | ⟨4,7,8⟩ | 164 | no | R/Q/Z | Perminov (tensor decomposition) | 51,3G | 500,5M | 63,3M | 63,3M |
| 2,8272 | ⟨5,6,10⟩ | 216 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 807,2M | 105,9M | 105,9M |
| 2,8272 | ⟨5,5,12⟩ | 216 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8272 | ⟨7,10,12⟩ | 570 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 83,3G | 10,5G | 10,5G |
| 2,8272 | ⟨9,11,11⟩ | 728 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 448,1G | 57,1G | 57,1G |
| 2,8274 | ⟨7,7,8⟩ | 278 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 24,3G | 3,1G | 3,1G |
| 2,8274 | ⟨10,10,12⟩ | 798 | no | R/Q/Z | AlphaTensor 2022 | 571,7G | 360,6G | 45,1G | 45,1G |
| 2,8274 | ⟨3,3,11⟩ | 76 | no | R/Q/Z | Perminov (tensor decomposition) | 115,3M | 1,0M | 179,6k | 179,6k |
| 2,8274 | ⟨6,7,7⟩ | 212 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 2,5G | 318,4M | 318,4M |
| 2,8274 | ⟨9,9,10⟩ | 551 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 659,9G | 83,5G | 83,5G |
| 2,8275 | ⟨7,11,12⟩ | 624 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 65,7G | 8,4G | 8,4G |
| 2,8276 | ⟨8,10,11⟩ | 596 | no | R/Q/Z | AlphaTensor 2022 | 1048,0G | 358,5G | 45,1G | 45,1G |
| 2,8279 | ⟨10,11,12⟩ | 874 | no | R/Q/Z | AlphaTensor 2022 | 400,2G | 284,2G | 35,8G | 35,8G |
| 2,8279 | ⟨2,3,4⟩ | 20 | no | F2 | Hopcroft-Kerr 1971 | 1,1M | 2,1k | 308 | 308 |
| 2,8279 | ⟨2,3,4⟩ | 20 | no | R/Q/Z | AlphaTensor 2022 | 1,1M | 2,1k | 308 | 308 |
| 2,8280 | ⟨7,10,11⟩ | 526 | no | R/Q/Z | Perminov (tensor decomposition) | 733,6G | 143,7G | 18,3G | 18,3G |
| 2,8280 | ⟨7,9,12⟩ | 517 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8282 | ⟨6,9,9⟩ | 341 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 39,0G | 4,9G | 4,9G |
| 2,8282 | ⟨7,7,11⟩ | 376 | no | R/Q/Z | Fmm-lille | 513,5G | 33,2G | 4,3G | 4,3G |
| 2,8282 | ⟨8,8,12⟩ | 525 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 119,9G | 15,0G | 15,0G |
| 2,8282 | ⟨5,10,12⟩ | 416 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 6,2G | 819,6M | 819,6M |
| 2,8282 | ⟨5,8,8⟩ | 230 | no | R/Q/Z | Perminov (tensor decomposition) | 238,2G | 4,5G | 596,1M | 596,1M |
| 2,8282 | ⟨4,8,10⟩ | 230 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 73,3G | 2,2G | 270,8M | 270,8M |
| 2,8284 | ⟨5,6,9⟩ | 196 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,3G | 738,3M | 97,5M | 97,5M |
| 2,8284 | ⟨3,3,9⟩ | 63 | no | R/Q/Z | Perminov (tensor decomposition) | 185,3M | 1,2M | 208,3k | 208,3k |
| 2,8286 | ⟨3,5,9⟩ | 102 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,8G | 16,2M | 2,5M | 2,5M |
| 2,8286 | ⟨4,11,11⟩ | 340 | no | R/Q/Z | Perminov (tensor decomposition) | 35,9G | 2,3G | 295,7M | 295,7M |
| 2,8287 | ⟨5,9,12⟩ | 377 | no | R/Q/Z | Perminov (tensor decomposition) | 102,3G | 5,7G | 754,2M | 754,2M |
| 2,8287 | ⟨4,4,11⟩ | 131 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 14,1M | 1,8M | 1,8M |
| 2,8287 | ⟨3,3,6⟩ | 43 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 62,9M | 78,6k | 13,4k | 13,4k |
| 2,8287 | ⟨4,7,7⟩ | 145 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 200,6M | 25,7M | 25,7M |
| 2,8289 | ⟨2,5,5⟩ | 40 | no | F2 | Hopcroft-Kerr 1971 | 53,0M | 102,4k | 14,1k | 14,1k |
| 2,8289 | ⟨2,5,5⟩ | 40 | no | R/Q/Z | AlphaTensor 2022 | 53,0M | 102,4k | 14,1k | 14,1k |
| 2,8289 | ⟨7,9,12⟩ | 518 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8289 | ⟨3,4,8⟩ | 74 | no | R/Q/Z | Alphaevolve | 768,8M | 2,8M | 410,3k | 410,3k |
| 2,8289 | ⟨7,11,12⟩ | 626 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 65,7G | 8,4G | 8,4G |
| 2,8290 | ⟨4,6,9⟩ | 159 | no | R/Q/Z | Perminov (tensor decomposition) | 31,5G | 203,0M | 25,5M | 25,5M |
| 2,8292 | ⟨7,7,8⟩ | 279 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 24,3G | 3,1G | 3,1G |
| 2,8292 | ⟨7,11,11⟩ | 577 | no | R/Q/Z | Fmm-lille | 513,5G | 113,2G | 14,5G | 14,5G |
| 2,8295 | ⟨4,5,10⟩ | 148 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8295 | ⟨4,5,9⟩ | 134 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 13,1G | 59,5M | 7,9M | 7,9M |
| 2,8295 | ⟨3,6,10⟩ | 134 | no | R/Q/Z | Fmm-lille | 6,0G | 60,5M | 8,8M | 8,8M |
| 2,8295 | ⟨6,6,12⟩ | 306 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 83,3G | 1,3G | 157,4M | 157,4M |
| 2,8296 | ⟨5,7,12⟩ | 298 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 1,4G | 191,6M | 191,6M |
| 2,8296 | ⟨6,9,9⟩ | 342 | no | R/Q/Z | AlphaTensor 2022 | 723,5G | 39,0G | 4,9G | 4,9G |
| 2,8296 | ⟨6,9,9⟩ | 342 | no | R/Q/Z | Perminov (tensor decomposition) | 723,5G | 39,0G | 4,9G | 4,9G |
| 2,8296 | ⟨5,6,10⟩ | 217 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 807,2M | 105,9M | 105,9M |
| 2,8296 | ⟨5,5,12⟩ | 217 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8297 | ⟨7,9,11⟩ | 478 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 131,4G | 16,8G | 16,8G |
| 2,8297 | ⟨7,8,9⟩ | 354 | no | R/Q/Z | AlphaTensor 2022 | 1179,0G | 96,1G | 12,2G | 12,2G |
| 2,8297 | ⟨2,5,7⟩ | 55 | no | R/Q/Z | meta-flip-graph search | 233,2M | 1,4M | 181,4k | 181,4k |
| 2,8298 | ⟨7,9,10⟩ | 437 | no | R/Q/Z | Perminov (tensor decomposition) | 1179,0G | 166,7G | 21,2G | 21,2G |
| 2,8299 | ⟨6,7,7⟩ | 213 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 2,5G | 318,4M | 318,4M |
| 2,8299 | ⟨5,5,11⟩ | 200 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 26,5G | 186,4M | 25,9M | 25,9M |
| 2,8299 | ⟨6,7,9⟩ | 270 | no | R/Q/Z | AlphaTensor 2022 | 450,2G | 9,9G | 1,3G | 1,3G |
| 2,8299 | ⟨4,4,10⟩ | 120 | no | R/Q/Z | Perminov (tensor decomposition) | 3,6G | 17,9M | 2,2M | 2,2M |
| 2,8301 | ⟨7,7,10⟩ | 345 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 42,1G | 5,4G | 5,4G |
| 2,8301 | ⟨4,8,9⟩ | 209 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 82,5G | 2,0G | 249,1M | 249,1M |
| 2,8302 | ⟨7,9,9⟩ | 396 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 152,5G | 19,5G | 19,5G |
| 2,8302 | ⟨5,6,12⟩ | 258 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 34,7G | 369,0M | 48,4M | 48,4M |
| 2,8302 | ⟨4,5,12⟩ | 176 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,5G | 29,7M | 3,9M | 3,9M |
| 2,8306 | ⟨7,9,11⟩ | 479 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 131,4G | 16,8G | 16,8G |
| 2,8306 | ⟨7,9,12⟩ | 520 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8308 | ⟨4,7,12⟩ | 242 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 19,6G | 396,8M | 50,2M | 50,2M |
| 2,8308 | ⟨7,7,11⟩ | 378 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 33,2G | 4,3G | 4,3G |
| 2,8309 | ⟨5,6,11⟩ | 238 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 636,2M | 84,1M | 84,1M |
| 2,8310 | ⟨8,9,9⟩ | 450 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 380,5G | 48,1G | 48,1G |
| 2,8310 | ⟨6,9,12⟩ | 450 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 245,6G | 19,5G | 2,5G | 2,5G |
| 2,8310 | ⟨4,8,11⟩ | 253 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 1,7G | 214,9M | 214,9M |
| 2,8311 | ⟨6,11,11⟩ | 501 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 29,0G | 3,7G | 3,7G |
| 2,8311 | ⟨5,6,9⟩ | 197 | no | R/Q/Z | Perminov (tensor decomposition) | 102,3G | 738,3M | 97,5M | 97,5M |
| 2,8311 | ⟨3,7,7⟩ | 111 | no | R/Q/Z | meta-flip-graph search | 7,7G | 54,7M | 8,2M | 8,2M |
| 2,8312 | ⟨9,9,11⟩ | 608 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8313 | ⟨3,6,6⟩ | 83 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,3G | 3,6M | 521,7k | 521,7k |
| 2,8313 | ⟨3,4,9⟩ | 83 | no | R/Q/Z | Perminov (tensor decomposition) | 864,9M | 4,5M | 654,7k | 654,7k |
| 2,8313 | ⟨3,6,7⟩ | 96 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 14,0M | 2,1M | 2,1M |
| 2,8314 | ⟨5,7,8⟩ | 204 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 166,7G | 1,8G | 241,7M | 241,7M |
| 2,8314 | ⟨7,10,11⟩ | 530 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 143,7G | 18,3G | 18,3G |
| 2,8315 | ⟨7,7,10⟩ | 346 | no | R/Q/Z | Perminov (tensor decomposition) | 733,6G | 42,1G | 5,4G | 5,4G |
| 2,8315 | ⟨4,11,11⟩ | 342 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 2,3G | 295,7M | 295,7M |
| 2,8315 | ⟨7,11,11⟩ | 580 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 113,2G | 14,5G | 14,5G |
| 2,8315 | ⟨6,6,11⟩ | 283 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 2,2G | 273,2M | 273,2M |
| 2,8316 | ⟨5,5,9⟩ | 166 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 42,6G | 216,3M | 30,0M | 30,0M |
| 2,8316 | ⟨7,9,11⟩ | 480 | no | R/Q/Z | Perminov (tensor decomposition) | 825,3G | 131,4G | 16,8G | 16,8G |
| 2,8316 | ⟨3,7,8⟩ | 126 | no | R/Q/Z | Fmm-lille | 11,0G | 136,5M | 20,1M | 20,1M |
| 2,8317 | ⟨5,10,11⟩ | 386 | no | R/Q/Z | Perminov (tensor decomposition) | 166,7G | 10,8G | 1,4G | 1,4G |
| 2,8319 | ⟨5,12,12⟩ | 498 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 34,7G | 2,9G | 374,7M | 374,7M |
| 2,8320 | ⟨9,9,9⟩ | 504 | yes | R/Q/Z | Solven-strassen-2026 | 2131,7G | 603,6G | 76,8G | -- |
| 2,8320 | ⟨5,8,10⟩ | 286 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 238,2G | 7,9G | 1,0G | 1,0G |
| 2,8321 | ⟨5,5,12⟩ | 218 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8321 | ⟨7,8,10⟩ | 393 | no | R/Q/Z | AlphaTensor 2022 | 1048,0G | 105,1G | 13,3G | 13,3G |
| 2,8321 | ⟨5,8,11⟩ | 313 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 166,7G | 6,2G | 820,4M | 820,4M |
| 2,8323 | ⟨8,9,12⟩ | 592 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 190,2G | 23,9G | 23,9G |
| 2,8323 | ⟨6,6,7⟩ | 185 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 636,5M | 80,5M | 80,5M |
| 2,8323 | ⟨3,7,9⟩ | 141 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 216,5M | 32,1M | 32,1M |
| 2,8324 | ⟨6,8,10⟩ | 340 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 26,9G | 3,4G | 3,4G |
| 2,8324 | ⟨5,8,12⟩ | 340 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 3,6G | 472,6M | 472,6M |
| 2,8324 | ⟨3,5,6⟩ | 70 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 954,0M | 1,0M | 160,5k | 160,5k |
| 2,8325 | ⟨4,6,9⟩ | 160 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 31,5G | 203,0M | 25,5M | 25,5M |
| 2,8325 | ⟨3,6,12⟩ | 160 | no | R/Q/Z | Fmm-lille | 2,3G | 27,7M | 4,0M | 4,0M |
| 2,8325 | ⟨7,9,11⟩ | 481 | no | R/Q/Z | AlphaTensor 2022 | 825,3G | 131,4G | 16,8G | 16,8G |
| 2,8326 | ⟨5,5,11⟩ | 201 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 26,5G | 186,4M | 25,9M | 25,9M |
| 2,8326 | ⟨6,8,12⟩ | 404 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 12,3G | 1,5G | 1,5G |
| 2,8327 | ⟨4,8,9⟩ | 210 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 82,5G | 2,0G | 249,1M | 249,1M |
| 2,8327 | ⟨3,5,7⟩ | 81 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,7G | 4,1M | 635,0k | 635,0k |
| 2,8328 | ⟨3,6,9⟩ | 122 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 55,4M | 8,1M | 8,1M |
| 2,8329 | ⟨6,7,10⟩ | 300 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 10,8G | 1,4G | 1,4G |
| 2,8329 | ⟨5,7,12⟩ | 300 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 1,4G | 191,6M | 191,6M |
| 2,8329 | ⟨4,11,11⟩ | 343 | no | R/Q/Z | AlphaTensor 2022 | 35,9G | 2,3G | 295,7M | 295,7M |
| 2,8330 | ⟨7,11,11⟩ | 582 | no | R/Q/Z | AlphaTensor 2022 | 513,5G | 113,2G | 14,5G | 14,5G |
| 2,8332 | ⟨4,5,11⟩ | 163 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 8,2G | 51,3M | 6,8M | 6,8M |
| 2,8333 | ⟨8,8,11⟩ | 489 | no | R/Q/Z | AlphaTensor 2022 | 1048,0G | 206,7G | 26,0G | 26,0G |
| 2,8333 | ⟨4,7,11⟩ | 224 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 684,1M | 87,1M | 87,1M |
| 2,8333 | ⟨4,5,10⟩ | 149 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8335 | ⟨5,5,10⟩ | 184 | no | R/Q/Z | Perminov (tensor decomposition) | 37,9G | 236,5M | 32,6M | 32,6M |
| 2,8335 | ⟨7,9,11⟩ | 482 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 131,4G | 16,8G | 16,8G |
| 2,8335 | ⟨5,11,11⟩ | 424 | no | R/Q/Z | Fmm-lille | 116,7G | 8,5G | 1,1G | 1,1G |
| 2,8335 | ⟨3,5,8⟩ | 92 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 10,2M | 1,6M | 1,6M |
| 2,8335 | ⟨3,4,10⟩ | 92 | no | R/Q/Z | Perminov (tensor decomposition) | 768,8M | 4,9M | 711,5k | 711,5k |
| 2,8336 | ⟨10,10,11⟩ | 746 | no | R/Q/Z | AlphaTensor 2022 | 1048,0G | 621,7G | 78,2G | 78,2G |
| 2,8337 | ⟨2,7,7⟩ | 76 | no | R/Q/Z | meta-flip-graph search | 1,0G | 18,2M | 2,3M | 2,3M |
| 2,8337 | ⟨5,9,12⟩ | 381 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,3G | 5,7G | 754,2M | 754,2M |
| 2,8337 | ⟨6,8,11⟩ | 373 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 21,2G | 2,7G | 2,7G |
| 2,8337 | ⟨5,8,11⟩ | 314 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 166,7G | 6,2G | 820,4M | 820,4M |
| 2,8337 | ⟨5,7,9⟩ | 229 | no | R/Q/Z | Perminov (tensor decomposition) | 187,6G | 2,9G | 385,7M | 385,7M |
| 2,8337 | ⟨7,9,9⟩ | 399 | no | R/Q/Z | AlphaTensor 2022 | 1326,4G | 152,5G | 19,5G | 19,5G |
| 2,8337 | ⟨7,9,9⟩ | 399 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 152,5G | 19,5G | 19,5G |
| 2,8338 | ⟨5,8,10⟩ | 287 | no | R/Q/Z | AlphaTensor 2022 | 238,2G | 7,9G | 1,0G | 1,0G |
| 2,8338 | ⟨5,8,10⟩ | 287 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 238,2G | 7,9G | 1,0G | 1,0G |
| 2,8338 | ⟨10,10,10⟩ | 682 | yes | R/Q/Z | AlphaTensor 2022 | 1497,2G | 788,9G | 98,6G | -- |
| 2,8338 | ⟨6,11,11⟩ | 504 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 29,0G | 3,7G | 3,7G |
| 2,8338 | ⟨4,5,9⟩ | 135 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 13,1G | 59,5M | 7,9M | 7,9M |
| 2,8338 | ⟨5,10,12⟩ | 421 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 6,2G | 819,6M | 819,6M |
| 2,8340 | ⟨7,9,10⟩ | 441 | no | R/Q/Z | AlphaTensor 2022 | 1179,0G | 166,7G | 21,2G | 21,2G |
| 2,8340 | ⟨5,7,8⟩ | 205 | no | R/Q/Z | Perminov (tensor decomposition) | 166,7G | 1,8G | 241,7M | 241,7M |
| 2,8341 | ⟨5,7,11⟩ | 277 | no | R/Q/Z | Perminov (tensor decomposition) | 116,7G | 2,5G | 332,6M | 332,6M |
| 2,8341 | ⟨5,8,9⟩ | 260 | no | R/Q/Z | Perminov (tensor decomposition) | 268,0G | 7,2G | 951,3M | 951,3M |
| 2,8341 | ⟨5,6,12⟩ | 260 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 34,7G | 369,0M | 48,4M | 48,4M |
| 2,8342 | ⟨5,11,12⟩ | 461 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 4,9G | 650,4M | 650,4M |
| 2,8342 | ⟨7,7,9⟩ | 315 | no | R/Q/Z | Fmm-lille | 825,3G | 38,5G | 5,0G | 5,0G |
| 2,8342 | ⟨4,6,11⟩ | 194 | no | R/Q/Z | Perminov (tensor decomposition) | 19,6G | 174,9M | 22,0M | 22,0M |
| 2,8343 | ⟨8,10,12⟩ | 657 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 208,0G | 26,0G | 26,0G |
| 2,8346 | ⟨7,11,11⟩ | 584 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 113,2G | 14,5G | 14,5G |
| 2,8347 | ⟨10,10,12⟩ | 812 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 571,7G | 360,6G | 45,1G | 45,1G |
| 2,8348 | ⟨6,7,7⟩ | 215 | no | R/Q/Z | Perminov (tensor decomposition) | 280,1G | 2,5G | 318,4M | 318,4M |
| 2,8349 | ⟨4,4,10⟩ | 121 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 3,6G | 17,9M | 2,2M | 2,2M |
| 2,8349 | ⟨5,5,9⟩ | 167 | no | R/Q/Z | Perminov (tensor decomposition) | 42,6G | 216,3M | 30,0M | 30,0M |
| 2,8349 | ⟨3,6,11⟩ | 148 | no | R/Q/Z | Fmm-lille | 4,2G | 47,7M | 7,0M | 7,0M |
| 2,8350 | ⟨11,11,11⟩ | 896 | yes | R/Q/Z | AlphaTensor 2022 | 513,5G | 386,1G | 49,3G | -- |
| 2,8351 | ⟨4,8,11⟩ | 255 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 1,7G | 214,9M | 214,9M |
| 2,8351 | ⟨8,8,11⟩ | 491 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1048,0G | 206,7G | 26,0G | 26,0G |
| 2,8352 | ⟨5,5,11⟩ | 202 | no | R/Q/Z | Perminov (tensor decomposition) | 26,5G | 186,4M | 25,9M | 25,9M |
| 2,8353 | ⟨5,6,11⟩ | 240 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 636,2M | 84,1M | 84,1M |
| 2,8354 | ⟨7,7,7⟩ | 249 | yes | R/Q/Z | Solven-strassen-2026 | 513,5G | 9,7G | 1,3G | -- |
| 2,8355 | ⟨5,8,10⟩ | 288 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 238,2G | 7,9G | 1,0G | 1,0G |
| 2,8355 | ⟨3,4,11⟩ | 101 | no | R/Q/Z | Perminov (tensor decomposition) | 538,1M | 3,8M | 564,6k | 564,6k |
| 2,8356 | ⟨4,7,11⟩ | 225 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 684,1M | 87,1M | 87,1M |
| 2,8356 | ⟨5,9,10⟩ | 322 | no | R/Q/Z | Perminov (tensor decomposition) | 268,0G | 12,5G | 1,6G | 1,6G |
| 2,8357 | ⟨3,5,10⟩ | 114 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 17,7M | 2,7M | 2,7M |
| 2,8358 | ⟨7,7,9⟩ | 316 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 38,5G | 5,0G | 5,0G |
| 2,8358 | ⟨5,7,10⟩ | 254 | no | R/Q/Z | Perminov (tensor decomposition) | 166,7G | 3,2G | 419,1M | 419,1M |
| 2,8358 | ⟨4,4,12⟩ | 144 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,4G | 8,2M | 1,0M | 1,0M |
| 2,8360 | ⟨5,5,5⟩ | 96 | yes | F2 | AlphaTensor 2022 | 6,0G | 4,1M | 592,7k | 102,3k |
| 2,8361 | ⟨2,6,6⟩ | 57 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 305,3M | 1,2M | 149,1k | 149,1k |
| 2,8361 | ⟨3,4,6⟩ | 57 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 293,5M | 288,3k | 42,0k | 42,0k |
| 2,8361 | ⟨3,3,8⟩ | 57 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 164,7M | 767,2k | 130,5k | 130,5k |
| 2,8362 | ⟨2,4,8⟩ | 51 | no | R/Q/Z | Alphaevolve | 102,5M | 937,7k | 117,2k | 117,2k |
| 2,8363 | ⟨10,10,10⟩ | 686 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 1497,2G | 788,9G | 98,6G | -- |
| 2,8364 | ⟨4,5,12⟩ | 178 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,5G | 29,7M | 3,9M | 3,9M |
| 2,8365 | ⟨4,4,5⟩ | 63 | no | F2 | AlphaTensor 2022 | 570,8M | 309,8k | 40,7k | 40,7k |
| 2,8365 | ⟨4,4,5⟩ | 63 | no | R/Q/Z | AlphaTensor 2022 | 570,8M | 309,8k | 40,7k | 40,7k |
| 2,8365 | ⟨2,5,8⟩ | 63 | no | R/Q/Z | meta-flip-graph search | 333,1M | 3,4M | 447,6k | 447,6k |
| 2,8365 | ⟨7,10,11⟩ | 536 | no | R/Q/Z | AlphaTensor 2022 | 733,6G | 143,7G | 18,3G | 18,3G |
| 2,8365 | ⟨5,9,11⟩ | 353 | no | R/Q/Z | Perminov (tensor decomposition) | 187,6G | 9,8G | 1,3G | 1,3G |
| 2,8366 | ⟨5,10,11⟩ | 390 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 166,7G | 10,8G | 1,4G | 1,4G |
| 2,8366 | ⟨5,7,8⟩ | 206 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 166,7G | 1,8G | 241,7M | 241,7M |
| 2,8366 | ⟨4,7,10⟩ | 206 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 868,0M | 109,8M | 109,8M |
| 2,8367 | ⟨2,6,7⟩ | 66 | no | R/Q/Z | meta-flip-graph search | 559,7M | 4,7M | 589,7k | 589,7k |
| 2,8368 | ⟨5,11,11⟩ | 427 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 116,7G | 8,5G | 1,1G | 1,1G |
| 2,8368 | ⟨3,7,10⟩ | 157 | no | R/Q/Z | Fmm-lille | 11,0G | 236,7M | 34,9M | 34,9M |
| 2,8369 | ⟨5,5,12⟩ | 220 | no | R/Q/Z | Perminov (tensor decomposition) | 14,5G | 108,1M | 14,9M | 14,9M |
| 2,8370 | ⟨10,10,10⟩ | 687 | yes | R/Q/Z | Solven-strassen-2026 | 1497,2G | 788,9G | 98,6G | -- |
| 2,8370 | ⟨2,4,7⟩ | 45 | no | R/Q/Z | Alphaevolve | 71,8M | 375,8k | 47,5k | 47,5k |
| 2,8370 | ⟨2,4,7⟩ | 45 | no | R/Q/Z | AlphaEvolve 2025 | 71,8M | 375,8k | 47,5k | 47,5k |
| 2,8370 | ⟨7,7,10⟩ | 350 | no | R/Q/Z | AlphaTensor 2022 | 733,6G | 42,1G | 5,4G | 5,4G |
| 2,8371 | ⟨4,8,11⟩ | 256 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 1,7G | 214,9M | 214,9M |
| 2,8371 | ⟨4,5,10⟩ | 150 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8371 | ⟨8,10,11⟩ | 609 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1048,0G | 358,5G | 45,1G | 45,1G |
| 2,8372 | ⟨5,9,10⟩ | 323 | no | R/Q/Z | AlphaTensor 2022 | 268,0G | 12,5G | 1,6G | 1,6G |
| 2,8372 | ⟨5,9,10⟩ | 323 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 268,0G | 12,5G | 1,6G | 1,6G |
| 2,8372 | ⟨6,9,11⟩ | 420 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 33,6G | 4,3G | 4,3G |
| 2,8374 | ⟨6,10,12⟩ | 504 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 218,3G | 21,3G | 2,7G | 2,7G |
| 2,8374 | ⟨4,4,9⟩ | 110 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,0G | 16,4M | 2,1M | 2,1M |
| 2,8374 | ⟨6,9,10⟩ | 384 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 42,6G | 5,4G | 5,4G |
| 2,8375 | ⟨7,7,7⟩ | 250 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 9,7G | 1,3G | -- |
| 2,8377 | ⟨3,6,7⟩ | 97 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 14,0M | 2,1M | 2,1M |
| 2,8377 | ⟨2,6,8⟩ | 75 | no | R/Q/Z | meta-flip-graph search | 799,5M | 11,6M | 1,5M | 1,5M |
| 2,8379 | ⟨5,9,11⟩ | 354 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 187,6G | 9,8G | 1,3G | 1,3G |
| 2,8379 | ⟨4,7,11⟩ | 226 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 684,1M | 87,1M | 87,1M |
| 2,8380 | ⟨5,8,9⟩ | 262 | no | R/Q/Z | AlphaTensor 2022 | 268,0G | 7,2G | 951,3M | 951,3M |
| 2,8381 | ⟨3,6,10⟩ | 136 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 60,5M | 8,8M | 8,8M |
| 2,8381 | ⟨3,5,12⟩ | 136 | no | R/Q/Z | Perminov (tensor decomposition) | 954,0M | 8,1M | 1,2M | 1,2M |
| 2,8381 | ⟨4,7,9⟩ | 187 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 57,7G | 793,9M | 101,0M | 101,0M |
| 2,8382 | ⟨5,11,12⟩ | 465 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 63,7G | 4,9G | 650,4M | 650,4M |
| 2,8382 | ⟨5,9,9⟩ | 293 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 301,5G | 11,4G | 1,5G | 1,5G |
| 2,8383 | ⟨5,7,9⟩ | 231 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 187,6G | 2,9G | 385,7M | 385,7M |
| 2,8383 | ⟨7,7,11⟩ | 384 | no | R/Q/Z | AlphaTensor 2022 | 513,5G | 33,2G | 4,3G | 4,3G |
| 2,8384 | ⟨5,8,11⟩ | 317 | no | R/Q/Z | AlphaTensor 2022 | 166,7G | 6,2G | 820,4M | 820,4M |
| 2,8388 | ⟨9,10,12⟩ | 742 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 643,1G | 329,8G | 41,5G | 41,5G |
| 2,8389 | ⟨7,7,9⟩ | 318 | no | R/Q/Z | AlphaTensor 2022 | 825,3G | 38,5G | 5,0G | 5,0G |
| 2,8390 | ⟨3,4,9⟩ | 84 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 864,9M | 4,5M | 654,7k | 654,7k |
| 2,8390 | ⟨3,3,12⟩ | 84 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 62,9M | 608,3k | 103,5k | 103,5k |
| 2,8391 | ⟨2,4,6⟩ | 39 | no | R/Q/Z | Perminov (tensor decomposition) | 39,1M | 96,1k | 12,0k | 12,0k |
| 2,8392 | ⟨7,8,12⟩ | 474 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 48,0G | 6,1G | 6,1G |
| 2,8392 | ⟨4,7,10⟩ | 207 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 51,3G | 868,0M | 109,8M | 109,8M |
| 2,8392 | ⟨3,3,9⟩ | 64 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 185,3M | 1,2M | 208,3k | 208,3k |
| 2,8392 | ⟨2,3,5⟩ | 25 | no | F2 | Hopcroft-Kerr 1971 | 3,5M | 7,7k | 1,2k | 1,2k |
| 2,8392 | ⟨2,3,5⟩ | 25 | no | R/Q/Z | AlphaTensor 2022 | 3,5M | 7,7k | 1,2k | 1,2k |
| 2,8392 | ⟨4,7,12⟩ | 246 | no | R/Q/Z | Perminov (tensor decomposition) | 19,6G | 396,8M | 50,2M | 50,2M |
| 2,8393 | ⟨5,9,11⟩ | 355 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 187,6G | 9,8G | 1,3G | 1,3G |
| 2,8394 | ⟨6,9,9⟩ | 349 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 39,0G | 4,9G | 4,9G |
| 2,8394 | ⟨3,6,12⟩ | 162 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,3G | 27,7M | 4,0M | 4,0M |
| 2,8395 | ⟨4,5,12⟩ | 179 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,5G | 29,7M | 3,9M | 3,9M |
| 2,8395 | ⟨5,7,11⟩ | 280 | no | R/Q/Z | AlphaTensor 2022 | 116,7G | 2,5G | 332,6M | 332,6M |
| 2,8397 | ⟨4,4,10⟩ | 122 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 3,6G | 17,9M | 2,2M | 2,2M |
| 2,8398 | ⟨3,8,8⟩ | 145 | no | R/Q/Z | Fmm-lille | 15,7G | 340,7M | 49,7M | 49,7M |
| 2,8398 | ⟨4,4,12⟩ | 145 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,4G | 8,2M | 1,0M | 1,0M |
| 2,8398 | ⟨7,7,10⟩ | 352 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 42,1G | 5,4G | 5,4G |
| 2,8399 | ⟨6,7,11⟩ | 333 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 8,5G | 1,1G | 1,1G |
| 2,8400 | ⟨4,5,11⟩ | 165 | no | R/Q/Z | Perminov (tensor decomposition) | 8,2G | 51,3M | 6,8M | 6,8M |
| 2,8402 | ⟨4,7,11⟩ | 227 | no | R/Q/Z | Perminov (tensor decomposition) | 35,9G | 684,1M | 87,1M | 87,1M |
| 2,8403 | ⟨3,4,10⟩ | 93 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 768,8M | 4,9M | 711,5k | 711,5k |
| 2,8404 | ⟨3,7,10⟩ | 158 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 236,7M | 34,9M | 34,9M |
| 2,8404 | ⟨3,7,9⟩ | 143 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 216,5M | 32,1M | 32,1M |
| 2,8404 | ⟨3,5,9⟩ | 104 | no | R/Q/Z | Perminov (tensor decomposition) | 2,8G | 16,2M | 2,5M | 2,5M |
| 2,8405 | ⟨6,7,12⟩ | 362 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 152,8G | 4,9G | 622,7M | 622,7M |
| 2,8405 | ⟨10,10,11⟩ | 758 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1048,0G | 621,7G | 78,2G | 78,2G |
| 2,8406 | ⟨3,5,7⟩ | 82 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,7G | 4,1M | 635,0k | 635,0k |
| 2,8406 | ⟨3,7,11⟩ | 173 | no | R/Q/Z | Fmm-lille | 7,7G | 186,6M | 27,7M | 27,7M |
| 2,8407 | ⟨7,8,11⟩ | 438 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 82,8G | 10,5G | 10,5G |
| 2,8408 | ⟨3,7,8⟩ | 128 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 136,5M | 20,1M | 20,1M |
| 2,8409 | ⟨4,5,10⟩ | 151 | no | R/Q/Z | Perminov (tensor decomposition) | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8409 | ⟨3,5,10⟩ | 115 | no | R/Q/Z | Perminov (tensor decomposition) | 2,5G | 17,7M | 2,7M | 2,7M |
| 2,8410 | ⟨9,9,12⟩ | 675 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 723,5G | 301,7G | 38,2G | 38,2G |
| 2,8410 | ⟨12,12,12⟩ | 1164 | yes | R/Q/Z | Solven-strassen-2026 | 83,3G | 75,4G | 9,4G | -- |
| 2,8410 | ⟨4,7,9⟩ | 188 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 57,7G | 793,9M | 101,0M | 101,0M |
| 2,8410 | ⟨3,7,12⟩ | 188 | no | R/Q/Z | Fmm-lille | 4,2G | 108,2M | 16,0M | 16,0M |
| 2,8416 | ⟨3,5,11⟩ | 126 | no | R/Q/Z | Perminov (tensor decomposition) | 1,7G | 14,0M | 2,2M | 2,2M |
| 2,8416 | ⟨3,4,11⟩ | 102 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 538,1M | 3,8M | 564,6k | 564,6k |
| 2,8416 | ⟨5,9,9⟩ | 295 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 301,5G | 11,4G | 1,5G | 1,5G |
| 2,8417 | ⟨5,9,10⟩ | 326 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 268,0G | 12,5G | 1,6G | 1,6G |
| 2,8418 | ⟨4,4,11⟩ | 134 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 14,1M | 1,8M | 1,8M |
| 2,8418 | ⟨5,7,10⟩ | 257 | no | R/Q/Z | AlphaTensor 2022 | 166,7G | 3,2G | 419,1M | 419,1M |
| 2,8418 | ⟨8,8,9⟩ | 412 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 239,9G | 30,2G | 30,2G |
| 2,8419 | ⟨3,7,7⟩ | 113 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 54,7M | 8,2M | 8,2M |
| 2,8419 | ⟨3,3,10⟩ | 71 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 164,7M | 1,3M | 226,4k | 226,4k |
| 2,8421 | ⟨6,10,11⟩ | 469 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 36,7G | 4,6G | 4,6G |
| 2,8422 | ⟨2,7,7⟩ | 77 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,0G | 18,2M | 2,3M | 2,3M |
| 2,8423 | ⟨5,11,11⟩ | 432 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 116,7G | 8,5G | 1,1G | 1,1G |
| 2,8423 | ⟨4,5,9⟩ | 137 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 13,1G | 59,5M | 7,9M | 7,9M |
| 2,8423 | ⟨3,6,10⟩ | 137 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 60,5M | 8,8M | 8,8M |
| 2,8424 | ⟨3,6,9⟩ | 124 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 55,4M | 8,1M | 8,1M |
| 2,8424 | ⟨2,5,7⟩ | 56 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 233,2M | 1,4M | 181,4k | 181,4k |
| 2,8425 | ⟨3,6,11⟩ | 150 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 47,7M | 7,0M | 7,0M |
| 2,8425 | ⟨3,8,10⟩ | 180 | no | R/Q/Z | Perminov (tensor decomposition) | 15,7G | 590,8M | 86,2M | 86,2M |
| 2,8425 | ⟨4,5,12⟩ | 180 | no | R/Q/Z | Perminov (tensor decomposition) | 4,5G | 29,7M | 3,9M | 3,9M |
| 2,8426 | ⟨8,9,11⟩ | 558 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 327,9G | 41,5G | 41,5G |
| 2,8427 | ⟨4,8,9⟩ | 214 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 82,5G | 2,0G | 249,1M | 249,1M |
| 2,8429 | ⟨3,8,9⟩ | 163 | no | R/Q/Z | Perminov (tensor decomposition) | 17,7G | 540,3M | 79,3M | 79,3M |
| 2,8429 | ⟨3,4,12⟩ | 111 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 293,5M | 2,2M | 325,2k | 325,2k |
| 2,8432 | ⟨7,7,11⟩ | 388 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 33,2G | 4,3G | 4,3G |
| 2,8433 | ⟨5,9,11⟩ | 358 | no | R/Q/Z | AlphaTensor 2022 | 187,6G | 9,8G | 1,3G | 1,3G |
| 2,8433 | ⟨5,9,9⟩ | 296 | no | R/Q/Z | AlphaTensor 2022 | 301,5G | 11,4G | 1,5G | 1,5G |
| 2,8436 | ⟨2,4,5⟩ | 33 | no | F2 | Hopcroft-Kerr 1971 | 16,3M | 28,2k | 3,7k | 3,7k |
| 2,8436 | ⟨2,4,5⟩ | 33 | no | R/Q/Z | AlphaTensor 2022 | 16,3M | 28,2k | 3,7k | 3,7k |
| 2,8437 | ⟨3,8,8⟩ | 146 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 340,7M | 49,7M | 49,7M |
| 2,8440 | ⟨7,8,10⟩ | 403 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1048,0G | 105,1G | 13,3G | 13,3G |
| 2,8441 | ⟨9,10,10⟩ | 632 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1684,3G | 721,5G | 90,7G | 90,7G |
| 2,8443 | ⟨3,9,12⟩ | 240 | no | R/Q/Z | Fmm-lille | 6,7G | 428,4M | 62,8M | 62,8M |
| 2,8443 | ⟨7,9,12⟩ | 536 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 450,2G | 76,2G | 9,7G | 9,7G |
| 2,8443 | ⟨3,3,11⟩ | 78 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 115,3M | 1,0M | 179,6k | 179,6k |
| 2,8446 | ⟨4,5,10⟩ | 152 | no | R/Q/Z | AlphaTensor 2022 | 11,7G | 65,0M | 8,5M | 8,5M |
| 2,8449 | ⟨5,5,7⟩ | 134 | no | R/Q/Z | AlphaTensor 2022 | 26,5G | 54,7M | 7,6M | 7,6M |
| 2,8450 | ⟨5,7,9⟩ | 234 | no | R/Q/Z | AlphaTensor 2022 | 187,6G | 2,9G | 385,7M | 385,7M |
| 2,8451 | ⟨7,7,12⟩ | 423 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 19,3G | 2,5G | 2,5G |
| 2,8451 | ⟨3,9,9⟩ | 183 | no | R/Q/Z | Fmm-lille | 19,9G | 857,0M | 126,5M | 126,5M |
| 2,8452 | ⟨3,8,11⟩ | 198 | no | R/Q/Z | Perminov (tensor decomposition) | 11,0G | 465,6M | 68,4M | 68,4M |
| 2,8455 | ⟨9,10,11⟩ | 694 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 568,6G | 72,0G | 72,0G |
| 2,8456 | ⟨3,8,10⟩ | 181 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 590,8M | 86,2M | 86,2M |
| 2,8457 | ⟨9,9,10⟩ | 574 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1894,9G | 659,9G | 83,5G | 83,5G |
| 2,8460 | ⟨5,9,11⟩ | 360 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 187,6G | 9,8G | 1,3G | 1,3G |
| 2,8461 | ⟨3,5,10⟩ | 116 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 17,7M | 2,7M | 2,7M |
| 2,8463 | ⟨3,8,9⟩ | 164 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 17,7G | 540,3M | 79,3M | 79,3M |
| 2,8463 | ⟨3,5,9⟩ | 105 | no | R/Q/Z | AlphaTensor 2022 | 2,8G | 16,2M | 2,5M | 2,5M |
| 2,8464 | ⟨2,5,10⟩ | 79 | no | R/Q/Z | Perminov (tensor decomposition) | 333,1M | 5,9M | 776,2k | 776,2k |
| 2,8465 | ⟨2,6,8⟩ | 76 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 799,5M | 11,6M | 1,5M | 1,5M |
| 2,8466 | ⟨3,6,6⟩ | 85 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,3G | 3,6M | 521,7k | 521,7k |
| 2,8466 | ⟨3,4,9⟩ | 85 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 864,9M | 4,5M | 654,7k | 654,7k |
| 2,8466 | ⟨3,9,11⟩ | 222 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 738,5M | 109,1M | 109,1M |
| 2,8467 | ⟨2,4,11⟩ | 70 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 71,8M | 1,3M | 161,3k | 161,3k |
| 2,8467 | ⟨2,7,8⟩ | 88 | no | R/Q/Z | Perminov (tensor decomposition) | 1,5G | 45,5M | 5,8M | 5,8M |
| 2,8468 | ⟨3,7,12⟩ | 190 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 108,2M | 16,0M | 16,0M |
| 2,8469 | ⟨2,6,7⟩ | 67 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 4,7M | 589,7k | 589,7k |
| 2,8470 | ⟨3,7,11⟩ | 175 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 186,6M | 27,7M | 27,7M |
| 2,8470 | ⟨3,5,8⟩ | 94 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,5G | 10,2M | 1,6M | 1,6M |
| 2,8470 | ⟨3,4,10⟩ | 94 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 768,8M | 4,9M | 711,5k | 711,5k |
| 2,8470 | ⟨2,6,10⟩ | 94 | no | R/Q/Z | Perminov (tensor decomposition) | 799,5M | 20,2M | 2,5M | 2,5M |
| 2,8470 | ⟨2,5,12⟩ | 94 | no | R/Q/Z | Perminov (tensor decomposition) | 127,2M | 2,7M | 354,8k | 354,8k |
| 2,8470 | ⟨3,3,7⟩ | 51 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 115,3M | 307,4k | 52,9k | 52,9k |
| 2,8471 | ⟨7,8,9⟩ | 367 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 96,1G | 12,2G | 12,2G |
| 2,8472 | ⟨3,9,10⟩ | 203 | no | R/Q/Z | Fmm-lille | 17,7G | 937,0M | 137,5M | 137,5M |
| 2,8472 | ⟨2,5,8⟩ | 64 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 333,1M | 3,4M | 447,6k | 447,6k |
| 2,8472 | ⟨2,4,10⟩ | 64 | no | R/Q/Z | Perminov (tensor decomposition) | 102,5M | 1,6M | 203,3k | 203,3k |
| 2,8474 | ⟨2,3,6⟩ | 30 | no | R/Q/Z | Perminov (tensor decomposition) | 8,4M | 26,2k | 3,8k | 3,8k |
| 2,8474 | ⟨2,8,8⟩ | 100 | no | R/Q/Z | Perminov (tensor decomposition) | 2,1G | 113,6M | 14,2M | 14,2M |
| 2,8474 | ⟨3,7,10⟩ | 160 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 236,7M | 34,9M | 34,9M |
| 2,8476 | ⟨3,4,11⟩ | 103 | no | R/Q/Z | AlphaTensor 2022 | 538,1M | 3,8M | 564,6k | 564,6k |
| 2,8476 | ⟨2,6,11⟩ | 103 | no | R/Q/Z | Perminov (tensor decomposition) | 559,7M | 15,9M | 2,0M | 2,0M |
| 2,8476 | ⟨3,8,12⟩ | 216 | no | R/Q/Z | Perminov (tensor decomposition) | 6,0G | 270,1M | 39,4M | 39,4M |
| 2,8476 | ⟨5,11,11⟩ | 437 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 116,7G | 8,5G | 1,1G | 1,1G |
| 2,8481 | ⟨5,10,12⟩ | 434 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 90,9G | 6,2G | 819,6M | 819,6M |
| 2,8482 | ⟨4,11,11⟩ | 354 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 2,3G | 295,7M | 295,7M |
| 2,8483 | ⟨3,6,8⟩ | 112 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 34,9M | 5,1M | 5,1M |
| 2,8483 | ⟨2,6,12⟩ | 112 | no | R/Q/Z | Perminov (tensor decomposition) | 305,3M | 9,2M | 1,2M | 1,2M |
| 2,8483 | ⟨3,4,12⟩ | 112 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 293,5M | 2,2M | 325,2k | 325,2k |
| 2,8483 | ⟨3,3,8⟩ | 58 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 164,7M | 767,2k | 130,5k | 130,5k |
| 2,8483 | ⟨2,4,9⟩ | 58 | no | R/Q/Z | Perminov (tensor decomposition) | 115,3M | 1,5M | 187,0k | 187,0k |
| 2,8483 | ⟨3,7,9⟩ | 145 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 216,5M | 32,1M | 32,1M |
| 2,8484 | ⟨3,5,7⟩ | 83 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,7G | 4,1M | 635,0k | 635,0k |
| 2,8485 | ⟨7,9,10⟩ | 455 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1179,0G | 166,7G | 21,2G | 21,2G |
| 2,8486 | ⟨3,8,10⟩ | 182 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 590,8M | 86,2M | 86,2M |
| 2,8486 | ⟨5,10,11⟩ | 400 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 166,7G | 10,8G | 1,4G | 1,4G |
| 2,8488 | ⟨5,5,5⟩ | 98 | yes | R/Q/Z | Sedoglavic-Smirnov 2021 | 6,0G | 4,1M | 592,7k | 102,3k |
| 2,8489 | ⟨7,9,9⟩ | 412 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 152,5G | 19,5G | 19,5G |
| 2,8493 | ⟨7,10,12⟩ | 599 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 400,2G | 83,3G | 10,5G | 10,5G |
| 2,8494 | ⟨7,9,11⟩ | 499 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 131,4G | 16,8G | 16,8G |
| 2,8496 | ⟨3,10,12⟩ | 268 | no | R/Q/Z | Fmm-lille | 6,0G | 468,4M | 68,3M | 68,3M |
| 2,8497 | ⟨5,9,12⟩ | 394 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,3G | 5,7G | 754,2M | 754,2M |
| 2,8498 | ⟨3,3,9⟩ | 65 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 185,3M | 1,2M | 208,3k | 208,3k |
| 2,8498 | ⟨3,9,10⟩ | 204 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 17,7G | 937,0M | 137,5M | 137,5M |
| 2,8503 | ⟨2,5,11⟩ | 87 | no | R/Q/Z | Perminov (tensor decomposition) | 233,2M | 4,7M | 615,9k | 615,9k |
| 2,8504 | ⟨8,11,11⟩ | 687 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 282,5G | 35,8G | 35,8G |
| 2,8504 | ⟨3,6,7⟩ | 99 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 14,0M | 2,1M | 2,1M |
| 2,8504 | ⟨2,7,9⟩ | 99 | no | R/Q/Z | Perminov (tensor decomposition) | 1,6G | 72,2M | 9,2M | 9,2M |
| 2,8506 | ⟨3,8,11⟩ | 200 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 465,6M | 68,4M | 68,4M |
| 2,8506 | ⟨6,11,11⟩ | 523 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 280,1G | 29,0G | 3,7G | 3,7G |
| 2,8507 | ⟨4,5,9⟩ | 139 | no | R/Q/Z | AlphaTensor 2022 | 13,1G | 59,5M | 7,9M | 7,9M |
| 2,8507 | ⟨3,6,10⟩ | 139 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 60,5M | 8,8M | 8,8M |
| 2,8507 | ⟨3,5,12⟩ | 139 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 954,0M | 8,1M | 1,2M | 1,2M |
| 2,8507 | ⟨3,9,12⟩ | 243 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 428,4M | 62,8M | 62,8M |
| 2,8508 | ⟨3,5,11⟩ | 128 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,7G | 14,0M | 2,2M | 2,2M |
| 2,8510 | ⟨3,10,10⟩ | 226 | no | R/Q/Z | Fmm-lille | 15,7G | 1,0G | 149,4M | 149,4M |
| 2,8511 | ⟨3,9,9⟩ | 185 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 19,9G | 857,0M | 126,5M | 126,5M |
| 2,8512 | ⟨2,5,9⟩ | 72 | no | R/Q/Z | Perminov (tensor decomposition) | 374,8M | 5,4M | 714,2k | 714,2k |
| 2,8512 | ⟨3,3,10⟩ | 72 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 164,7M | 1,3M | 226,4k | 226,4k |
| 2,8514 | ⟨3,9,11⟩ | 224 | no | R/Q/Z | Perminov (tensor decomposition) | 12,4G | 738,5M | 109,1M | 109,1M |
| 2,8514 | ⟨7,10,11⟩ | 554 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 143,7G | 18,3G | 18,3G |
| 2,8515 | ⟨3,8,8⟩ | 148 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 340,7M | 49,7M | 49,7M |
| 2,8516 | ⟨3,12,12⟩ | 320 | no | R/Q/Z | Fmm-lille | 2,3G | 214,1M | 31,2M | 31,2M |
| 2,8518 | ⟨9,9,11⟩ | 637 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1326,4G | 520,1G | 66,2G | 66,2G |
| 2,8518 | ⟨2,9,9⟩ | 126 | no | R/Q/Z | Fmm-lille | 2,7G | 285,7M | 36,1M | 36,1M |
| 2,8518 | ⟨3,6,9⟩ | 126 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 55,4M | 8,1M | 8,1M |
| 2,8519 | ⟨12,12,12⟩ | 1196 | yes | R/Q/Z | Dis09-Q | 83,3G | 75,4G | 9,4G | -- |
| 2,8522 | ⟨3,10,11⟩ | 248 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 807,4M | 118,6M | 118,6M |
| 2,8524 | ⟨3,9,10⟩ | 205 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 17,7G | 937,0M | 137,5M | 137,5M |
| 2,8524 | ⟨3,7,7⟩ | 115 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 54,7M | 8,2M | 8,2M |
| 2,8525 | ⟨3,7,12⟩ | 192 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 108,2M | 16,0M | 16,0M |
| 2,8527 | ⟨3,3,11⟩ | 79 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 115,3M | 1,0M | 179,6k | 179,6k |
| 2,8528 | ⟨3,9,12⟩ | 244 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 428,4M | 62,8M | 62,8M |
| 2,8532 | ⟨11,11,11⟩ | 936 | yes | R/Q/Z | Solven-strassen-2026 | 513,5G | 386,1G | 49,3G | -- |
| 2,8532 | ⟨3,7,11⟩ | 177 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 186,6M | 27,7M | 27,7M |
| 2,8533 | ⟨4,5,11⟩ | 169 | no | R/Q/Z | AlphaTensor 2022 | 8,2G | 51,3M | 6,8M | 6,8M |
| 2,8533 | ⟨3,10,10⟩ | 227 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 1,0G | 149,4M | 149,4M |
| 2,8534 | ⟨3,10,12⟩ | 270 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 468,4M | 68,3M | 68,3M |
| 2,8535 | ⟨2,6,11⟩ | 104 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 15,9M | 2,0M | 2,0M |
| 2,8536 | ⟨2,7,10⟩ | 110 | no | R/Q/Z | Perminov (tensor decomposition) | 1,5G | 78,9M | 10,0M | 10,0M |
| 2,8537 | ⟨3,6,8⟩ | 113 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 34,9M | 5,1M | 5,1M |
| 2,8537 | ⟨2,8,9⟩ | 113 | no | R/Q/Z | Perminov (tensor decomposition) | 2,4G | 180,1M | 22,6M | 22,6M |
| 2,8537 | ⟨2,6,12⟩ | 113 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 305,3M | 9,2M | 1,2M | 1,2M |
| 2,8537 | ⟨2,3,7⟩ | 35 | no | R/Q/Z | Perminov (tensor decomposition) | 15,4M | 102,5k | 15,1k | 15,1k |
| 2,8537 | ⟨3,9,11⟩ | 225 | no | R/Q/Z | AlphaTensor 2022 | 12,4G | 738,5M | 109,1M | 109,1M |
| 2,8539 | ⟨7,11,11⟩ | 610 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 113,2G | 14,5G | 14,5G |
| 2,8540 | ⟨3,11,12⟩ | 296 | no | R/Q/Z | Fmm-lille | 4,2G | 369,1M | 54,2M | 54,2M |
| 2,8540 | ⟨3,3,3⟩ | 23 | yes | R/Q/Z | Alphaevolve | 1,7M | 1,7k | 343 | 84 |
| 2,8540 | ⟨3,3,3⟩ | 23 | yes | F2 | Laderman 1976 | 1,7M | 1,7k | 343 | 84 |
| 2,8541 | ⟨2,6,9⟩ | 86 | no | R/Q/Z | Perminov (tensor decomposition) | 899,5M | 18,5M | 2,3M | 2,3M |
| 2,8541 | ⟨3,3,12⟩ | 86 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 62,9M | 608,3k | 103,5k | 103,5k |
| 2,8541 | ⟨2,8,10⟩ | 125 | no | R/Q/Z | Perminov (tensor decomposition) | 2,1G | 196,9M | 24,6M | 24,6M |
| 2,8543 | ⟨3,10,11⟩ | 249 | no | R/Q/Z | Perminov (tensor decomposition) | 11,0G | 807,4M | 118,6M | 118,6M |
| 2,8544 | ⟨2,7,12⟩ | 131 | no | R/Q/Z | Perminov (tensor decomposition) | 559,7M | 36,1M | 4,6M | 4,6M |
| 2,8546 | ⟨2,5,10⟩ | 80 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 333,1M | 5,9M | 776,2k | 776,2k |
| 2,8548 | ⟨2,9,10⟩ | 140 | no | R/Q/Z | Fmm-lille | 2,4G | 312,3M | 39,3M | 39,3M |
| 2,8548 | ⟨3,6,10⟩ | 140 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,0G | 60,5M | 8,8M | 8,8M |
| 2,8548 | ⟨3,5,12⟩ | 140 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 954,0M | 8,1M | 1,2M | 1,2M |
| 2,8549 | ⟨2,5,7⟩ | 57 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 233,2M | 1,4M | 181,4k | 181,4k |
| 2,8550 | ⟨2,6,8⟩ | 77 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 799,5M | 11,6M | 1,5M | 1,5M |
| 2,8550 | ⟨2,4,12⟩ | 77 | no | R/Q/Z | Perminov (tensor decomposition) | 39,1M | 743,4k | 92,9k | 92,9k |
| 2,8555 | ⟨4,7,7⟩ | 152 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 35,9G | 200,6M | 25,7M | 25,7M |
| 2,8555 | ⟨2,2,4⟩ | 14 | no | F2 | Hopcroft-Kerr 1971 | 143,4k | 704 | 88 | 88 |
| 2,8555 | ⟨2,2,4⟩ | 14 | no | R/Q/Z | AlphaTensor 2022 | 143,4k | 704 | 88 | 88 |
| 2,8557 | ⟨3,10,10⟩ | 228 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 1,0G | 149,4M | 149,4M |
| 2,8557 | ⟨2,10,10⟩ | 155 | no | R/Q/Z | Fmm-lille | 2,1G | 341,5M | 42,7M | 42,7M |
| 2,8557 | ⟨10,11,11⟩ | 860 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 733,6G | 490,0G | 62,1G | 62,1G |
| 2,8561 | ⟨3,9,11⟩ | 226 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 738,5M | 109,1M | 109,1M |
| 2,8562 | ⟨3,7,9⟩ | 147 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 216,5M | 32,1M | 32,1M |
| 2,8562 | ⟨2,4,11⟩ | 71 | no | R/Q/Z | Perminov (tensor decomposition) | 71,8M | 1,3M | 161,3k | 161,3k |
| 2,8563 | ⟨9,11,11⟩ | 779 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 825,3G | 448,1G | 57,1G | 57,1G |
| 2,8564 | ⟨2,7,11⟩ | 121 | no | R/Q/Z | Perminov (tensor decomposition) | 1,0G | 62,2M | 7,9M | 7,9M |
| 2,8564 | ⟨3,10,11⟩ | 250 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 807,4M | 118,6M | 118,6M |
| 2,8565 | ⟨3,6,9⟩ | 127 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 55,4M | 8,1M | 8,1M |
| 2,8566 | ⟨3,6,7⟩ | 100 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 14,0M | 2,1M | 2,1M |
| 2,8566 | ⟨2,7,9⟩ | 100 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,6G | 72,2M | 9,2M | 9,2M |
| 2,8568 | ⟨3,11,11⟩ | 274 | no | R/Q/Z | Perminov (tensor decomposition) | 7,7G | 636,3M | 94,1M | 94,1M |
| 2,8569 | ⟨2,6,7⟩ | 68 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 4,7M | 589,7k | 589,7k |
| 2,8569 | ⟨3,9,9⟩ | 187 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 19,9G | 857,0M | 126,5M | 126,5M |
| 2,8571 | ⟨3,9,12⟩ | 246 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 6,7G | 428,4M | 62,8M | 62,8M |
| 2,8574 | ⟨3,11,12⟩ | 298 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 369,1M | 54,2M | 54,2M |
| 2,8574 | ⟨2,9,11⟩ | 154 | no | R/Q/Z | Fmm-lille | 1,6G | 246,2M | 31,2M | 31,2M |
| 2,8576 | ⟨3,9,10⟩ | 207 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 17,7G | 937,0M | 137,5M | 137,5M |
| 2,8578 | ⟨3,12,12⟩ | 324 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,3G | 214,1M | 31,2M | 31,2M |
| 2,8578 | ⟨2,5,8⟩ | 65 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 333,1M | 3,4M | 447,6k | 447,6k |
| 2,8578 | ⟨2,4,10⟩ | 65 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 102,5M | 1,6M | 203,3k | 203,3k |
| 2,8580 | ⟨3,10,10⟩ | 229 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 15,7G | 1,0G | 149,4M | 149,4M |
| 2,8584 | ⟨3,9,11⟩ | 227 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 738,5M | 109,1M | 109,1M |
| 2,8584 | ⟨3,10,11⟩ | 251 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 807,4M | 118,6M | 118,6M |
| 2,8587 | ⟨2,3,8⟩ | 40 | no | R/Q/Z | Perminov (tensor decomposition) | 22,0M | 255,7k | 37,3k | 37,3k |
| 2,8588 | ⟨4,5,8⟩ | 126 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,7G | 37,5M | 4,9M | 4,9M |
| 2,8588 | ⟨2,8,10⟩ | 126 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,1G | 196,9M | 24,6M | 24,6M |
| 2,8588 | ⟨3,7,8⟩ | 132 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 11,0G | 136,5M | 20,1M | 20,1M |
| 2,8588 | ⟨2,7,12⟩ | 132 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 36,1M | 4,6M | 4,6M |
| 2,8589 | ⟨2,8,11⟩ | 138 | no | R/Q/Z | Perminov (tensor decomposition) | 1,5G | 155,2M | 19,5M | 19,5M |
| 2,8590 | ⟨2,8,9⟩ | 114 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,4G | 180,1M | 22,6M | 22,6M |
| 2,8590 | ⟨2,6,12⟩ | 114 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 305,3M | 9,2M | 1,2M | 1,2M |
| 2,8591 | ⟨2,11,11⟩ | 187 | no | R/Q/Z | Fmm-lille | 1,0G | 212,1M | 26,9M | 26,9M |
| 2,8591 | ⟨2,8,12⟩ | 150 | no | R/Q/Z | Perminov (tensor decomposition) | 799,5M | 90,0M | 11,3M | 11,3M |
| 2,8597 | ⟨2,9,12⟩ | 168 | no | R/Q/Z | Fmm-lille | 899,5M | 142,8M | 18,0M | 18,0M |
| 2,8599 | ⟨2,10,11⟩ | 171 | no | R/Q/Z | Fmm-lille | 1,5G | 269,1M | 33,9M | 33,9M |
| 2,8603 | ⟨2,4,9⟩ | 59 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 115,3M | 1,5M | 187,0k | 187,0k |
| 2,8605 | ⟨2,10,12⟩ | 186 | no | R/Q/Z | Fmm-lille | 799,5M | 156,1M | 19,5M | 19,5M |
| 2,8605 | ⟨3,11,11⟩ | 276 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 636,3M | 94,1M | 94,1M |
| 2,8608 | ⟨3,11,12⟩ | 300 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 4,2G | 369,1M | 54,2M | 54,2M |
| 2,8609 | ⟨3,3,11⟩ | 80 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 115,3M | 1,0M | 179,6k | 179,6k |
| 2,8610 | ⟨2,7,8⟩ | 90 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,5G | 45,5M | 5,8M | 5,8M |
| 2,8613 | ⟨2,7,11⟩ | 122 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,0G | 62,2M | 7,9M | 7,9M |
| 2,8613 | ⟨2,11,12⟩ | 204 | no | R/Q/Z | Perminov (tensor decomposition) | 559,7M | 123,0M | 15,5M | 15,5M |
| 2,8615 | ⟨3,3,12⟩ | 87 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 62,9M | 608,3k | 103,5k | 103,5k |
| 2,8621 | ⟨2,12,12⟩ | 222 | no | R/Q/Z | Fmm-lille | 305,3M | 71,4M | 8,9M | 8,9M |
| 2,8628 | ⟨3,7,7⟩ | 117 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 54,7M | 8,2M | 8,2M |
| 2,8628 | ⟨2,7,9⟩ | 101 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,6G | 72,2M | 9,2M | 9,2M |
| 2,8629 | ⟨2,3,9⟩ | 45 | no | R/Q/Z | Perminov (tensor decomposition) | 24,7M | 405,6k | 59,5k | 59,5k |
| 2,8629 | ⟨2,8,12⟩ | 151 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 799,5M | 90,0M | 11,3M | 11,3M |
| 2,8631 | ⟨2,8,11⟩ | 139 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,5G | 155,2M | 19,5M | 19,5M |
| 2,8632 | ⟨2,7,12⟩ | 133 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 36,1M | 4,6M | 4,6M |
| 2,8637 | ⟨11,11,11⟩ | 960 | yes | R/Q/Z | Perminov (FastMatrixMultiplication) | 513,5G | 386,1G | 49,3G | -- |
| 2,8639 | ⟨3,7,9⟩ | 149 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 12,4G | 216,5M | 32,1M | 32,1M |
| 2,8642 | ⟨3,11,11⟩ | 278 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 7,7G | 636,3M | 94,1M | 94,1M |
| 2,8643 | ⟨2,8,9⟩ | 115 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,4G | 180,1M | 22,6M | 22,6M |
| 2,8645 | ⟨2,7,10⟩ | 112 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,5G | 78,9M | 10,0M | 10,0M |
| 2,8663 | ⟨2,10,11⟩ | 173 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,5G | 269,1M | 33,9M | 33,9M |
| 2,8663 | ⟨2,10,12⟩ | 188 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 799,5M | 156,1M | 19,5M | 19,5M |
| 2,8664 | ⟨2,3,10⟩ | 50 | no | R/Q/Z | Perminov (tensor decomposition) | 22,0M | 443,5k | 64,7k | 64,7k |
| 2,8665 | ⟨2,11,12⟩ | 206 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 123,0M | 15,5M | 15,5M |
| 2,8671 | ⟨2,9,10⟩ | 143 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,4G | 312,3M | 39,3M | 39,3M |
| 2,8681 | ⟨2,8,10⟩ | 128 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,1G | 196,9M | 24,6M | 24,6M |
| 2,8684 | ⟨2,9,11⟩ | 157 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,6G | 246,2M | 31,2M | 31,2M |
| 2,8689 | ⟨2,7,9⟩ | 102 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,6G | 72,2M | 9,2M | 9,2M |
| 2,8691 | ⟨2,11,12⟩ | 207 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 123,0M | 15,5M | 15,5M |
| 2,8694 | ⟨2,3,11⟩ | 55 | no | R/Q/Z | Perminov (tensor decomposition) | 15,4M | 349,5k | 51,3k | 51,3k |
| 2,8695 | ⟨2,8,9⟩ | 116 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,4G | 180,1M | 22,6M | 22,6M |
| 2,8695 | ⟨2,10,11⟩ | 174 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,5G | 269,1M | 33,9M | 33,9M |
| 2,8696 | ⟨2,9,12⟩ | 171 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 899,5M | 142,8M | 18,0M | 18,0M |
| 2,8711 | ⟨2,9,10⟩ | 144 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 2,4G | 312,3M | 39,3M | 39,3M |
| 2,8717 | ⟨2,11,12⟩ | 208 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 559,7M | 123,0M | 15,5M | 15,5M |
| 2,8720 | ⟨2,9,11⟩ | 158 | no | R/Q/Z | Perminov (FastMatrixMultiplication) | 1,6G | 246,2M | 31,2M | 31,2M |
| 2,8721 | ⟨2,3,12⟩ | 60 | no | R/Q/Z | Perminov (tensor decomposition) | 8,4M | 202,8k | 29,6k | 29,6k |
| 2,8739 | ⟨2,2,6⟩ | 21 | no | R/Q/Z | AlphaTensor 2022 | 1,1M | 8,7k | 1,1k | 1,1k |
| 2,8844 | ⟨2,2,8⟩ | 28 | no | R/Q/Z | AlphaTensor 2022 | 2,9M | 85,2k | 10,7k | 10,7k |
| 2,8865 | ⟨10,10,10⟩ | 770 | yes | R/Q/Z | Dis09-Q | 1497,2G | 788,9G | 98,6G | -- |
| 2,8914 | ⟨2,2,10⟩ | 35 | no | R/Q/Z | Perminov (tensor decomposition) | 2,9M | 147,8k | 18,5k | 18,5k |
| 2,8945 | ⟨2,2,5⟩ | 18 | no | F2 | Hopcroft-Kerr 1971 | 465,9k | 2,6k | 336 | 336 |
| 2,8945 | ⟨2,2,5⟩ | 18 | no | R/Q/Z | AlphaTensor 2022 | 465,9k | 2,6k | 336 | 336 |
| 2,8950 | ⟨2,2,3⟩ | 11 | no | F2 | Hopcroft-Kerr 1971 | 30,7k | 192 | 28 | 28 |
| 2,8950 | ⟨2,2,3⟩ | 11 | no | R/Q/Z | AlphaTensor 2022 | 30,7k | 192 | 28 | 28 |
| 2,8965 | ⟨2,2,12⟩ | 42 | no | R/Q/Z | Perminov (tensor decomposition) | 1,1M | 67,6k | 8,4k | 8,4k |
| 2,8980 | ⟨2,2,7⟩ | 25 | no | R/Q/Z | AlphaTensor 2022 | 2,1M | 34,2k | 4,3k | 4,3k |
| 2,9014 | ⟨2,2,9⟩ | 32 | no | R/Q/Z | Perminov (tensor decomposition) | 3,3M | 135,2k | 17,0k | 17,0k |
| 2,9044 | ⟨2,2,11⟩ | 39 | no | R/Q/Z | Perminov (tensor decomposition) | 2,1M | 116,5k | 14,7k | 14,7k |
| 2,9197 | ⟨11,11,11⟩ | 1098 | yes | R/Q/Z | Dis09-Q | 513,5G | 386,1G | 49,3G | -- |
| 2,9443 | ⟨8,8,8⟩ | 456 | yes | R/Q/Z | Dis09-Q | 1497,2G | 151,2G | 18,9G | -- |
| 2,9737 | ⟨9,9,9⟩ | 688 | yes | R/Q/Z | Dis09-Q | 2131,7G | 603,6G | 76,8G | -- |
| 3,0541 | ⟨6,6,6⟩ | 238 | yes | R/Q/Z | Dis09-Q | 83,3G | 162,8M | 20,3M | -- |
| 3,0660 | ⟨7,7,7⟩ | 390 | yes | R/Q/Z | Dis09-Q | 513,5G | 9,7G | 1,3G | -- |
| 3,2536 | ⟨5,5,5⟩ | 188 | yes | R/Q/Z | Dis09-Q | 6,0G | 4,1M | 592,7k | 102,3k |
| 3,3219 | ⟨4,4,4⟩ | 100 | yes | R/Q/Z | Dis09-Q | 175,6M | 85,2k | 10,6k | 2,0k |
