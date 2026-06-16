# Catalog discrepancies vs fmm-lille

Generated from `references/fmm-lille-catalog.json` (main fmm-lille catalog page).
Compares 10540 of our scheme files against fmm-lille's best-rank table.

- **Sub-optimal** (we have a higher rank than fmm-lille): 2531
- **Matched** (same rank): 2726
- **Better** (we have a lower rank than fmm-lille): 169
- **Missing** (fmm-lille has it; we don't): 0

## Segmentation by atom vs composed (#198)

`atom` = our best entry is a primitive (explicit import or a formula-constructor ref); `composed` = we built it by Kron/concat/recombination/etc. Losses concentrate in *composed* entries above `MATERIALISE_MAX_DIM=16` — the ≤16 band is atom-dominated and solved.

| band | our entry | tie | better | worse |
|---|---|---:|---:|---:|
| <=16 | atom | 653 | 8 | 3 |
| <=16 | composed | 2 | 0 | 0 |
| 17-32 | atom | 7 | 0 | 5 |
| 17-32 | composed | 2064 | 161 | 2523 |

## Sub-optimal (top 50)

Formats where fmm-lille knows of a strictly lower-rank algorithm than what we have.

| format | our rank | fmm rank | Δ | our count | fmm refs | fmm detail |
|---|---|---|---|---|---|---|
| ⟨29,30,30⟩ | 13902 | 12588 | +1314 | 1 | — | [link](https://fmm.univ-lille.fr/29x30x30.html) |
| ⟨27,28,28⟩ | 11718 | 10413 | +1305 | 2 | — | [link](https://fmm.univ-lille.fr/27x28x28.html) |
| ⟨31,32,32⟩ | 16002 | 15006 | +996 | 1 | — | [link](https://fmm.univ-lille.fr/31x32x32.html) |
| ⟨25,26,26⟩ | 9506 | 8552 | +954 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x26.html) |
| ⟨29,30,31⟩ | 14302 | 13458 | +844 | 2 | — | [link](https://fmm.univ-lille.fr/29x30x31.html) |
| ⟨27,28,29⟩ | 12010 | 11169 | +841 | 2 | — | [link](https://fmm.univ-lille.fr/27x28x29.html) |
| ⟨27,28,32⟩ | 13124 | 12339 | +785 | 1 | — | [link](https://fmm.univ-lille.fr/27x28x32.html) |
| ⟨29,29,30⟩ | 13309 | 12526 | +783 | 3 | — | [link](https://fmm.univ-lille.fr/29x29x30.html) |
| ⟨31,31,32⟩ | 15700 | 14940 | +760 | 2 | — | [link](https://fmm.univ-lille.fr/31x31x32.html) |
| ⟨30,32,32⟩ | 15595 | 14874 | +721 | 2 | drevet:2011a | [link](https://fmm.univ-lille.fr/30x32x32.html) |
| ⟨29,32,32⟩ | 15432 | 14738 | +694 | 1 | drevet:2011a | [link](https://fmm.univ-lille.fr/29x32x32.html) |
| ⟨25,26,27⟩ | 9878 | 9202 | +676 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x27.html) |
| ⟨29,30,32⟩ | 14576 | 13908 | +668 | 2 | — | [link](https://fmm.univ-lille.fr/29x30x32.html) |
| ⟨27,28,31⟩ | 12739 | 12072 | +667 | 1 | — | [link](https://fmm.univ-lille.fr/27x28x31.html) |
| ⟨28,30,30⟩ | 13088 | 12464 | +624 | 3 | drevet:2011a | [link](https://fmm.univ-lille.fr/28x30x30.html) |
| ⟨25,26,28⟩ | 10154 | 9540 | +614 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x28.html) |
| ⟨25,28,28⟩ | 10803 | 10206 | +597 | 1 | drevet:2011a | [link](https://fmm.univ-lille.fr/25x28x28.html) |
| ⟨25,25,26⟩ | 9082 | 8498 | +584 | 1 | — | [link](https://fmm.univ-lille.fr/25x25x26.html) |
| ⟨27,30,30⟩ | 12880 | 12336 | +544 | 2 | drevet:2011a | [link](https://fmm.univ-lille.fr/27x30x30.html) |
| ⟨26,28,28⟩ | 10858 | 10326 | +532 | 3 | drevet:2011a | [link](https://fmm.univ-lille.fr/26x28x28.html) |
| ⟨25,26,30⟩ | 10746 | 10215 | +531 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x30.html) |
| ⟨23,27,27⟩ | 9369 | 8856 | +513 | 1 | — | [link](https://fmm.univ-lille.fr/23x27x27.html) |
| ⟨27,28,30⟩ | 12062 | 11561 | +501 | 2 | — | [link](https://fmm.univ-lille.fr/27x28x30.html) |
| ⟨25,26,29⟩ | 10486 | 9986 | +500 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x29.html) |
| ⟨24,26,26⟩ | 8918 | 8444 | +474 | 1 | drevet:2011a | [link](https://fmm.univ-lille.fr/24x26x26.html) |
| ⟨24,26,27⟩ | 9302 | 8840 | +462 | 1 | — | [link](https://fmm.univ-lille.fr/24x26x27.html) |
| ⟨29,29,32⟩ | 14259 | 13802 | +457 | 2 | — | [link](https://fmm.univ-lille.fr/29x29x32.html) |
| ⟨23,26,27⟩ | 9077 | 8651 | +426 | 1 | — | [link](https://fmm.univ-lille.fr/23x26x27.html) |
| ⟨27,27,28⟩ | 10800 | 10384 | +416 | 1 | — | [link](https://fmm.univ-lille.fr/27x27x28.html) |
| ⟨29,29,31⟩ | 13777 | 13367 | +410 | 3 | — | [link](https://fmm.univ-lille.fr/29x29x31.html) |
| ⟨23,27,28⟩ | 9682 | 9278 | +404 | 1 | — | [link](https://fmm.univ-lille.fr/23x27x28.html) |
| ⟨27,29,29⟩ | 12352 | 11952 | +400 | 2 | — | [link](https://fmm.univ-lille.fr/27x29x29.html) |
| ⟨21,22,22⟩ | 5871 | 5476 | +395 | 1 | — | [link](https://fmm.univ-lille.fr/21x22x22.html) |
| ⟨25,26,31⟩ | 11078 | 10684 | +394 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x31.html) |
| ⟨25,26,32⟩ | 11343 | 10954 | +389 | 1 | — | [link](https://fmm.univ-lille.fr/25x26x32.html) |
| ⟨22,27,28⟩ | 9380 | 8997 | +383 | 1 | — | [link](https://fmm.univ-lille.fr/22x27x28.html) |
| ⟨23,26,26⟩ | 8711 | 8332 | +379 | 1 | drevet:2011a | [link](https://fmm.univ-lille.fr/23x26x26.html) |
| ⟨22,26,27⟩ | 8783 | 8408 | +375 | 1 | — | [link](https://fmm.univ-lille.fr/22x26x27.html) |
| ⟨27,27,29⟩ | 11474 | 11113 | +361 | 3 | — | [link](https://fmm.univ-lille.fr/27x27x29.html) |
| ⟨29,31,32⟩ | 15060 | 14705 | +355 | 2 | drevet:2011a | [link](https://fmm.univ-lille.fr/29x31x32.html) |
| ⟨25,27,28⟩ | 10518 | 10164 | +354 | 1 | — | [link](https://fmm.univ-lille.fr/25x27x28.html) |
| ⟨24,26,28⟩ | 9590 | 9240 | +350 | 1 | — | [link](https://fmm.univ-lille.fr/24x26x28.html) |
| ⟨27,29,30⟩ | 12651 | 12305 | +346 | 1 | drevet:2011a | [link](https://fmm.univ-lille.fr/27x29x30.html) |
| ⟨28,30,31⟩ | 13639 | 13304 | +335 | 3 | — | [link](https://fmm.univ-lille.fr/28x30x31.html) |
| ⟨29,31,31⟩ | 14692 | 14357 | +335 | 2 | — | [link](https://fmm.univ-lille.fr/29x31x31.html) |
| ⟨23,27,29⟩ | 9955 | 9625 | +330 | 2 | — | [link](https://fmm.univ-lille.fr/23x27x29.html) |
| ⟨27,29,32⟩ | 13510 | 13187 | +323 | 1 | — | [link](https://fmm.univ-lille.fr/27x29x32.html) |
| ⟨26,28,31⟩ | 12243 | 11930 | +313 | 1 | — | [link](https://fmm.univ-lille.fr/26x28x31.html) |
| ⟨23,26,28⟩ | 9356 | 9054 | +302 | 1 | — | [link](https://fmm.univ-lille.fr/23x26x28.html) |
| ⟨24,26,29⟩ | 9878 | 9576 | +302 | 1 | — | [link](https://fmm.univ-lille.fr/24x26x29.html) |

## Cases where we have a lower rank than fmm-lille

Likely indicates fmm-lille is missing a recent result — worth flagging upstream.

| format | our rank | fmm rank | Δ | our count |
|---|---|---|---|---|
| ⟨23,25,31⟩ | 9667 | 9776 | -109 | 2 |
| ⟨24,25,32⟩ | 10176 | 10268 | -92 | 1 |
| ⟨19,30,30⟩ | 9248 | 9329 | -81 | 1 |
| ⟨11,26,32⟩ | 5343 | 5415 | -72 | 1 |
| ⟨23,25,32⟩ | 9948 | 10018 | -70 | 1 |
| ⟨15,20,31⟩ | 5247 | 5316 | -69 | 2 |
| ⟨18,29,30⟩ | 8469 | 8535 | -66 | 1 |
| ⟨22,29,32⟩ | 10980 | 11046 | -66 | 1 |
| ⟨17,29,30⟩ | 8139 | 8204 | -65 | 2 |
| ⟨24,25,31⟩ | 9909 | 9972 | -63 | 1 |
| ⟨26,32,32⟩ | 13896 | 13957 | -61 | 2 |
| ⟨13,26,32⟩ | 6242 | 6296 | -54 | 1 |
| ⟨18,30,31⟩ | 8988 | 9042 | -54 | 2 |
| ⟨21,24,32⟩ | 8640 | 8694 | -54 | 1 |
| ⟨23,23,32⟩ | 9102 | 9156 | -54 | 1 |
| ⟨17,31,31⟩ | 8832 | 8885 | -53 | 2 |
| ⟨12,26,32⟩ | 5586 | 5636 | -50 | 2 |
| ⟨22,25,31⟩ | 9424 | 9473 | -49 | 2 |
| ⟨11,26,31⟩ | 5218 | 5264 | -46 | 2 |
| ⟨14,15,29⟩ | 3591 | 3636 | -45 | 2 |
| ⟨21,29,32⟩ | 10612 | 10657 | -45 | 1 |
| ⟨13,27,32⟩ | 6494 | 6536 | -42 | 1 |
| ⟨21,21,21⟩ | 5202 | 5240 | -38 | 3 |
| ⟨28,31,32⟩ | 14439 | 14477 | -38 | 3 |
| ⟨11,25,32⟩ | 5157 | 5193 | -36 | 1 |
| ⟨13,25,32⟩ | 6012 | 6048 | -36 | 1 |
| ⟨13,26,31⟩ | 6101 | 6137 | -36 | 1 |
| ⟨11,27,32⟩ | 5577 | 5612 | -35 | 1 |
| ⟨12,13,16⟩ | 1509 | 1544 | -35 | 6 |
| ⟨11,25,31⟩ | 5012 | 5045 | -33 | 2 |

