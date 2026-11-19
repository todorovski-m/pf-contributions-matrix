# Screening of 15 MATPOWER cases under rule Q3 version 2.0 (brief T16, S5)

Earlier outcome: the scan of brief T14f (rule Q3 version 1.0). New outcome: `benchmark/timing/scan.py` of this version on the exports of `benchmark/timing/export_large.sh` (MATPOWER 7.1). Status per case and kind: usable (with the spectral radius of N when the digraph is cyclic) or the error identifier with the spectral radius or the buses named in the message. No case was tuned. Made by `reports/T16_scripts/s5_table.py` (outside this repository).

Memory limit for the dense tables: 4294967296 bytes.

| case | kind | earlier outcome (rule Q3 v1.0, T14f) | new outcome (rule Q3 v2.0) | changed | dense table bytes (new) | within the memory limit |
|---|---|---|---|---|---|---|
| case145 | apfc | pfc:apfc:spectralRadius, rho = 1.0000000000000011 | pfc:apfc:spectralRadius, rho = 1.0000000000000011 | no | 686760 | True |
| case145 | rpfc | pfc:rpfc:spectralRadius, rho = 1.6956089156006535 | usable (cyclic, rho(N) = 0.58087477604349269) | yes | 8736768 | True |
| case_ACTIVSg200 | apfc | usable | usable | no | 279072 | True |
| case_ACTIVSg200 | rpfc | usable | usable | no | 2771280 | True |
| case1888rte | apfc | usable (cyclic, rho(N) = 0.15079631981127514) | usable (cyclic, rho(N) = 0.15079631981127514) | no | 25833456 | True |
| case1888rte | rpfc | usable (cyclic, rho(N) = 0.97928887465785142) | usable (cyclic, rho(N) = 0.83808745727218592) | no | 309297664 | True |
| case1951rte | apfc | usable (cyclic, rho(N) = 0.10090282235667325) | usable (cyclic, rho(N) = 0.10090282235667325) | no | 32928768 | True |
| case1951rte | rpfc | pfc:rpfc:spectralRadius, rho = 1.0182112890707311 | usable (cyclic, rho(N) = 0.87139170553102063) | yes | 319359680 | True |
| case2736sp | apfc | usable | usable | no | 25609272 | True |
| case2736sp | rpfc | usable (cyclic, rho(N) = 0.55692767858691872) | usable (cyclic, rho(N) = 0.60941553331621789) | no | 546612480 | True |
| case2737sop | apfc | usable | usable | no | 22007520 | True |
| case2737sop | rpfc | usable (cyclic, rho(N) = 0.6735753424072547) | usable (cyclic, rho(N) = 0.60587890241535025) | no | 548052920 | True |
| case2746wop | apfc | usable | usable | no | 44522720 | True |
| case2746wop | rpfc | usable (cyclic, rho(N) = 0.66252953469562814) | usable (cyclic, rho(N) = 0.63819122050522703) | no | 575719680 | True |
| case2746wp | apfc | usable | usable | no | 46222704 | True |
| case2746wp | rpfc | usable (cyclic, rho(N) = 0.68962469382218405) | usable (cyclic, rho(N) = 0.62622485110809456) | no | 579377040 | True |
| case2848rte | apfc | usable | usable | no | 67014416 | True |
| case2848rte | rpfc | pfc:rpfc:zeroInflow, buses 686, 689, | usable (cyclic, rho(N) = 0.78764050840134947) | yes | 676260816 | True |
| case3012wp | apfc | pfc:apfc:negativeLoss | pfc:apfc:negativeLoss | no | 45421328 | True |
| case3012wp | rpfc | pfc:rpfc:spectralRadius, rho = 1.9440186823873005 | usable (cyclic, rho(N) = 0.56386563792772848) | yes | 661586064 | True |
| case3375wp | apfc | pfc:apfc:negativeLoss | pfc:apfc:negativeLoss | no | 64638176 | True |
| case3375wp | rpfc | usable (cyclic, rho(N) = 0.55818663287206272) | usable (cyclic, rho(N) = 0.58574382149780446) | no | 896593600 | True |
| case6468rte | apfc | usable (cyclic, rho(N) = 0.95909992048608861) | usable (cyclic, rho(N) = 0.95909992048608861) | no | 161028000 | True |
| case6468rte | rpfc | pfc:rpfc:zeroInflow, buses 1354, | usable (cyclic, rho(N) = 0.85244411117667229) | yes | 3480685536 | True |
| case6470rte | apfc | usable (cyclic, rho(N) = 0.9590995996596432) | usable (cyclic, rho(N) = 0.9590995996596432) | no | 257334672 | True |
| case6470rte | rpfc | pfc:rpfc:spectralRadius, rho = 4.2960374900025702 | usable (cyclic, rho(N) = 0.79346138584134074) | yes | 3552915664 | True |
| case6495rte | apfc | usable (cyclic, rho(N) = 0.95909949199418298) | usable (cyclic, rho(N) = 0.95909949199418298) | no | 244553280 | True |
| case6495rte | rpfc | pfc:rpfc:spectralRadius, rho = 1.1951851984685702 | usable (cyclic, rho(N) = 0.79357231587679622) | yes | 3568524960 | True |
| case6515rte | apfc | usable (cyclic, rho(N) = 0.95909936842954902) | usable (cyclic, rho(N) = 0.95909936842954902) | no | 243846976 | True |
| case6515rte | rpfc | pfc:rpfc:spectralRadius, rho = 1.2418094435886058 | usable (cyclic, rho(N) = 0.79138328894829679) | yes | 3569903736 | True |

combinations whose outcome changed: 8
