# Getting Started

1. Modify line 13 to point to your local irt path.
2. Run main.m, which should execute quickly. PERK results should be the same as in the paper, but VPM and PGPM results are poorer because coarse grid search and zero iterations are used by default. 
3. To replicate paper results, comment out lines 139, 142, 155 and uncomment lines 140, 143, 156. 
4. To replicate PERK performance degradation results in the case of mismatch between acquisition design versus sampling distribution support, comment out lines 209, 214 and uncomment lines 210, 215.
