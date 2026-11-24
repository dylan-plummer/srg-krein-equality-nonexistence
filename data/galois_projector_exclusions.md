# Further BCN arrays excluded through rational projector sums (Section 8)

Among the 5475 arrays in the BCN files, the rational projector sums give nine exclusions beyond those of single
rational eigenspaces. Each array would produce a strongly regular graph that is already excluded. Every failing sum
below omits the principal idempotent E_0, so its first moment vanishes, and it is checked with Theorem 3.1 and
Lemma 3.2: (gamma_i) is the primitive integral vector proportional to the sum, indexed by distance; F and m are its
frame constant and rank. Recompute with `python -B list_projector_sums.py` (about one minute).

| n | intersection array | failing sum E_θ + E_θ' (+ …) | (γ_i) | F | m | cubic | would give |
|---|---|---|---|---|---|---|---|
| 300 | {39,30,4;1,5,36} | θ = 9, −6 | (91, 1, 1, −14) | 150 | 182 | no | (300,26,4,2), excluded by Kovaľ |
| 352 | {25,24,2,1;1,2,24,25} | θ = 5, 3; also 3, −5 | (65, 10, −1, −1, 10) | 110 | 208 | no | (352,26,0,2), Corollary 7.3 |
| 784 | {70,69,56,10;1,14,60,70} | θ = −14, −70 | (29, −6, 1, 1, −6) | 196 | 116 | yes | Kre(4) |
| 784 | {115,96,20,1;1,20,96,115} | θ = 23, −25 | (29, 1, 1, −6, −6) | 196 | 116 | yes | Kre(4) |
| 1500 | {15,14,14,11;1,1,4,15} | θ = 5, 0, −15; also 0, −5, −15 | (529, 73, −11, −2, 4) | 750 | 1058 | no | halved: (750,210,55,60), excluded by Kovaľ |
| 12675 | {594,550,60,1;1,30,550,594} | θ = 54, 9; also −11, −56 | (1012, 37, −2, −2, 583) | 2535 | 5060 | no | folded: (4225,594,43,90), Proposition 4.4 |
| 13872 | {285,280,6,1;1,6,280,285} | θ = 19, −21; also −15, −21 | (715, 1, 1, −50, −50) | 1734 | 5720 | no | folded: (6936,285,4,12), Theorem 7.2 |
| 14196 | {165,160,10,1;1,2,160,165} | θ = −11, −17 | (1275, −90, 1, 1, −90) | 2366 | 7650 | no | distances 1 ∪ 4: (14196,170,4,2), Theorem 7.2 |
| 14536 | {169,168,2,1;1,2,168,169} | θ = 13, 11; also 11, −13 | (1105, 78, −1, −1, 78) | 2054 | 7820 | no | folded: (7268,169,0,4), Corollary 7.3 |

Where two sums are listed, the data shown are for the first; `list_projector_sums.py` prints both.
