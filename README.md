# Characteristic vectors and strongly regular graphs at Krein equality

Draft preprint (`krein_unified.tex`, `krein_unified.pdf`) and the scripts and data behind every count in it.

A primitive strongly regular graph with integral eigenvalues and vanishing Krein parameter q^s_ss gives
an integral tight frame whose first and third moments vanish. The vanishing third moment lifts a
congruence for characteristic vectors from modulus 8 to 16, which forces every maximal integral lattice
containing the frame to be unimodular at 2. Van der Blij's congruence then gives a congruence modulo 4
in the parameters (Theorem 1.1). Consequences include:

- no triangle-free Krein graph Kre(t) for t ≡ 4 (mod 8), in particular no srg(784,116,0,20), and no
  3-(t(t^2+3t+1), t(t+1), t-1) design for these t (e.g. no 3-(116,20,3) design);
- no srg((t^2+2t-1)(t^2+3t+1), t^2(t+2), 0, t^2) for t ≡ 4 (mod 8), in particular no srg(667,96,0,16):
  a Norton-algebra argument shows that every such graph comes from a design as above;
- two further infinite series, starting with srg(392,115,18,40) and srg(1080,364,88,140);
- further families from a projection away from one vertex and from the second moment (including
  Biggs's triangle-free families with mu = 2, 4, 6 and a Pell family), and exclusions of
  distance-regular graphs.

Environment: Python 3 with sympy, numpy and matplotlib; PARI/GP for the lattice controls.

## Smith census (Section 9, Table 1, Figure 2)

    python -B scan_complete.py 135 smith_complete_r135.tsv      # about 3 s

Writes one row per Smith parameter set with r <= 135 (2320 rows) and prints the summary:

| quantity | value |
|---|---|
| Theorem 1.1 (excl / pass / undet / na) | 68 / 132 / 986 / 1134 |
| Theorem 6.1 (excl / pass / undet / na) | 12 / 64 / 488 / 1756 |
| excluded by 1.1 or 6.1 | 77 |
| ... passing the feasibility checks (`checks` column) | 71 |
| ... of these with a Prop. 4.4 exclusion (graph or subconstituent) | 28 |
| ... with no Prop. 4.4 exclusion | 43 |
| Prop. 4.4 exclusions missed by 1.1/6.1 | 0 |
| realized rows (table snapshot) / excluded / failing the checks | 8 / 0 / 0 |
| tight 4-design sets / excluded | 135 / 0 |

For the 68 members of Kre(t) whose L = t has a prime divisor ≡ 3 (mod 4), Delta is computed from the
neighbor basis, as in the proof of Theorem 1.2(a); 11 of them are excluded.

`feasible` is parent feasibility (integrality, Krein conditions, absolute bound). `checks` is the
predicate used in the paper ("passes the feasibility checks"): `feasible`, integral forced local gaps
d and e-d, and each nontrivial subconstituent feasible or a disjoint union of equal cliques (or the
complement of one). `S2P`/`S2S` record a Prop. 4.4 exclusion of the graph or its complement, or of
either subconstituent or its complement.

An independently written enumeration of all Smith parameter sets with n <= 10^7:

    python -B scan_krein.py 10000000 400 hits.tsv

Figure 2: `cd figures && MPLCONFIGDIR=<writable dir> python -B reach_plot.py` (reads the TSV).

## Brouwer's table

`data/brouwer_srg_table_20261118.csv` is a snapshot (18 November 2026) of the status, parameters and
comments of the 4538 rows of A. E. Brouwer's table of strongly regular graphs,
<https://aeb.win.tue.nl/graphs/srg/>. Prop. 4.4 fails directly on 124 rows and in either orientation
on 154 rows (77 complementary pairs); all are marked nonexistent. All tests of Theorems 1.1 and 6.1 and
Prop. 4.4 together hit 156 nonexistent rows plus the six open orientations of 392/784/1080;
Theorem 1.2(d) adds the two open rows on 667 vertices.

## Second subconstituent of Kre(t) (Section 5.3, Theorem 1.2(d))

    python -B second_subconstituent.py   # about 1 minute

Checks the homomorphism expansion of the sum Phi against brute force on five triangle-free strongly
regular graphs, the identities of the proof of Theorem 5.1 in Q(t), and, for t = 2, recovers the
3-(22,6,1) design from a shuffled srg(77,16,0,4) through the Norton algebra.

## Distance-regular graphs (Section 8, Tables 2 and 3)

The machine-readable BCN tables (files IA3P, IA4a, IA5a, IA6b, ..., IA17b) are available from Brouwer's
pages, <https://www.win.tue.nl/~aeb/drg/>. They are not redistributed here; put them in a directory
`BCN_DIR` and check them against `data/bcn_files.sha256` (`cd BCN_DIR && sha256sum -c .../bcn_files.sha256`).

    python -B frame_check.py BCN_DIR > drg_hits.tsv
      # controls: 223 arrays, 769 rational eigenspaces, 0 contradictions
      # arrays 5475, failing eigenspaces 23, failing arrays 13
    python -B audit_projector_sums.py BCN_DIR
      # 720282 single-residual sums + 240 split sums = 720522; 0 additional failing sums

`frame_check.py` is an independent implementation of Theorem 3.1 and Lemma 3.2; `distance_arrays.py`
and `frame_test.py` are the original implementation used by the projector-sum audit. Control inventory
(`frame_check.controls`): Hamming, Johnson, Grassmann, dual polar, bilinear forms, halved cubes,
generalized hexagons and octagons, and ten sporadic graphs.

Derived graphs in Table 2: for bipartite arrays the halved graph is (A^2 - kI)/c_2 on one half; its
complement at 1014 / 5046 / 11774 is (507,44,1,4) / (2523,182,1,14) / (5887,216,5,8), each also failing
Prop. 4.4. At 750, A_1 + A_3 is (750,210,55,60). At 616 a fusion gives (616,165,44,44), which is open;
Prop. 4.4 does not apply to it or its complement.

## Families (Section 7)

    python -B srg_check.py       # Pell (12 members), Theorem 7.2 and Corollary 7.3 for q < 3000

## Lattice controls (Section 9)

    cd controls && python -B graphs_gram.py
    echo 'run(G,2)'    | gp -q G_hs.gp twoadic.gp            # Higman-Sims
    echo 'run(G,21/5)' | gp -q -s 2G G_q5.gp twoadic.gp      # Q^-(5,5), a few minutes

Each computes a maximal 2-integral lattice containing the vectors of b(rI-A+uJ), confirms that it is
unimodular at 2, and prints w^2, n/F and m+1-Delta modulo 4 (2, 2, 2 and 1, 1, 1).
