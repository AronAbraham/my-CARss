Perfect — it's genuinely sparse. `issparse = 1`, `nnz = 118664` (≈4 nonzeros per row on average, exactly what a cortical graph Laplacian should be). This makes T3 completely tractable — the graph term is a cheap sparse operation, no 7 GB dense matrix, no slow inversions.

[Certain] So the plan is clean: load the sparse L, convert to scipy.sparse, write graph-regularized ALS where the spatial-factor update includes the `μL` penalty via a sparse solve. Let me write T3.

One thing I need to handle: [Certain] the Laplacian is in `Laplacian.mat`, which is old-format (scipy-loadable). scipy loads MATLAB sparse matrices directly as `scipy.sparse` arrays — so this is smooth.Stop — I need to be honest before you run this, because writing T3 surfaced a genuine problem I should have caught earlier, and it's important.

[Certain] **There's a dimension mismatch that breaks the graph regularization as designed.** The graph Laplacian `L` is 29666×29666 — it's over **source space (dipoles)**. But your wavelet tensor is `75 × 300 × 70` — its spatial mode is **75 channels (sensor space)**, not 29666 dipoles. The paper's `Tr(IᵀLI)` term applies to the *source-space* spatial factor. You can't apply a 29666-dipole Laplacian to a 75-channel factor — the dimensions don't match, and it's not just a technicality, it's a modeling issue.

The script detects this and falls back to plain CPD so it still runs, but I owe you the real explanation of what this means:

[Certain] **The paper builds its tensor in source space, or maps to source space before the graph term.** Look back at the tensor-CARSS paper's structure: it does CPD to get spatial *components*, and the graph regularization enforces smoothness on those components *in the brain*. For that to work, the spatial factor must be over dipoles (29666), not channels (75). Your current pipeline builds a *sensor-space* tensor (channels × time × freq), so its spatial factor is channels — the graph term can't apply.

[Likely] So there are two ways the paper's method actually works, and I need to figure out which before T3 is meaningful:

1. **The tensor is built in source space** — i.e., the EEG is first mapped to the 29666 dipoles (via some inverse), *then* the tensor is channels→dipoles. But that's circular (you'd need source localization before the tensor).

2. **More likely: the CPD is on the sensor tensor (channels), giving channel-space spatial components, and CARSS then maps *those* to source space** — and the graph regularization is applied during the CARSS/inverse step, not the CPD. In this reading, the graph term regularizes the *final source estimate*, not the CPD factors.

[Guessing] I genuinely am not certain which the paper intends — the equation (5) shows `Tr(IᵀLI)` on the CPD factor `I`, suggesting reading 1, but that requires source-space tensors which the paper doesn't clearly describe building. This is exactly the kind of gap where the paper states the objective but not the full pipeline.

So here's the honest situation: **T3 as literally written (graph term on the CPD spatial factor) requires a source-space spatial factor, which your sensor-space tensor doesn't have.** I don't want to hand you code that "runs" via a fallback but doesn't actually do what the paper means — that's the kind of thing that falls apart in a viva.

Before proceeding, I need to resolve this by reading the paper more carefully. The key question: **does the paper's CPD produce channel-space or source-space spatial components, and where exactly does the Laplacian apply?**

Let me pull the relevant section from the paper you uploaded:Now I have the definitive answer, and it confirms the problem I flagged — but also resolves it. This is important, so let me be precise.

[Certain] Read equation 5 carefully: the spatial factor `Iₖ` has size **M × 1**, where the paper defines **M as the number of voxels** (from line 76: "M is the number of voxels"). So the paper's tensor `Γ` is built in **source space** — its first mode is voxels (M ≈ 29666), NOT channels (75). The graph Laplacian (M×M) applies to the voxel-space spatial factor. That matches: L is 29666×29666, Iₖ is 29666×1. Dimensions agree.

**So your pipeline and the paper diverge at a fundamental point:**

[Certain] Your T1 tensor is `75 × 300 × 70` — **sensor space** (channels). The paper's tensor `Γ(r, τ, f)` is **source space** — `voxels × time × freq`. Look at eq 3: `Γ(r, τ, f)` where `r` is spatial *position* (source), not channel. The paper builds the tensor over *source locations*, not electrodes.

This is the real modeling difference, and it changes the whole architecture. Here's what it means:

[Likely] **The paper's pipeline is:**
1. Map EEG to source space first (some inverse/imaging step) → get source-space time series
2. Wavelet-transform → source × time × freq tensor
3. Graph-regularized CPD on that (graph term applies, dimensions match)
4. CARSS on the resulting source-space spatial components

[Certain] **Your pipeline is:**
1. Wavelet-transform the *sensor* EEG → channel × time × freq tensor
2. CPD → channel-space spatial components
3. CARSS localizes those to source space

These are genuinely different. In the paper, source localization happens partly *before* the tensor (to get source space) and the CPD works on already-localized data. In yours, the CPD works on sensor data and CARSS does the localization.

[Guessing] Honestly, the paper is a bit circular/underspecified here — if you need source space *before* the CPD, how do you get it without already having done source localization? The paper doesn't clearly say. This is a real gap in the paper, not just your understanding.

**So here's the honest situation and your options:**

[Certain] **The graph-regularized CPD (T3) as the paper specifies requires a source-space tensor, which requires source localization before the CPD.** Your current sensor-space pipeline can't apply the graph term directly. So you have three honest paths:

**Path 1 — build the source-space tensor (faithful to paper).** Add a step: map EEG to source space (e.g., via a simple inverse like minimum-norm), then build the voxel×time×freq tensor, then graph-CPD applies correctly. But this is a big architectural change, it's partly circular, and the source-space tensor is huge (29666×300×70 = 622M entries — memory-heavy).

**Path 2 — apply the graph term in the CARSS/inverse step instead.** Keep your sensor-space CPD, and add the graph regularization where it naturally fits — in the source estimation. This is arguably more sensible than the paper's formulation and avoids the circularity.

**Path 3 — document the difference, proceed with plain CPD → CARSS.** Your sensor-space CPD + CARSS is a *coherent, working* pipeline. Note honestly that the paper's graph term applies in source space and your architecture places localization differently, so you implement the tensor-CARSS *core* (wavelet + CPD + CARSS) without the source-space graph term, and discuss it as a deviation.

[Likely] My honest recommendation: **Path 3 for now, consider Path 1 later if time allows.** Reasoning: you have a working, validated wavelet→CPD→CARSS pipeline within reach (T1 ✓, T2 ✓, CARSS ✓). The graph regularization requires a source-space tensor that's a major architectural change with real memory/circularity issues. Getting the core pipeline complete and working — then honestly documenting that graph regularization requires source-space tensors which you note as a limitation/future-work — is a defensible, complete project. Chasing the source-space tensor risks the STSI trap: a big, memory-heavy, underspecified addition that could eat your timeline.

[Certain] This is exactly the kind of decision where I won't push you into the risky path. You've built something real and working. The graph term is a refinement the paper specifies ambiguously, and implementing it faithfully requires re-architecting to source space. That's a genuine research decision, not a quick add.

So the honest question for you: **do you need the graph regularization to be faithful to the paper (Path 1, big change, risky), or is the tensor-CARSS core (wavelet + CPD + CARSS) sufficient with graph regularization documented as a deviation (Path 3, working now)?**

This partly depends on your guide — did they specifically require the graph-regularized part, or the tensor-CARSS approach generally? And it depends on your timeline. But I'd strongly caution against Path 1 unless there's a clear reason, because it re-introduces exactly the architectural risk you escaped.

What's your read — must-have graph regularization (accept the source-space rebuild), or core tensor-CARSS with the graph term as documented future work?