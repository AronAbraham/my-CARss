Step 1 passes — cleanly and unambiguously. That eigenvalue spectrum is textbook:

- Eigenvalue 1: exactly 1.0
- Eigenvalue 2 onward: ~1e-16 (machine zero)

[Certain] That's a *perfect* single-source signature. One dominant eigenvalue, everything else at numerical precision zero — meaning the covariance is rank-1, exactly what one coherent source with one time course produces. Estimated n = 1, correct. The subspace split works. The v7.3 fallback to mat73 worked too. No issues at all.

One thing worth noting for your understanding: [Likely] the drop to 1e-16 is *because* this is noise-free data. When you later add noise, eigenvalues 2+ will rise from machine-zero to small-but-nonzero, and the "cliff" becomes a "step." The source-counting heuristic (threshold at 1% of max) is designed for exactly that noisy case. On clean data it's trivially obvious; on noisy data it earns its keep. Good to know the clean case is the easy validation, not the real test.

Now **Step 2 — the MUSIC localizer**, which is the real proof. This scans all 29666 dipoles, computes how well each explains the signal subspace, and its peak should land exactly on dipole 29466 (your source). *This* is the gate that proves CARSS's core actually localizes.

Let me write it.Step 2 is in your outputs. This is the real test — run it, and it scans all 29666 dipoles with the MUSIC localizer.

**What you're validating — this is the gate that proves CARSS's core works:**

[Certain] The peak dipole should be **29466** (your source, 0-based) or a dipole essentially at the same position, with **localization error near 0 mm**. The script prints PASS if error < 5 mm. If the peak lands on the true source, MUSIC localization works and CARSS's foundation is proven — we proceed to Step 3 (the recursion for multiple sources).

Two honest predictions about what you'll see:

[Likely] **On clean single-source data, MUSIC should nail it — error ~0 mm, mu at peak ~1.0.** MUSIC is exact in the noise-free single-source case; that's its theoretical property. So I expect a clean PASS. If it *doesn't* land on 29466, that signals a real bug (leadfield/data mismatch, orientation issue) and we debug before going further.

[Guessing] **One thing that could trip it:** the scalar leadfield approximation (Option A) throws away orientation. For a single source this usually still localizes correctly because the *location's* leadfield pattern dominates. But if the source's orientation is unusual, the scalar collapse could blur it slightly. If you get error of, say, 5–15 mm instead of ~0, that's likely the scalar approximation, not a bug — and it's the signal that we need the full 3-orientation version sooner rather than later. Watch for that.

The script also prints the **top 10 dipoles** by localizer value and their distance to truth — that's diagnostic. Even if the #1 peak is right, seeing whether #2–10 are *near* the source (good — they're the neighbors) or scattered (concerning) tells you about the localizer's sharpness.

Run it and paste the full output — the peak dipole, the error, mu at peak, and the top-10 list. Three outcomes:

1. **Error ~0 mm, PASS** → MUSIC works, CARSS core validated → we build Step 3 (RAP-MUSIC recursion for the two-source case).
2. **Error 5–15 mm** → probably the scalar approximation; we discuss whether to upgrade to 3-orientation now or proceed and upgrade later.
3. **Error large / wrong location** → a real bug; we stop and debug before anything else.

This is the moment that tells you whether the whole CARSS path is sound. Step 1 proved the subspace math; Step 2 proves it can actually find a source. Run it.