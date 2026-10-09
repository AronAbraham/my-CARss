Everything's there — `EEG_clean` (75×300, noise-free), `GT_ez_dip` (29467), `GT_ez_pos_mm`, `GT_Location` (3×29666). That's exactly the clean single-source test case we need. And scalar-first is confirmed by your going ahead.

Let me write Step 1 — the covariance and subspace split — in Python. This is the foundation: load your real leadfield, build the scalar version, compute the measurement covariance, eigendecompose, and **show you the eigenvalue drop-off** so you can see it identify one source. That eigenvalue plot is your first validation gate.Step 1 is in your outputs. Run it (edit the two paths at top if needed). It loads your real leadfield, builds the scalar version, computes the covariance, and produces the eigenvalue plot.

**What you're validating — read this carefully, it's the whole point of Step 1:**

[Certain] For a **clean, single source**, the covariance should have essentially **one dominant eigenvalue, then a sharp drop to near-zero.** The plot (`carss_step1_eigenvalues.png`) shows this two ways — log scale (dramatic cliff) and linear (first 20 values). If you see one big value then a cliff, and the script prints "Estimated number of sources: 1" — **Step 1 works, subspace splitting is correct, and we proceed to Step 2 (the MUSIC localizer that actually finds where the source is).**

If it prints something other than 1, or the eigenvalues don't show a clean drop, stop and paste the output — that tells us something's off in the covariance or the data, and we fix it before building anything on top.

Three honest cautions:

[Guessing] **Your leadfield might be a v7.3 .mat file** (it's large — 88998 columns). If `scipy.io.loadmat` fails with a NotImplementedError, the script tells you to `pip install mat73` and it'll auto-retry. Just install it and rerun. Not a real problem, just a format branch.

[Certain] **This is scalar leadfield (Option A) — the validation version, not the final faithful CARSS.** Once the whole pipeline works scalar, we upgrade to the 3-orientation subspace version. The eigenvalue split itself (Step 1) is identical either way, so this step is valid regardless.

[Likely] **One dominant eigenvalue is expected here because your source is a single coherent patch with one time course.** If you'd used the *two-source* data, you'd expect *two* dominant eigenvalues. That's actually a nice secondary check later — the eigenvalue count should match the number of sources. For now, single source → one eigenvalue.

Run it, paste the output (the printed eigenvalues + the estimated n), and tell me what the plot shows. If it's one dominant value and n=1, we go straight to Step 2 — the localizer that scans all 29666 dipoles and should peak exactly at dipole 29467, your known source. That Step 2 peak landing on the true source is the *real* proof CARSS is working, and it's the next gate.