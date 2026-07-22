Here is the implementation audit for the Phase-1 codebase with exact paths, defects, and concrete fixes.

### **1. CACHE/RESUME/FINGERPRINT SEMANTICS (S1)**
**[CORRECTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:225`
**Issue:** `check_fp_local` silently ignores fields in the user's current fingerprint if they are missing in the cached `S.fp` (schema-older path). This allows a stale cache built before a new field (like `x0Elems` or `gain`) was introduced to be silently accepted for a completely different run.
**Fix:** Change the loop to explicitly hard-error on missing fields:
```matlab
    if ~isfield(S.fp, f)
        error('run_cr3bp_geo:fpSchemaOlder', 'checkpoint %s missing field ''%s'' -- stale cache config, change runName or rerun.', file, f);
    elseif ~isequal(S.fp.(f), fp.(f))
        error('run_cr3bp_geo:fpMismatch', ...);
    end
```

**[CORRECTNESS]** `orbit_transfer/earth_elliptic_to_geo/direct/core/homotopy_mee.m:132` (and identical blocks in `bridge_mu_continuation.m:335`, `solve_cr3bp_minfuel.m:233`)
**Issue:** The schema-older path explicitly uses `warning(...)` and `continue;` allowing potentially fatal tag collisions if configs drift between script edits.
**Fix:** Upgrade the `warning(...); continue;` statements in all three files to `error(...)` to guarantee configuration hygiene.

**[CORRECTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:72`
**Issue:** The struct `fp` omits the user's `x0Elems`, `ctf`, and `maxIter`. If the user updates the start state or node density but doesn't bump the `.mat` tag names, the cache is still matched.
**Fix:** Embed missing components so parameter edits forcibly break existing caches:
```matlab
if isempty(x0Elems), x0_fp = []; else, x0_fp = x0Elems(:).'; end
fp = struct(..., 'ctf',ctf, 'maxIter',maxIter, 'x0Elems',x0_fp);
```

**[CORRECTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:98` and `115`
**Issue:** Partial cache deletes blindly assume dimension safety. If `ckSeed` is missing but `ckE`/`ckG` are loaded, `mee_seed` will silently regenerate `sigma`, potentially at a different `N`. The mismatch (`len(sigma)` vs `len(S.o.X)`) crashes CasADi later.
**Fix:** Add `N` mesh size safety checks right after `S = load(ckE)` and `S = load(ckG)`:
```matlab
assert(size(o.X, 2) == numel(sigma)+1, 'run_cr3bp_geo:meshMismatch', 'Partial delete caused spatial mismatch (N=%d vs %d). Rerun required.', size(o.X,2), numel(sigma)+1);
```

### **2. BASIN CONTROL (S2)**
**[CORRECTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:88`
**Issue:** `seedOpts.N = max(60, round(nodesPerRev * infoP.nRev))` silently resizes the optimal grid density set by the campaign recipe, fundamentally biasing which basin is acquired. It also misses the certified $6.5 \leq nRev \leq 9$ window check altogether.
**Fix:** Use exact recipe logic and enforce the window:
```matlab
assert(infoP.nRev >= 6.5 && infoP.nRev <= 9, 'run_cr3bp_geo:seedWindow', 'Seed nRev (%.2f) outside [6.5, 9] window.', infoP.nRev);
seedOpts.N = round(nodesPerRev * infoP.nRev);
```

**[INTERFACE]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:144`
**Issue:** If the user specifies the standard Phase-1 config but hits a local minimum instead of the certified campaign basin (e.g., losing 1.84 kg), the script remains silent.
**Fix:** After getting `best` from the pipeline, add:
```matlab
if gain == 1 && epsMin == 0 && exist('cert','var') && ~isempty(cert)
    dmf = cert.m_f_kg - par.m0kg*best.mf;
    if dmf > 0.5   % Alert if we fall materially short
        warning('run_cr3bp_geo:suboptimal', 'Result materially worse than optimal baseline (-%.2f kg).', dmf);
    end
end
```

### **3. FRONT-DOOR ROBUSTNESS (S3)**
**[ROBUSTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:50`
**Issue:** Zero parameter validation at the beginning of the monolithic pipeline. Illegal lunar gains, negative epsilons, or mis-shaped matrices crash brutally deep inside CasADi iteration logs.
**Fix:** Add strict bounds checks upfront:
```matlab
assert(isscalar(epsMin) && epsMin >= 0 && epsMin <= 1, 'run_cr3bp_geo:invalidParam', 'epsMin must be in [0,1]');
assert(isscalar(gain) && gain >= 0 && gain <= 1, 'run_cr3bp_geo:invalidParam', 'gain must be in [0,1]');
assert(isempty(x0Elems) || numel(x0Elems) == 5, 'run_cr3bp_geo:invalidParam', 'x0Elems must be 5x1');
assert(isempty(xfElems) || numel(xfElems) == 5, 'run_cr3bp_geo:invalidParam', 'xfElems must be 5x1');
```

**[ROBUSTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:67`
**Issue:** A blank `catch` wrapper hides non-trivial bugs (like a missing path dependency or `.m` syntax errors inside the table file) while safely defaulting non-table thrust rungs.
**Fix:** Filter by known exceptions and throw everything else:
```matlab
catch ME
    if ~strcmp(ME.identifier, 'MATLAB:UndefinedFunction') && ~contains(ME.message, 'not found')
        rethrow(ME);
    end
end
```

### **4. PERT BRANCH AND ORACLES (S4)**
**[CORRECTNESS]** `orbit_transfer/earth_elliptic_to_geo/direct/core/lt_mee_rhs.m:86`
**Issue:** The `$d_3$` distance expression is artificially padded: `d3 = (dx^2 + dy^2 + dz^2 + 1e-12)^1.5;`. By design this formulation holds physical bodies ($d \geq 8$ LU separation guaranteed). CasADi `pow` cleanly derivatives $(>0)^{1.5}$. Because of the `+1e-12`, CasADi's resultant equations carry a systemic bias distinct from the pure dynamics implemented inside `test_lt_mee_rhs_pert.m`.
**Fix:** Simplify the guard to exact mathematical equivalency without float conditioning: 
```matlab
    d3 = (dx^2 + dy^2 + dz^2)^1.5;
```

### **5. DATA PRODUCTS (S5)**
**[STYLE]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:18`
**Issue:** The meaning and domain orientation of `lamDef` are wholly undocumented in the front door.
**Fix:** Update the description block to clarify what constraints the multipliers bind to:
```matlab
%                              U, discrete costates lamDef (IPOPT multipliers for dynamics defects, 7x(N-1)), sigma/L mesh, dL,
```

**[ROBUSTNESS]** `orbit_transfer/earth_elliptic_to_geo_CR3BP/direct/run_cr3bp_geo.m:169`
**Issue:** Blind Cartesian mapping extracts `tD = cart.fuel.X(8,:)...` but time indexing natively in MEE arrays behaves as row `7`. While unplotted, executing `X(8,:)` against a raw mapped target risks a sudden indexing trap if internal state logic omits time.
**Fix:** Remove unused `tD` altogether. The plot section natively uses exact MEE state parameters (`par.m0kg*X(6,:)` and `tTU_s`); `tD` serves no purpose other than potential fault. Replace the mapping string with:
```matlab
r = cart.fuel.X(1:3,:);  sTh = cart.fuel.U(4,:); % tD extraction removed
```

---

### **OVERALL VERDICT:** 
**SHIP WITH FIXES**

**Top-3 Priority List:**
1. **Schema-older Cache Bypass (S1):** Change the silent-drop cache loading functions in `check_cache_fp`/`check_fp_local` to HARD ERRORS. Right now doing a `runName` collision while editing a fundamental setting (`maxIter`, `ctf`) seamlessly uses foreign dynamics.
2. **Deterministic Mesh Logic (S2/S1):** Remove the uncertified `max(60,...)` limiter in `run_cr3bp_geo.m` forcing deterministic density alignment alongside injecting the core safety `infoP.nRev >= 6.5` window checks. (Fixes silent mis-allocation to bad basins). 
3. **True Third-Body Distances (S4):** Remove the `1e-12` perturbation guard from `d^3` in `lt_mee_rhs.m` to eliminate numeric drift against the uncoupled `test_cr3bp_consistency.m` oracle limits.
