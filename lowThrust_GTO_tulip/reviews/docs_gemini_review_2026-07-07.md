1. **`gto_tulip_mintime_theory.tex`, Section 4 ("The direct route: transcription to an NLP")**
   - **Severity**: Critical
   - **Issue**: The theory document defined the NLP decision vector using a 4-dimensional control per node (3D direction + throttle $s$), claiming `nZ = 33,012`. However, the tutorial correctly instructs the learner to fix the throttle at $1$ to prevent an interior-point solver stall, dropping the decision vector to 3-dimensional ($W \in \R^{3\times (N+1)}$) with `nZ = 30,011`. This direct contradiction would derail a learner implementing transcription sizes.
   - **Fix**: Redefined the control in the theory document as simply the unit vector $W \in \R^{3\times (N+1)}$ and constrained it via $\mathbf{w}_k\T\mathbf{w}_k = 1$. The NLP scale descriptions were properly synchronized to `$nZ = 10(N{+}1)+1 = 30{,}011$`. Also updated the $\partial d_k/\partial \mathbf{w}_{(\cdot)}$ Jacobian notation.

2. **`gto_tulip_mintime_theory.tex`, Section 4 ("The direct route: transcription to an NLP")**
   - **Severity**: Major
   - **Issue**: The pedagogical "ballast exploit" Pitfall box warned against relaxing $\mathbf{w}_k\T\mathbf{w}_k \le s_k^2$ based on the (now removed) independent $s_k$ throttle variable.
   - **Fix**: Rewrote the pitfall accurately for the fixed-throttle formulation. It now warns against relaxing the spherical bound to $\mathbf{w}_k\T\mathbf{w}_k \le 1$. Emphasized that since fixing throttle implies maximum uninterrupted mass flow unconditionally, lowering $\|\mathbf{w}\|$ allows the optimizer to exploit the physics by dumping propellant without thrusting to achieve unphysical terminal acceleration.

3. **`building_the_gto_tulip_solvers.tex`, Section 2 ("Phase A: the augmented PMP dynamics")**
   - **Severity**: Minor
   - **Issue**: The code skeleton explicitly requested the output signature `function [yDot, Ht, S, aThrust] = lt_pmp_eom(...)`, but `aThrust` was never defined in the instructions, derived, or used in the checkpoint. A learner would be confused about what variable mapping satisfies this signature.
   - **Fix**: Removed `aThrust` from the boilerplate signature, mapping it elegantly to the variables actually described and verified in the instructions: `[yDot, Ht, S]`.


***

### The 3 Highest-Priority Edits Made

1. **Synchronizing Transcription Size:** Replaced the length-4 control arrays $U$ with length-3 direction controls $W$ in the theory note and corrected mathematical claims on Jacobian block sizes from $33,012$ string to $30{,}011$ to match Phase D of the coding tutorial.
2. **Rewriting the Ballast Exploit:** Rephrased the physics justification of the pitfall in the theory note to establish *why* replacing equality with $\mathbf{w}\T\mathbf{w} \le 1$ causes mass-dumping under a fixed throttle NLP model. 
3. **Removing Extraneous Signature Outputs:** Stripped the undefined `aThrust` output from the Phase A `lt_pmp_eom` skeleton in the tutorial to unblock learners. 

*(Math, Physics, and Checkpoints—including complex sensitivity claims, costate gradient norms, transversality consistency logic with $\lambda_m > 0$ driving down to $0$, and relative variable conversions—were rigorously verified to be highly consistent and accurate!)*
