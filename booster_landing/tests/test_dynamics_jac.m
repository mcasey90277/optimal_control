% TEST_DYNAMICS_JAC  Complex-step vs analytic Jacobians of pdg_dynamics.
%
% Perturbs each of the 10 inputs (7 state + 3 thrust) with h=1e-30i and
% compares imag(f)/h against the analytic A and B columns. Run twice:
% vacuum (Phase 1 default) and with the drag branch forced on, at a state
% with nonzero velocity (the ||v||v drag term is not differentiable at
% v=0; we never linearize there).
%
% INPUTS: none   OUTPUTS: none (throws on failure)
here_ = fileparts(mfilename('fullpath'));
addpath(fullfile(here_, '..'));  setup_paths;

P  = booster_params();
x0 = [400; 80; 1500; -25; 5; -140; 28000];
T0 = [2e4; -1e4; 5e5];
h  = 1e-30;

for dragOn = [false true]
    P.drag.on = dragOn;
    [~, A, B] = pdg_dynamics(x0, T0, P);
    Acs = zeros(7,7);  Bcs = zeros(7,3);
    for k = 1:7
        xp = complex(x0);  xp(k) = xp(k) + 1i*h;
        Acs(:,k) = imag(pdg_dynamics(xp, complex(T0), P)) / h;
    end
    for k = 1:3
        Tp = complex(T0);  Tp(k) = Tp(k) + 1i*h;
        Bcs(:,k) = imag(pdg_dynamics(complex(x0), Tp, P)) / h;
    end
    errA = max(abs(A(:) - Acs(:))) / max(1, max(abs(Acs(:))));
    errB = max(abs(B(:) - Bcs(:))) / max(1, max(abs(Bcs(:))));
    assert(errA < 1e-12, 'A mismatch (drag=%d): %.2e', dragOn, errA);
    assert(errB < 1e-12, 'B mismatch (drag=%d): %.2e', dragOn, errB);
end
fprintf('test_dynamics_jac PASS (vacuum + drag)\n');
