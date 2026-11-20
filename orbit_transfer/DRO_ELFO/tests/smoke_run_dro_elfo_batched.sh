#!/bin/bash
# SMOKE_RUN_DRO_ELFO_BATCHED  Smoke test of catalog/run_dro_elfo_batched.sh
# with a stub MATLAB (MATLAB_BIN) in a mktemp -d folder; no MATLAB runs.
#
# Cases: done-marker stop; a foreign / mismatched marker is ignored;
# no-progress stop; MATLAB error abort; hang kill (OS watchdog); bad
# arguments refused before MATLAB; lock refusal (live and stale holder, the
# holder's lock left intact); SIGTERM to the driver (stub child gone, lock
# removed). Prints PASS/FAIL per case, cleans up, exits nonzero on any FAIL.
#
# Usage: ./smoke_run_dro_elfo_batched.sh        (about 1 minute)
#
# M. Casey                                                   (c) 11/19/2026
# Copyright Coorbital Inc.

set -u
HERE="$(cd "$(dirname "$0")" && pwd)"
DRIVER="${SMOKE_DRIVER:-$HERE/../catalog/run_dro_elfo_batched.sh}"   # override: mutation checks
TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT
NFAIL=0

pass() { echo "PASS  $1"; }
fail() { echo "FAIL  $1 -- $2"; NFAIL=$((NFAIL+1)); }

# --- the stub: behaviour by STUB_MODE; counts calls; records its pid -------
cat > "$TMP/stub.sh" <<'EOF'
#!/bin/bash
CMD="$2"
RID=$(printf '%s' "$CMD" | sed -n "s/.*'runId', '\([^']*\)'.*/\1/p")
N=$(cat "$STUB_DIR/count" 2>/dev/null || echo 0); N=$((N+1)); echo "$N" > "$STUB_DIR/count"
LOG="$OUTDIR/combined_progress.txt"
PROG="$OUTDIR/dro_elfo_tau1_raan000_progress.txt"
case "$STUB_MODE" in
  complete)
    echo "SHEETS SELECTED [run $RID]: dro_elfo_tau1_raan000" >> "$LOG"
    echo "line $N" >> "$PROG"
    [ "$N" -ge 2 ] && echo "SELECTED SHEETS COMPLETE [run $RID]: dro_elfo_tau1_raan000" >> "$LOG"
    exit 0;;
  foreign)
    echo "SHEETS SELECTED [run $RID]: dro_elfo_tau1_raan000" >> "$LOG"
    echo "SELECTED SHEETS COMPLETE [run 1-1]: dro_elfo_tau1_raan000" >> "$LOG"
    echo "SELECTED SHEETS COMPLETE [run $RID]: dro_tulip_tau1" >> "$LOG"
    [ "$N" -eq 1 ] && echo "line $N" >> "$PROG"
    exit 0;;
  dry) exit 0;;
  fail) echo "Error using run_dro_elfo_torus"; exit 1;;
  hang) echo $$ > "$STUB_DIR/stub.pid"; exec sleep 300;;
esac
EOF
chmod +x "$TMP/stub.sh"

# run_case <name> <mode> [driver args...]: fresh OUTDIR, driver in foreground
run_case() {
  local name=$1 mode=$2; shift 2
  export STUB_DIR="$TMP/$name" OUTDIR="$TMP/$name/out" MATLAB_BIN="$TMP/stub.sh" STUB_MODE="$mode"
  mkdir -p "$OUTDIR"; rm -f "$STUB_DIR/count"
  "$DRIVER" "$@" > "$STUB_DIR/stdout.txt" 2>&1; RC=$?
  CALLS=$(cat "$STUB_DIR/count" 2>/dev/null || echo 0)
  CLOG="$OUTDIR/combined_progress.txt"
}

# 1 done-marker stop
run_case done complete 1 2 60
if [ $RC -eq 0 ] && [ "$CALLS" -eq 2 ] && grep -q "completion marker of run" "$CLOG" && [ ! -d "$OUTDIR/.dro_elfo_driver.lock" ]; then
  pass "done-marker stop (2 MATLAB calls, lock released)"
else fail "done-marker stop" "rc=$RC calls=$CALLS"; fi

# 2 foreign / mismatched marker ignored (stops by the no-progress rule instead)
run_case foreign foreign 1 2 60
if [ $RC -eq 0 ] && ! grep -q "completion marker" "$CLOG" && grep -q "did not grow for 3 passes" "$CLOG" && [ "$CALLS" -eq 4 ]; then
  pass "foreign marker ignored (4 calls, wedge stop)"
else fail "foreign marker ignored" "rc=$RC calls=$CALLS"; fi

# 3 no-progress stop
run_case dry dry 1 2 60
if [ $RC -eq 0 ] && [ "$CALLS" -eq 3 ] && grep -q "did not grow for 3 passes" "$CLOG"; then
  pass "no-progress stop (3 calls)"
else fail "no-progress stop" "rc=$RC calls=$CALLS"; fi

# 4 MATLAB error abort
run_case err fail 1 2 60
if [ $RC -eq 1 ] && [ "$CALLS" -eq 1 ] && grep -q "ABORTING" "$CLOG" && [ ! -d "$OUTDIR/.dro_elfo_driver.lock" ]; then
  pass "error abort (exit 1 after 1 call, lock released)"
else fail "error abort" "rc=$RC calls=$CALLS"; fi

# 5 hang kill by the OS watchdog (BSEC 0, KILL_MARGIN 2)
KILL_MARGIN=2 run_case hang hang 1 2 0
SPID=$(cat "$STUB_DIR/stub.pid" 2>/dev/null)
NKILL=$(grep -c "timed out" "$CLOG" 2>/dev/null)
if [ $RC -eq 0 ] && [ "$NKILL" -eq 3 ] && [ -n "$SPID" ] && ! kill -0 "$SPID" 2>/dev/null; then
  pass "hang kill (3 watchdog kills, stub gone)"
else fail "hang kill" "rc=$RC kills=$NKILL stubpid=$SPID"; fi

# 6 bad arguments refused before MATLAB
BADOK=1
for args in "1;quit|2|60" "1|x|60" "1|2|1e3" "1|2|60|[15 12];exit" "1|2|60||0.3;exit"; do
  IFS='|' read -r a1 a2 a3 a4 a5 <<< "$args"
  run_case bad dry "$a1" "$a2" "$a3" "${a4:-}" "${a5:-}"
  if [ $RC -ne 2 ] || [ "$CALLS" -ne 0 ]; then BADOK=0; echo "      bad args '$args': rc=$RC calls=$CALLS"; fi
done
KILL_MARGIN="9;x" run_case badk dry 1 2 60
if [ $RC -ne 2 ] || [ "$CALLS" -ne 0 ]; then BADOK=0; echo "      bad KILL_MARGIN: rc=$RC"; fi
if [ $BADOK -eq 1 ]; then pass "bad arguments refused (exit 2, MATLAB never called)"; else fail "bad arguments" "see above"; fi

# 7 lock refusal: a live holder (this shell) and a stale one (a dead pid)
export STUB_DIR="$TMP/lock" OUTDIR="$TMP/lock/out" MATLAB_BIN="$TMP/stub.sh" STUB_MODE=dry
mkdir -p "$OUTDIR/.dro_elfo_driver.lock"; rm -f "$STUB_DIR/count"
echo "$$ 2026-11-19 00:00:00" > "$OUTDIR/.dro_elfo_driver.lock/pid"
OUT1=$("$DRIVER" 1 2 60 2>&1); RC1=$?
( exit 0 ) & DEAD=$!; wait $DEAD
echo "$DEAD 2026-11-19 00:00:00" > "$OUTDIR/.dro_elfo_driver.lock/pid"
OUT2=$("$DRIVER" 1 2 60 2>&1); RC2=$?
if [ $RC1 -eq 1 ] && [ $RC2 -eq 1 ] && printf '%s' "$OUT1" | grep -q "pid $$ .* is ALIVE" && \
   printf '%s' "$OUT2" | grep -q "is NOT running (stale lock)" && printf '%s' "$OUT2" | grep -q "rm -rf" && \
   [ -f "$OUTDIR/.dro_elfo_driver.lock/pid" ] && [ ! -f "$STUB_DIR/count" ]; then
  pass "lock refusal (live and stale holder named, holder's lock intact)"
else fail "lock refusal" "rc=$RC1/$RC2 out1='$OUT1' out2='$OUT2'"; fi

# 8 SIGTERM to the driver: the MATLAB stub child is stopped, the lock removed
export STUB_DIR="$TMP/term" OUTDIR="$TMP/term/out" MATLAB_BIN="$TMP/stub.sh" STUB_MODE=hang
mkdir -p "$OUTDIR"; rm -f "$STUB_DIR/stub.pid"
"$DRIVER" 1 2 60 > "$STUB_DIR/stdout.txt" 2>&1 &
DPID=$!
for _ in $(seq 1 50); do [ -s "$STUB_DIR/stub.pid" ] && break; sleep 0.2; done
SPID=$(cat "$STUB_DIR/stub.pid" 2>/dev/null)
LOCKPID=$(cut -d' ' -f1 "$OUTDIR/.dro_elfo_driver.lock/pid" 2>/dev/null)
kill -TERM "$DPID"
wait "$DPID"; RC=$?
sleep 0.5
if [ -n "$SPID" ] && ! kill -0 "$SPID" 2>/dev/null && [ ! -d "$OUTDIR/.dro_elfo_driver.lock" ] && \
   [ "$LOCKPID" = "$DPID" ] && [ $RC -eq 143 ] && grep -q "received SIGTERM" "$OUTDIR/combined_progress.txt"; then
  pass "SIGTERM (exit 143, stub child gone, lock removed, pid file named the driver)"
else
  fail "SIGTERM" "rc=$RC stubpid=$SPID alive=$(kill -0 "$SPID" 2>/dev/null && echo yes || echo no) lockpid=$LOCKPID driver=$DPID"
  [ -n "$SPID" ] && kill -9 "$SPID" 2>/dev/null
fi

if [ $NFAIL -eq 0 ]; then echo "smoke_run_dro_elfo_batched: ALL PASS"; exit 0; fi
echo "smoke_run_dro_elfo_batched: $NFAIL FAILED"; exit 1
