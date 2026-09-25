#!/bin/zsh
# make_draft_voice.sh -- a TIMING DRAFT of the narration with the Mac's built-in
# text-to-speech, one WAV per scene, into audio/.
#
#   ./make_draft_voice.sh                 # voice Samantha, 165 words/min
#   VOICE=Daniel RATE=160 ./make_draft_voice.sh
#
# Better free voices: System Settings > Accessibility > Spoken Content >
# System Voice > Manage Voices, download an "Enhanced"/"Premium" voice
# (e.g. "Ava (Premium)"), then VOICE="Ava (Premium)" ./make_draft_voice.sh
#
# The FINAL voice: put ElevenLabs (or any) files in audio/ as scene01.mp3 ...
# scene07.mp3. The director prefers .mp3 over .wav, so no code changes.
set -e
cd "$(dirname "$0")"
VOICE=${VOICE:-Samantha}
RATE=${RATE:-165}
mkdir -p audio
for f in narration/scene*.txt; do
    b=$(basename "$f" .txt)
    say -v "$VOICE" -r "$RATE" --data-format=LEI16@44100 -o "audio/$b.wav" -f "$f"
    echo "wrote audio/$b.wav"
done
