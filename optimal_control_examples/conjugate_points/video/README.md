# Making the conjugate-point video

The video is **generated**, not screen-recorded. A MATLAB "director" drives the
explorer headless through a shot list, timed to the narration audio, and writes
a 1920×1080 MP4. You supply the voice; everything else is scripted and
repeatable.

```
narration/sceneNN.txt   what the narrator says, one file per scene (written for speech)
video_scenes.m          the shot list: which view, caption, preset at which phrase
make_draft_voice.sh     a timing-draft voice from the Mac's built-in speech
make_conjugate_video.m  the director: renders video + narration track, joins them
mux_video.sh            re-join video and audio without re-rendering
audio/                  sceneNN.wav (draft) or sceneNN.mp3 (final voice)
out/                    the results (not in git)
```

## Step by step

### 1. Draft voice (1 minute)
```
cd optimal_control_examples/conjugate_points/video
./make_draft_voice.sh
```
Writes `audio/scene01.wav … scene07.wav` with the Mac's "Samantha" voice. It's
robotic, but it fixes the timing so you can check the whole video.

### 2. Render and watch the draft (about 10 min)
In MATLAB:
```matlab
cd optimal_control_examples/conjugate_points/video
make_conjugate_video                          % full quality
make_conjugate_video(struct('preview', true)) % quick 12-fps look (~3 min)
make_conjugate_video(struct('scenes', 3))     % one scene only
```
Open `out/conjugate_points_video.mp4` in QuickTime. If ffmpeg is not installed,
you get `out/conjugate_points_silent.mp4` + `out/narration.wav`: install ffmpeg
(`brew install ffmpeg`) and run `./mux_video.sh`, or join the two in iMovie.

### 3. The real voice (ElevenLabs)
1. Make a free account at elevenlabs.io (about 10,000 characters a month free;
   the whole narration is about 4,000).
2. Pick a voice. Calm, clear "narrator" voices suit this; try a few on one scene.
3. For each `narration/sceneNN.txt`: paste the text, generate, and download the
   MP3. Save it as `audio/sceneNN.mp3`, keeping the same number.
4. Re-run `make_conjugate_video`. It prefers `.mp3` over the draft `.wav`, and
   the pictures re-time themselves to the new voice automatically.

A scene that sounds wrong can be regenerated on its own: replace that one MP3
and re-render.

### 4. Upload
1. studio.youtube.com → **Create → Upload video** → choose
   `out/conjugate_points_video.mp4`.
2. Title, description and chapters: copy them from `../doc/youtube_script.md`.
3. Visibility **Unlisted** first; watch it once on YouTube; then Public.
4. Subtitles → upload the narration text (the scene files, concatenated) →
   **Auto-sync**.
5. For a custom thumbnail, verify the channel by phone once. A good thumbnail is
   the oscillator's fan crossing at π with the words "Conjugate Points".
6. Under "Altered content", tick yes (AI narration).

## Changing the video

- **Wording:** edit `narration/sceneNN.txt`, regenerate that scene's audio,
  and re-render. If you change a sentence that a beat is anchored to, update
  the anchor in `video_scenes.m`; the director stops with a clear error if an
  anchor phrase is missing.
- **Pictures:** edit `video_scenes.m`. Each beat is `anchor phrase → list of
  operations` (preset, extremal, delta, view, caption, card, shrink). The
  header lists every operation.
- **Timing:** beats fire where their phrase falls in the text, by character
  position, scaled to the audio length. That's accurate to about a second. To
  nudge a beat, anchor it a few words earlier or later.
- Beat times of the last render are in `out/beat_times.csv`.
