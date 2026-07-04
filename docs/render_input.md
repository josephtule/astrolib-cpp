# Render Input

This page documents the current first-pass Raylib render-loop keyboard and mouse input.

## Simulation

These inputs use no modifier key.

| Input | Action |
| --- | --- |
| `Space` | Pause/resume simulation stepping. |
| `Up` | Increase `WorldStepperConfig::ticks`. |
| `Down` | Decrease `WorldStepperConfig::ticks`, clamped to at least `1`. |
| `Right` | Increase `WorldStepperConfig::substeps`. |
| `Left` | Decrease `WorldStepperConfig::substeps`, clamped to at least `1`. |

## Render Toggles

Hold `Left Shift` for render/display toggles.

| Input | Action |
| --- | --- |
| `Left Shift + ,` | Toggle body axes. |
| `Left Shift + .` | Toggle inertial axes. |
| `Left Shift + /` | Toggle colored axes. |
| `Left Shift + 1` | Toggle XY grid. |
| `Left Shift + 2` | Toggle XZ grid. |
| `Left Shift + 3` | Toggle ZY grid. |
| `Left Shift + `` | Toggle FPS overlay. |
| `Left Shift + Backspace` | Reset render draw toggles to defaults. |

## Camera Settings

Hold `Right Shift` for camera settings.

| Input | Action |
| --- | --- |
| `Right Shift + 1` | Set camera mode to locked. |
| `Right Shift + 2` | Set camera mode to target. |
| `Right Shift + 3` | Set camera mode to origin. |
| `Right Shift + 4` | Set camera mode to free. |
| `Right Shift + mouse wheel` | Adjust camera zoom/FOV step rate. |
| `Right Shift + `` | Toggle mouse-wheel inversion. |
| `Right Shift + -` | Decrease camera FOV. |
| `Right Shift + =` | Increase camera FOV. |
| `Right Shift + \` | Reset camera FOV and zoom rate. |
| `Right Shift + Backspace` | Reset camera config to defaults. |

## Camera Speed Settings

Hold `Right Shift` and use arrow keys to adjust the active camera-speed setting.

In free mode:

| Input | Action |
| --- | --- |
| `Right Shift + Up` | Increase free-camera fly speed by `1000`. |
| `Right Shift + Right` | Increase free-camera fly speed by `100`. |
| `Right Shift + Down` | Decrease free-camera fly speed by `1000`, clamped to at least `100`. |
| `Right Shift + Left` | Decrease free-camera fly speed by `100`, clamped to at least `100`. |

In origin or target mode:

| Input | Action |
| --- | --- |
| `Right Shift + Up` | Increase orbit speed by `1.0`. |
| `Right Shift + Right` | Increase orbit speed by `0.1`. |
| `Right Shift + Down` | Decrease orbit speed by `1.0`, clamped to at least `0.1`. |
| `Right Shift + Left` | Decrease orbit speed by `0.1`, clamped to at least `0.1`. |

## Camera Movement

Movement uses no modifier key and depends on the current camera mode.

### Origin/Target Mode

| Input | Action |
| --- | --- |
| `W` / `S` | Orbit up/down around the pivot. |
| `A` / `D` | Orbit left/right around the pivot. |
| `Left mouse drag` | Orbit around the pivot. |

In origin mode, the pivot is the world origin. In target mode, the pivot is the target entity when available.

### Free Mode

| Input | Action |
| --- | --- |
| `W` / `S` | Move forward/back along the camera look direction. |
| `A` / `D` | Strafe left/right in the camera plane. |
| `Q` / `E` | Move down/up along the configured up axis. |
| `Left mouse drag` | Rotate camera look direction. |

### Locked Mode

Locked mode does not apply camera movement.

## Planned Input Work

- Replace hardcoded input checks with a key map.
- Add a key-map settings page after the UI shell exists.
- Move mouse sensitivity and scroll sensitivity into runtime render settings.
