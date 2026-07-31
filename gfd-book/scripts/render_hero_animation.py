#!/usr/bin/env python3
"""Generate the front-page hero animation from ch16's baroclinic life-cycle data.

Rerun (`python3 scripts/render_hero_animation.py` from the repo root) any
time the source data or visual parameters change -- this is the only
place site/static/media/hero-baroclinic.gif comes from. No ffmpeg
dependency: matplotlib writes the animated GIF directly through Pillow.
"""
import io

import matplotlib
import numpy as np
from PIL import Image

matplotlib.use("Agg")
import matplotlib.pyplot as plt

DATA = "notebooks/data/ch16_baroclinic-instability.npz"
OUT = "site/static/media/hero-baroclinic.gif"
START_FRAME = 38  # skip the near-zero-amplitude lead-in (tuned by eye)
FPS = 10
N_TILES = 4  # doubly-periodic field -> seamless horizontal tiling into a wide banner
HEIGHT_PIXELS = 200  # native raster height
WIDTH_PIXELS = HEIGHT_PIXELS * N_TILES  # must stay an exact N_TILES multiple of the height


def render():
    d = np.load(DATA)
    q1 = d["q1"][START_FRAME:]
    vmax = np.abs(q1).max()  # one global scale -- the growth itself is the story

    dpi = 100
    size_in = (WIDTH_PIXELS / dpi, HEIGHT_PIXELS / dpi)
    frames = []
    for field in q1:
        # tile axis 0: field.T (below) swaps axes for imshow, so tiling axis 0
        # here becomes the horizontal direction of the displayed banner
        tiled = np.tile(field, (N_TILES, 1))
        fig = plt.figure(figsize=size_in, dpi=dpi)
        ax = fig.add_axes((0, 0, 1, 1))
        ax.imshow(tiled.T, origin="lower", cmap="RdBu_r", vmin=-vmax, vmax=vmax, aspect="equal")
        ax.axis("off")
        buf = io.BytesIO()
        fig.savefig(buf, format="png", dpi=dpi)
        plt.close(fig)
        buf.seek(0)
        frames.append(Image.open(buf).convert("RGB"))

    # Ping-pong: forward, then back (excluding both endpoints) -- the only
    # cheap way to seamlessly loop a growth/breaking sequence whose start
    # and end states don't match.
    ping_pong = frames + frames[-2:0:-1]
    ping_pong[0].save(
        OUT,
        save_all=True,
        append_images=ping_pong[1:],
        duration=int(1000 / FPS),
        loop=0,
        optimize=True,
    )
    print(f"wrote {OUT}: {len(ping_pong)} frames, vmax={vmax:.2f}")


if __name__ == "__main__":
    render()
