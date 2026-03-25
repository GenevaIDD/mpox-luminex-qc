#!/usr/bin/env python3
"""Generate macOS .icns and Windows .ico app icons from the logo."""

import platform
import subprocess
import tempfile
from pathlib import Path

from PIL import Image

# Paths
PROJECT_ROOT = Path(__file__).parent.parent
LOGO_SRC = PROJECT_ROOT / "gdd_antibody_square_tighter.png"
OUTPUT_ICNS = PROJECT_ROOT / "assets" / "app_icon.icns"
OUTPUT_ICO = PROJECT_ROOT / "assets" / "app_icon.ico"


def pad_to_square(img: Image.Image, size: int) -> Image.Image:
    """Pad image to a square canvas, centered, with transparent background."""
    # Add some padding around the logo (10% margin)
    margin = int(size * 0.05)
    inner_size = size - 2 * margin

    # Scale image to fit within inner_size, preserving aspect ratio
    img_copy = img.copy()
    img_copy.thumbnail((inner_size, inner_size), Image.LANCZOS)

    # Create square transparent canvas
    canvas = Image.new("RGBA", (size, size), (0, 0, 0, 0))

    # Center the image
    x = (size - img_copy.width) // 2
    y = (size - img_copy.height) // 2
    canvas.paste(img_copy, (x, y), img_copy if img_copy.mode == "RGBA" else None)

    return canvas


def main():
    OUTPUT_ICNS.parent.mkdir(parents=True, exist_ok=True)

    logo = Image.open(LOGO_SRC).convert("RGBA")
    print(f"Source logo: {LOGO_SRC} ({logo.width}x{logo.height})")

    # macOS iconset requires these sizes (name -> pixel size)
    icon_sizes = {
        "icon_16x16.png": 16,
        "icon_16x16@2x.png": 32,
        "icon_32x32.png": 32,
        "icon_32x32@2x.png": 64,
        "icon_128x128.png": 128,
        "icon_128x128@2x.png": 256,
        "icon_256x256.png": 256,
        "icon_256x256@2x.png": 512,
        "icon_512x512.png": 512,
        "icon_512x512@2x.png": 1024,
    }

    with tempfile.TemporaryDirectory() as tmpdir:
        iconset_dir = Path(tmpdir) / "app_icon.iconset"
        iconset_dir.mkdir()

        for name, size in icon_sizes.items():
            img = pad_to_square(logo, size)
            img.save(iconset_dir / name, "PNG")
            print(f"  Created {name} ({size}x{size})")

        # Convert to .icns using macOS iconutil (only available on macOS)
        if platform.system() == "Darwin":
            subprocess.run(
                ["iconutil", "-c", "icns", str(iconset_dir), "-o", str(OUTPUT_ICNS)],
                check=True,
            )
            print(f"\nmacOS icon saved to: {OUTPUT_ICNS}")
        else:
            print("\nSkipping .icns generation (not on macOS)")

    # Windows .ico (sizes: 16, 32, 48, 64, 128, 256)
    ico_sizes = [16, 32, 48, 64, 128, 256]
    ico_images = [pad_to_square(logo, s) for s in ico_sizes]
    ico_images[0].save(
        OUTPUT_ICO,
        format="ICO",
        sizes=[(s, s) for s in ico_sizes],
        append_images=ico_images[1:],
    )
    print(f"Windows icon saved to: {OUTPUT_ICO}")


if __name__ == "__main__":
    main()
