#!/usr/bin/env python3
"""Create an interactive 3D kfPoint/eigenvector scatter plot from a gap JSON file.

The file format used by Gap_test.jsn stores eigenvectors as rows. The script
defaults to the eigenvector whose eigenvalue is the second closest to 1, which
matches the current plot. Use --rank 1 for the closest eigenvalue instead.
"""

from __future__ import annotations

import argparse
import html
import json
import math
import webbrowser
from collections.abc import Iterable
from pathlib import Path
from typing import Any

PLOTLY_CDN = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"


def flatten_numbers(value: Any) -> Iterable[float]:
    """Yield numbers from a possibly nested JSON array."""
    if isinstance(value, list):
        for item in value:
            yield from flatten_numbers(item)
    else:
        yield float(value)


def select_eigenpair(eigenvalues: list[float], rank: int) -> tuple[int, float]:
    """Return (zero-based index, value) for the requested closeness rank."""
    if not eigenvalues:
        raise ValueError("Eigenvalues is empty")
    if rank < 1 or rank > len(eigenvalues):
        raise ValueError(f"--rank must be between 1 and {len(eigenvalues)}")

    ordered = sorted(
        enumerate(eigenvalues),
        key=lambda pair: (abs(pair[1] - 1.0), pair[0]),
    )
    return ordered[rank - 1]


def extract_eigenvector(
    matrix: Any,
    eigenpair_index: int,
    point_count: int,
    orientation: str,
) -> list[float]:
    """Extract one eigenvector from a row- or column-oriented matrix."""
    if not isinstance(matrix, list) or not matrix:
        raise ValueError("Eigenvectors must be a non-empty 2D array")

    if orientation == "rows":
        if eigenpair_index >= len(matrix):
            raise ValueError("Eigenvector row is missing")
        vector = [float(value) for value in matrix[eigenpair_index]]
    else:
        try:
            vector = [float(row[eigenpair_index]) for row in matrix]
        except (IndexError, TypeError) as exc:
            raise ValueError("Eigenvector column is missing") from exc

    if len(vector) != point_count:
        raise ValueError(
            f"Selected eigenvector has {len(vector)} values, "
            f"but kfPoints has {point_count} points"
        )
    if not all(math.isfinite(value) for value in vector):
        raise ValueError("Selected eigenvector contains a non-finite value")
    return vector


def load_plot_data(
    input_path: Path,
    rank: int,
    orientation: str,
) -> tuple[list[list[float]], int, float]:
    """Load and validate points, then select the requested eigenvector."""
    with input_path.open("r", encoding="utf-8") as handle:
        data = json.load(handle)

    try:
        raw_points = data["kfPoints"]
        raw_eigenvalues = data["Eigenvalues"]
        raw_eigenvectors = data["Eigenvectors"]
    except KeyError as exc:
        raise ValueError(f"Input is missing required field: {exc.args[0]}") from exc

    points: list[list[float]] = []
    for index, point in enumerate(raw_points):
        if not isinstance(point, list) or len(point) != 3:
            raise ValueError(f"kfPoints[{index}] is not a three-component point")
        xyz = [float(component) for component in point]
        if not all(math.isfinite(component) for component in xyz):
            raise ValueError(f"kfPoints[{index}] contains a non-finite value")
        points.append(xyz)

    eigenvalues = list(flatten_numbers(raw_eigenvalues))
    eigenpair_index, eigenvalue = select_eigenpair(eigenvalues, rank)
    eigenvector = extract_eigenvector(
        raw_eigenvectors,
        eigenpair_index,
        len(points),
        orientation,
    )

    rows = [
        [round(x, 8), round(y, 8), round(z, 8), round(value, 12)]
        for (x, y, z), value in zip(points, eigenvector)
    ]
    return rows, eigenpair_index, eigenvalue


def make_html(
    rows: list[list[float]],
    eigenpair_index: int,
    eigenvalue: float,
    rank: int,
    marker_size: float,
    title: str,
) -> str:
    """Build a standalone, responsive Plotly document."""
    amplitudes = [row[3] for row in rows]
    max_abs = max((abs(value) for value in amplitudes), default=1.0) or 1.0
    compact_rows = json.dumps(rows, separators=(",", ":"))
    safe_title = html.escape(title)
    eigenvalue_text = f"{eigenvalue:.10g}"

    return f"""<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>{safe_title}</title>
  <script src="{PLOTLY_CDN}"></script>
  <style>
    :root {{ color-scheme: light dark; }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      padding: 20px;
      background: light-dark(#ffffff, #111827);
      color: light-dark(#172033, #e5e7eb);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }}
    header {{
      display: flex;
      align-items: baseline;
      justify-content: space-between;
      gap: 12px;
      flex-wrap: wrap;
      margin-bottom: 8px;
    }}
    h1 {{ margin: 0; font-size: 1.25rem; font-weight: 600; }}
    .meta {{ color: light-dark(#5b6472, #aeb7c5); font-variant-numeric: tabular-nums; }}
    #plot {{ width: 100%; height: min(82vh, 850px); min-height: 520px; }}
  </style>
</head>
<body>
  <header>
    <h1>{safe_title}</h1>
    <div class="meta">Eigenpair {eigenpair_index + 1} · λ = {eigenvalue_text} · closeness rank {rank} · {len(rows):,} points</div>
  </header>
  <div id="plot" role="img" aria-label="Interactive three-dimensional k-space scatter plot colored by eigenvector amplitude."></div>

  <script>
    const rows = {compact_rows};
    const maxAbs = {max_abs:.16g};
    const plot = document.getElementById("plot");
    const darkQuery = window.matchMedia("(prefers-color-scheme: dark)");

    function draw() {{
      const dark = darkQuery.matches;
      const foreground = dark ? "#e5e7eb" : "#172033";
      const muted = dark ? "#aeb7c5" : "#5b6472";
      const border = dark ? "#4b5563" : "#c9ced6";
      const background = dark ? "#111827" : "#ffffff";
      const colorscale = [
        [0.0, "#2166ac"],
        [0.5, background],
        [1.0, "#b2182b"]
      ];
      const camera = plot.layout?.scene?.camera || {{eye: {{x: 1.45, y: 1.45, z: 1.15}}}};
      const axis = name => ({{
        title: {{text: name, font: {{color: foreground, size: 12}}}},
        color: foreground,
        tickfont: {{color: muted, size: 11}},
        gridcolor: border,
        zerolinecolor: border,
        showbackground: false,
        showspikes: false
      }});

      Plotly.react(plot, [{{
        type: "scatter3d",
        mode: "markers",
        x: rows.map(row => row[0]),
        y: rows.map(row => row[1]),
        z: rows.map(row => row[2]),
        customdata: rows.map((row, index) => [index + 1, row[3]]),
        marker: {{
          size: {marker_size:.8g},
          opacity: 0.9,
          color: rows.map(row => row[3]),
          cmin: -maxAbs,
          cmax: maxAbs,
          cmid: 0,
          colorscale,
          colorbar: {{
            title: {{text: "Eigenvector<br>amplitude", side: "right"}},
            tickformat: ".3f",
            thickness: 14,
            len: 0.72,
            outlinecolor: border,
            outlinewidth: 1
          }}
        }},
        hovertemplate:
          "<b>Point %{{customdata[0]}}</b><br>" +
          "kₓ = %{{x:.5f}}<br>" +
          "kᵧ = %{{y:.5f}}<br>" +
          "k_z = %{{z:.5f}}<br>" +
          "amplitude = %{{customdata[1]:.6f}}<extra></extra>"
      }}], {{
        autosize: true,
        margin: {{l: 0, r: 18, t: 6, b: 0}},
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: {{color: foreground, size: 12}},
        hoverlabel: {{bgcolor: background, bordercolor: border, font: {{color: foreground}}}},
        scene: {{
          xaxis: axis("kₓ"),
          yaxis: axis("kᵧ"),
          zaxis: axis("k_z"),
          aspectmode: "data",
          camera
        }},
        showlegend: false,
        uirevision: "gap-eigenvector"
      }}, {{
        responsive: true,
        displaylogo: false,
        scrollZoom: true,
        modeBarButtonsToRemove: ["sendDataToCloud", "lasso2d", "select2d"]
      }});
    }}

    draw();
    darkQuery.addEventListener("change", draw);
  </script>
</body>
</html>
"""


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Plot kfPoints in 3D, colored by a selected eigenvector."
    )
    parser.add_argument("input", type=Path, help="Input .jsn or JSON file")
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        help="Output HTML path (default: beside the input file)",
    )
    parser.add_argument(
        "--rank",
        type=int,
        default=1,
        help="1-based rank by distance |eigenvalue - 1| (default: 1)",
    )
    parser.add_argument(
        "--eigenvectors-as",
        choices=("rows", "columns"),
        default="rows",
        help="How eigenvectors are stored in the matrix (default: rows)",
    )
    parser.add_argument(
        "--marker-size",
        type=float,
        default=4.0,
        help="3D marker size (default: 4)",
    )
    parser.add_argument(
        "--title",
        default="Gap eigenvector on k-space points",
        help="Plot title",
    )
    parser.add_argument(
        "--open",
        action="store_true",
        help="Open the generated plot in the default browser",
    )
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    input_path = args.input.expanduser().resolve()
    if not input_path.is_file():
        raise SystemExit(f"Input file not found: {input_path}")
    if args.marker_size <= 0:
        raise SystemExit("--marker-size must be positive")

    output_path = (
        args.output.expanduser().resolve()
        if args.output
        else input_path.with_name(
            f"{input_path.stem}_eigenvector_rank{args.rank}_3d.html"
        )
    )

    try:
        rows, eigenpair_index, eigenvalue = load_plot_data(
            input_path,
            args.rank,
            args.eigenvectors_as,
        )
        document = make_html(
            rows,
            eigenpair_index,
            eigenvalue,
            args.rank,
            args.marker_size,
            args.title,
        )
    except (OSError, json.JSONDecodeError, TypeError, ValueError) as exc:
        raise SystemExit(f"Could not create plot: {exc}") from exc

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(document, encoding="utf-8")

    print(f"Selected eigenpair: {eigenpair_index + 1}")
    print(f"Eigenvalue: {eigenvalue:.13g}")
    print(f"Distance from 1: {abs(eigenvalue - 1.0):.13g}")
    print(f"Wrote: {output_path}")

    if args.open:
        webbrowser.open(output_path.as_uri())


if __name__ == "__main__":
    main()
