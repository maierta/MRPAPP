#!/usr/bin/env julia

# Create an interactive 3D kfPoint/eigenvector scatter plot from a gap JSON file.
#
# The Gap_test.jsn format stores eigenvectors as rows. This script defaults to
# the eigenvector whose eigenvalue is the second closest to 1. Use --rank 1 for
# the closest eigenvalue instead.

using Printf
import JSON3

const PLOTLY_CDN = "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"


function usage(status::Int = 0)
    println(stderr, """
Usage:
  julia plot_gap_eigenvector.jl INPUT [options]

Options:
  -o, --output PATH              Output HTML path
  --rank N                       Rank by |eigenvalue - 1| (default: 2)
  --eigenvectors-as rows|columns Matrix orientation (default: rows)
  --marker-size NUMBER           Marker size (default: 8)
  --title TEXT                   Plot title
  --open                         Open the result in a browser
  -h, --help                     Show this help
""")
    exit(status)
end


function parse_args(args::Vector{String})
    input = nothing
    output = nothing
    rank = 2
    orientation = "rows"
    marker_size = 8.0
    title = "Gap eigenvector on k-space points"
    open_result = false

    i = 1
    while i <= length(args)
        arg = args[i]
        if arg in ("-h", "--help")
            usage(0)
        elseif arg in ("-o", "--output")
            i += 1
            i <= length(args) || error("$arg requires a value")
            output = args[i]
        elseif arg == "--rank"
            i += 1
            i <= length(args) || error("--rank requires a value")
            rank = tryparse(Int, args[i])
            rank === nothing && error("--rank must be an integer")
        elseif arg == "--eigenvectors-as"
            i += 1
            i <= length(args) || error("--eigenvectors-as requires a value")
            orientation = args[i]
        elseif arg == "--marker-size"
            i += 1
            i <= length(args) || error("--marker-size requires a value")
            marker_size = tryparse(Float64, args[i])
            marker_size === nothing && error("--marker-size must be a number")
        elseif arg == "--title"
            i += 1
            i <= length(args) || error("--title requires a value")
            title = args[i]
        elseif arg == "--open"
            open_result = true
        elseif startswith(arg, "-")
            error("Unknown option: $arg")
        elseif input === nothing
            input = arg
        else
            error("Unexpected positional argument: $arg")
        end
        i += 1
    end

    input === nothing && usage(1)
    rank >= 1 || error("--rank must be a positive integer")
    orientation in ("rows", "columns") ||
        error("--eigenvectors-as must be rows or columns")
    isfinite(marker_size) && marker_size > 0 ||
        error("--marker-size must be positive")

    return (
        input = String(input),
        output = output === nothing ? nothing : String(output),
        rank = rank,
        orientation = orientation,
        marker_size = marker_size,
        title = title,
        open_result = open_result,
    )
end


function flatten_numbers!(result::Vector{Float64}, value)
    if value isa AbstractArray
        for item in value
            flatten_numbers!(result, item)
        end
    else
        push!(result, Float64(value))
    end
    return result
end


function extract_plot_data(input_path::String, rank::Int, orientation::String)
    data = JSON3.read(read(input_path, String))
    haskey(data, :kfPoints) || error("Input is missing required field: kfPoints")
    haskey(data, :Eigenvalues) || error("Input is missing required field: Eigenvalues")
    haskey(data, :Eigenvectors) || error("Input is missing required field: Eigenvectors")

    raw_points = data[:kfPoints]
    point_count = length(raw_points)
    point_count > 0 || error("kfPoints is empty")

    points = Vector{NTuple{3, Float64}}(undef, point_count)
    for i in eachindex(raw_points)
        point = raw_points[i]
        length(point) == 3 || error("kfPoints[$i] is not a three-component point")
        xyz = (Float64(point[1]), Float64(point[2]), Float64(point[3]))
        all(isfinite, xyz) || error("kfPoints[$i] contains a non-finite value")
        points[i] = xyz
    end

    eigenvalues = flatten_numbers!(Float64[], data[:Eigenvalues])
    isempty(eigenvalues) && error("Eigenvalues is empty")
    rank <= length(eigenvalues) ||
        error("--rank cannot exceed the number of eigenvalues")

    ordering = sortperm(
        eachindex(eigenvalues),
        by = i -> (abs(eigenvalues[i] - 1.0), i),
    )
    eigenpair_index = ordering[rank]
    eigenvalue = eigenvalues[eigenpair_index]
    raw_eigenvectors = data[:Eigenvectors]

    eigenvector = if orientation == "rows"
        eigenpair_index <= length(raw_eigenvectors) || error("Eigenvector row is missing")
        Float64.(collect(raw_eigenvectors[eigenpair_index]))
    else
        [Float64(row[eigenpair_index]) for row in raw_eigenvectors]
    end

    length(eigenvector) == point_count || error(
        "Selected eigenvector length does not match the number of kfPoints. " *
        "Try the other --eigenvectors-as setting."
    )
    all(isfinite, eigenvector) ||
        error("Selected eigenvector contains a non-finite value")

    rows = Vector{NTuple{4, Float64}}(undef, point_count)
    for i in eachindex(points)
        x, y, z = points[i]
        rows[i] = (
            round(x; digits = 8),
            round(y; digits = 8),
            round(z; digits = 8),
            round(eigenvector[i]; digits = 12),
        )
    end
    return rows, eigenpair_index, eigenvalue
end


function html_escape(value::String)
    value = replace(value, "&" => "&amp;")
    value = replace(value, "<" => "&lt;")
    value = replace(value, ">" => "&gt;")
    return replace(value, "\"" => "&quot;")
end


function make_document(
    rows::Vector{NTuple{4, Float64}},
    eigenpair_index::Int,
    eigenvalue::Float64,
    rank::Int,
    marker_size::Float64,
    title::String,
)
    max_abs = maximum(abs(row[4]) for row in rows)
    max_abs == 0 && (max_abs = 1.0)
    compact_rows = JSON3.write(rows)

    template = """<!doctype html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>__TITLE__</title>
  <script src="__PLOTLY_CDN__"></script>
  <style>
    :root { color-scheme: light dark; }
    * { box-sizing: border-box; }
    body {
      margin: 0;
      padding: 20px;
      background: light-dark(#ffffff, #111827);
      color: light-dark(#172033, #e5e7eb);
      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
    }
    header { display: flex; align-items: baseline; justify-content: space-between; gap: 12px; flex-wrap: wrap; margin-bottom: 8px; }
    h1 { margin: 0; font-size: 1.25rem; font-weight: 600; }
    .meta { color: light-dark(#5b6472, #aeb7c5); font-variant-numeric: tabular-nums; }
    #plot { width: 100%; height: min(82vh, 850px); min-height: 520px; }
  </style>
</head>
<body>
  <header>
    <h1>__TITLE__</h1>
    <div class="meta">Eigenpair __EIGENPAIR__ · λ = __EIGENVALUE__ · closeness rank __RANK__ · __POINT_COUNT__ points</div>
  </header>
  <div id="plot" role="img" aria-label="Interactive three-dimensional k-space scatter plot colored by eigenvector amplitude."></div>

  <script>
    const rows = __ROWS__;
    const maxAbs = __MAX_ABS__;
    const markerSize = __MARKER_SIZE__;
    const plot = document.getElementById("plot");
    const darkQuery = window.matchMedia("(prefers-color-scheme: dark)");

    function draw() {
      const dark = darkQuery.matches;
      const foreground = dark ? "#e5e7eb" : "#172033";
      const muted = dark ? "#aeb7c5" : "#5b6472";
      const border = dark ? "#4b5563" : "#c9ced6";
      const background = dark ? "#111827" : "#ffffff";
      const colorscale = [[0, "#2166ac"], [0.5, background], [1, "#b2182b"]];
      const camera = plot.layout?.scene?.camera || {eye: {x: 1.45, y: 1.45, z: 1.15}};
      const axis = name => ({
        title: {text: name, font: {color: foreground, size: 16}},
        color: foreground,
        tickfont: {color: muted, size: 14},
        gridcolor: border,
        zerolinecolor: border,
        showbackground: false,
        showspikes: false
      });

      Plotly.react(plot, [{
        type: "scatter3d",
        mode: "markers",
        x: rows.map(row => row[0]),
        y: rows.map(row => row[1]),
        z: rows.map(row => row[2]),
        customdata: rows.map((row, index) => [index + 1, row[3]]),
        marker: {
          size: markerSize,
          opacity: 0.9,
          color: rows.map(row => row[3]),
          cmin: -maxAbs,
          cmax: maxAbs,
          cmid: 0,
          colorscale,
          colorbar: {
            title: {text: "Eigenvector<br>amplitude", side: "right"},
            tickformat: ".3f",
            thickness: 14,
            len: 0.72,
            outlinecolor: border,
            outlinewidth: 1
          }
        },
        hovertemplate:
          "<b>Point %{customdata[0]}</b><br>" +
          "kₓ = %{x:.5f}<br>" +
          "kᵧ = %{y:.5f}<br>" +
          "k_z = %{z:.5f}<br>" +
          "amplitude = %{customdata[1]:.6f}<extra></extra>"
      }], {
        autosize: true,
        margin: {l: 0, r: 18, t: 6, b: 0},
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: {color: foreground, size: 12},
        hoverlabel: {bgcolor: background, bordercolor: border, font: {color: foreground}},
        scene: {
          xaxis: axis("kₓ"),
          yaxis: axis("kᵧ"),
          zaxis: axis("k_z"),
          aspectmode: "data",
          camera
        },
        showlegend: false,
        uirevision: "gap-eigenvector"
      }, {
        responsive: true,
        displaylogo: false,
        scrollZoom: true,
        toImageButtonOptions: {
                  format: "png",
                  filename: "gap_eigenvector",
                  width: 1200,
                  height: 900,
                  scale: 3
                },
        modeBarButtonsToRemove: ["sendDataToCloud", "lasso2d", "select2d"]
      });
    }

    draw();
    darkQuery.addEventListener("change", draw);
  </script>
</body>
</html>
"""

    replacements = [
        "__TITLE__" => html_escape(title),
        "__PLOTLY_CDN__" => PLOTLY_CDN,
        "__EIGENPAIR__" => string(eigenpair_index),
        "__EIGENVALUE__" => @sprintf("%.10g", eigenvalue),
        "__RANK__" => string(rank),
        "__POINT_COUNT__" => string(length(rows)),
        "__ROWS__" => compact_rows,
        "__MAX_ABS__" => @sprintf("%.16g", max_abs),
        "__MARKER_SIZE__" => @sprintf("%.10g", marker_size),
    ]
    for (key, value) in replacements
        template = replace(template, key => value)
    end
    return template
end


function default_output_path(input_path::String, rank::Int)
    stem, _ = splitext(input_path)
    return "$(stem)_eigenvector_rank$(rank)_3d.html"
end


function open_in_browser(path::String)
    if Sys.isapple()
        run(`open $path`)
    elseif Sys.iswindows()
        run(`cmd /c start "" $path`)
    else
        run(`xdg-open $path`)
    end
end


function main()
    options = parse_args(ARGS)

    input_path = abspath(expanduser(options.input))
    isfile(input_path) || error("Input file not found: $input_path")
    output_path = options.output === nothing ?
        default_output_path(input_path, options.rank) :
        abspath(expanduser(options.output))

    rows, eigenpair_index, eigenvalue = extract_plot_data(
        input_path,
        options.rank,
        options.orientation,
    )
    document = make_document(
        rows,
        eigenpair_index,
        eigenvalue,
        options.rank,
        options.marker_size,
        options.title,
    )

    mkpath(dirname(output_path))
    write(output_path, document)

    println("Selected eigenpair: $eigenpair_index")
    println("Eigenvalue: ", @sprintf("%.13g", eigenvalue))
    println("Distance from 1: ", @sprintf("%.13g", abs(eigenvalue - 1.0)))
    println("Wrote: $output_path")

    options.open_result && open_in_browser(output_path)
end


try
    main()
catch exception
    println(stderr, "Could not create plot: ", sprint(showerror, exception))
    exit(1)
end
