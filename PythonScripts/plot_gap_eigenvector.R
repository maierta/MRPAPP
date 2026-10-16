#!/usr/bin/env Rscript

# Create an interactive 3D kfPoint/eigenvector scatter plot from a gap JSON file.
#
# The Gap_test.jsn format stores eigenvectors as rows. This script defaults to
# the eigenvector whose eigenvalue is the second closest to 1. Use --rank 1 for
# the closest eigenvalue instead.

plotly_cdn <- "https://cdn.jsdelivr.net/npm/plotly.js-dist-min@2.35.2/plotly.min.js"

usage <- function(status = 0L) {
  message(paste(
    "Usage:",
    "  Rscript plot_gap_eigenvector.R INPUT [options]",
    "",
    "Options:",
    "  -o, --output PATH              Output HTML path",
    "  --rank N                       Rank by |eigenvalue - 1| (default: 2)",
    "  --eigenvectors-as rows|columns Matrix orientation (default: rows)",
    "  --marker-size NUMBER            Marker size (default: 4)",
    "  --title TEXT                    Plot title",
    "  --open                          Open the result in a browser",
    "  -h, --help                      Show this help",
    sep = "\n"
  ))
  quit(status = status)
}

parse_args <- function(args) {
  options <- list(
    input = NULL,
    output = NULL,
    rank = 2L,
    orientation = "rows",
    marker_size = 4,
    title = "Gap eigenvector on k-space points",
    open = FALSE
  )

  i <- 1L
  while (i <= length(args)) {
    arg <- args[[i]]
    if (arg %in% c("-h", "--help")) {
      usage(0L)
    } else if (arg %in% c("-o", "--output")) {
      i <- i + 1L
      if (i > length(args)) stop(arg, " requires a value", call. = FALSE)
      options$output <- args[[i]]
    } else if (arg == "--rank") {
      i <- i + 1L
      if (i > length(args)) stop("--rank requires a value", call. = FALSE)
      options$rank <- suppressWarnings(as.integer(args[[i]]))
    } else if (arg == "--eigenvectors-as") {
      i <- i + 1L
      if (i > length(args)) stop("--eigenvectors-as requires a value", call. = FALSE)
      options$orientation <- args[[i]]
    } else if (arg == "--marker-size") {
      i <- i + 1L
      if (i > length(args)) stop("--marker-size requires a value", call. = FALSE)
      options$marker_size <- suppressWarnings(as.numeric(args[[i]]))
    } else if (arg == "--title") {
      i <- i + 1L
      if (i > length(args)) stop("--title requires a value", call. = FALSE)
      options$title <- args[[i]]
    } else if (arg == "--open") {
      options$open <- TRUE
    } else if (startsWith(arg, "-")) {
      stop("Unknown option: ", arg, call. = FALSE)
    } else if (is.null(options$input)) {
      options$input <- arg
    } else {
      stop("Unexpected positional argument: ", arg, call. = FALSE)
    }
    i <- i + 1L
  }

  if (is.null(options$input)) usage(1L)
  if (is.na(options$rank) || options$rank < 1L) {
    stop("--rank must be a positive integer", call. = FALSE)
  }
  if (!(options$orientation %in% c("rows", "columns"))) {
    stop("--eigenvectors-as must be rows or columns", call. = FALSE)
  }
  if (!is.finite(options$marker_size) || options$marker_size <= 0) {
    stop("--marker-size must be positive", call. = FALSE)
  }
  options
}

html_escape <- function(value) {
  value <- gsub("&", "&amp;", value, fixed = TRUE)
  value <- gsub("<", "&lt;", value, fixed = TRUE)
  value <- gsub(">", "&gt;", value, fixed = TRUE)
  value <- gsub('"', "&quot;", value, fixed = TRUE)
  value
}

replace_literal <- function(text, pattern, replacement) {
  parts <- strsplit(text, pattern, fixed = TRUE)[[1L]]
  paste(parts, collapse = replacement)
}

default_output_path <- function(input_path, rank) {
  directory <- dirname(input_path)
  filename <- basename(input_path)
  stem <- sub("\\.[^.]*$", "", filename)
  file.path(directory, paste0(stem, "_eigenvector_rank", rank, "_3d.html"))
}

make_document <- function(rows, eigenpair_index, eigenvalue, rank, marker_size, title) {
  amplitudes <- rows[, 4L]
  max_abs <- max(abs(amplitudes))
  if (!is.finite(max_abs) || max_abs == 0) max_abs <- 1

  compact_rows <- as.character(jsonlite::toJSON(
    unname(rows),
    matrix = "rows",
    digits = 12,
    pretty = FALSE,
    auto_unbox = TRUE
  ))

  template <- paste0(
    '<!doctype html>\n',
    '<html lang="en">\n',
    '<head>\n',
    '  <meta charset="utf-8">\n',
    '  <meta name="viewport" content="width=device-width, initial-scale=1">\n',
    '  <title>__TITLE__</title>\n',
    '  <script src="__PLOTLY_CDN__"></script>\n',
    '  <style>\n',
    '    :root { color-scheme: light dark; }\n',
    '    * { box-sizing: border-box; }\n',
    '    body {\n',
    '      margin: 0; padding: 20px;\n',
    '      background: light-dark(#ffffff, #111827);\n',
    '      color: light-dark(#172033, #e5e7eb);\n',
    '      font-family: ui-sans-serif, system-ui, -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;\n',
    '    }\n',
    '    header { display: flex; align-items: baseline; justify-content: space-between; gap: 12px; flex-wrap: wrap; margin-bottom: 8px; }\n',
    '    h1 { margin: 0; font-size: 1.25rem; font-weight: 600; }\n',
    '    .meta { color: light-dark(#5b6472, #aeb7c5); font-variant-numeric: tabular-nums; }\n',
    '    #plot { width: 100%; height: min(82vh, 850px); min-height: 520px; }\n',
    '  </style>\n',
    '</head>\n',
    '<body>\n',
    '  <header>\n',
    '    <h1>__TITLE__</h1>\n',
    '    <div class="meta">Eigenpair __EIGENPAIR__ · λ = __EIGENVALUE__ · closeness rank __RANK__ · __POINT_COUNT__ points</div>\n',
    '  </header>\n',
    '  <div id="plot" role="img" aria-label="Interactive three-dimensional k-space scatter plot colored by eigenvector amplitude."></div>\n',
    '  <script>\n',
    '    const rows = __ROWS__;\n',
    '    const maxAbs = __MAX_ABS__;\n',
    '    const markerSize = __MARKER_SIZE__;\n',
    '    const plot = document.getElementById("plot");\n',
    '    const darkQuery = window.matchMedia("(prefers-color-scheme: dark)");\n',
    '    function draw() {\n',
    '      const dark = darkQuery.matches;\n',
    '      const foreground = dark ? "#e5e7eb" : "#172033";\n',
    '      const muted = dark ? "#aeb7c5" : "#5b6472";\n',
    '      const border = dark ? "#4b5563" : "#c9ced6";\n',
    '      const background = dark ? "#111827" : "#ffffff";\n',
    '      const colorscale = [[0, "#2166ac"], [0.5, background], [1, "#b2182b"]];\n',
    '      const camera = plot.layout?.scene?.camera || {eye: {x: 1.45, y: 1.45, z: 1.15}};\n',
    '      const axis = name => ({\n',
    '        title: {text: name, font: {color: foreground, size: 12}},\n',
    '        color: foreground, tickfont: {color: muted, size: 11},\n',
    '        gridcolor: border, zerolinecolor: border, showbackground: false, showspikes: false\n',
    '      });\n',
    '      Plotly.react(plot, [{\n',
    '        type: "scatter3d", mode: "markers",\n',
    '        x: rows.map(row => row[0]), y: rows.map(row => row[1]), z: rows.map(row => row[2]),\n',
    '        customdata: rows.map((row, index) => [index + 1, row[3]]),\n',
    '        marker: {\n',
    '          size: markerSize, opacity: 0.9, color: rows.map(row => row[3]),\n',
    '          cmin: -maxAbs, cmax: maxAbs, cmid: 0, colorscale,\n',
    '          colorbar: {\n',
    '            title: {text: "Eigenvector<br>amplitude", side: "right"},\n',
    '            tickformat: ".3f", thickness: 14, len: 0.72, outlinecolor: border, outlinewidth: 1\n',
    '          }\n',
    '        },\n',
    '        hovertemplate: "<b>Point %{customdata[0]}</b><br>" +\n',
    '          "kₓ = %{x:.5f}<br>kᵧ = %{y:.5f}<br>k_z = %{z:.5f}<br>" +\n',
    '          "amplitude = %{customdata[1]:.6f}<extra></extra>"\n',
    '      }], {\n',
    '        autosize: true, margin: {l: 0, r: 18, t: 6, b: 0},\n',
    '        paper_bgcolor: "rgba(0,0,0,0)", plot_bgcolor: "rgba(0,0,0,0)",\n',
    '        font: {color: foreground, size: 12},\n',
    '        hoverlabel: {bgcolor: background, bordercolor: border, font: {color: foreground}},\n',
    '        scene: {xaxis: axis("kₓ"), yaxis: axis("kᵧ"), zaxis: axis("k_z"), aspectmode: "data", camera},\n',
    '        showlegend: false, uirevision: "gap-eigenvector"\n',
    '      }, {\n',
    '        responsive: true, displaylogo: false, scrollZoom: true,\n',
    '        modeBarButtonsToRemove: ["sendDataToCloud", "lasso2d", "select2d"]\n',
    '      });\n',
    '    }\n',
    '    draw();\n',
    '    darkQuery.addEventListener("change", draw);\n',
    '  </script>\n',
    '</body>\n',
    '</html>\n'
  )

  replacements <- list(
    "__TITLE__" = html_escape(title),
    "__PLOTLY_CDN__" = plotly_cdn,
    "__EIGENPAIR__" = as.character(eigenpair_index),
    "__EIGENVALUE__" = format(eigenvalue, digits = 10L, scientific = FALSE, trim = TRUE),
    "__RANK__" = as.character(rank),
    "__POINT_COUNT__" = format(nrow(rows), big.mark = ",", scientific = FALSE),
    "__ROWS__" = compact_rows,
    "__MAX_ABS__" = format(max_abs, digits = 16L, scientific = TRUE, trim = TRUE),
    "__MARKER_SIZE__" = format(marker_size, digits = 10L, scientific = FALSE, trim = TRUE)
  )
  for (key in names(replacements)) {
    template <- replace_literal(template, key, replacements[[key]])
  }
  template
}

main <- function() {
  options <- parse_args(commandArgs(trailingOnly = TRUE))

  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop(
      "The jsonlite package is required. Install it with install.packages('jsonlite').",
      call. = FALSE
    )
  }

  input_path <- normalizePath(options$input, mustWork = TRUE)
  output_path <- if (is.null(options$output)) {
    default_output_path(input_path, options$rank)
  } else {
    path.expand(options$output)
  }

  data <- jsonlite::fromJSON(input_path, simplifyVector = TRUE)
  required <- c("kfPoints", "Eigenvalues", "Eigenvectors")
  missing <- required[!required %in% names(data)]
  if (length(missing)) {
    stop("Input is missing required field(s): ", paste(missing, collapse = ", "), call. = FALSE)
  }

  points <- as.matrix(data$kfPoints)
  eigenvalues <- as.numeric(data$Eigenvalues)
  eigenvectors <- as.matrix(data$Eigenvectors)

  if (ncol(points) != 3L) stop("kfPoints must have exactly three columns", call. = FALSE)
  if (any(!is.finite(points))) stop("kfPoints contains a non-finite value", call. = FALSE)
  if (options$rank > length(eigenvalues)) {
    stop("--rank cannot exceed the number of eigenvalues", call. = FALSE)
  }

  ordering <- order(abs(eigenvalues - 1), seq_along(eigenvalues))
  eigenpair_index <- ordering[[options$rank]]
  eigenvalue <- eigenvalues[[eigenpair_index]]
  eigenvector <- if (options$orientation == "rows") {
    as.numeric(eigenvectors[eigenpair_index, ])
  } else {
    as.numeric(eigenvectors[, eigenpair_index])
  }

  if (length(eigenvector) != nrow(points)) {
    stop(
      "Selected eigenvector length does not match the number of kfPoints. ",
      "Try the other --eigenvectors-as setting.",
      call. = FALSE
    )
  }
  if (any(!is.finite(eigenvector))) {
    stop("Selected eigenvector contains a non-finite value", call. = FALSE)
  }

  rows <- cbind(points, eigenvector)
  document <- make_document(
    rows,
    eigenpair_index,
    eigenvalue,
    options$rank,
    options$marker_size,
    options$title
  )

  output_dir <- dirname(output_path)
  if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
  writeLines(document, output_path, useBytes = TRUE)
  output_path <- normalizePath(output_path, mustWork = TRUE)

  message("Selected eigenpair: ", eigenpair_index)
  message("Eigenvalue: ", format(eigenvalue, digits = 13L))
  message("Distance from 1: ", format(abs(eigenvalue - 1), digits = 13L))
  message("Wrote: ", output_path)

  if (options$open) browseURL(paste0("file://", URLencode(output_path)))
}

tryCatch(
  main(),
  error = function(error) {
    message("Could not create plot: ", conditionMessage(error))
    quit(status = 1L)
  }
)
