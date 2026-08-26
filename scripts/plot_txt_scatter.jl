#!/usr/bin/env julia

const _RUNNING_AS_SCRIPT = abspath(PROGRAM_FILE) == @__FILE__

if _RUNNING_AS_SCRIPT
    import Pkg
    Pkg.activate(joinpath(@__DIR__, ".."); io = devnull)
end

using Cairo
using Fontconfig

const SCALE = 2.5
const WIDTH = 1600
const HEIGHT = 1600
const LEFT = 88.0 * SCALE
const RIGHT = 20.0 * SCALE
const TOP = 24.0 * SCALE
const BOTTOM = 84.0 * SCALE
const GRID_LINE_WIDTH = 1.0 * SCALE
const AXIS_LINE_WIDTH = 2.0 * SCALE
const POINT_RADIUS = 3.2 * SCALE
const TICK_LABEL_SIZE = 32.0 * SCALE
const X_TICK_LABEL_OFFSET = 24.0 * SCALE
const Y_TICK_LABEL_OFFSET = 12.0 * SCALE
const POINT_RE = r"^\(?\s*([+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?)\s*,\s*([+-]?(?:(?:\d+(?:\.\d*)?)|(?:\.\d+))(?:[eE][+-]?\d+)?)\s*\)?$"
const DENSITY_X_LIMITS = (0.0, 0.2)
const DOMAIN_SIZE_X_LIMITS = (0.0, 60.0)
const Y_LIMITS = (-0.05, 1.05)
const DENSITY_X_TICKS = [(0.0, "0"), (0.05, "0.05"), (0.1, "0.1"), (0.15, "0.15"), (0.2, "0.2")]
const DOMAIN_SIZE_X_TICKS = [(0.0, "0"), (15.0, "15"), (30.0, "30"), (45.0, "45"), (60.0, "60")]
const Y_TICKS = [(0.0, "0"), (0.25, "0.25"), (0.5, "0.5"), (0.75, "0.75"), (1.0, "1.0")]

function usage()
    return """
    Usage:
      julia --project=. scripts/plot_txt_scatter.jl input.txt [output.png]

    Reads one point per line as either "(x, y)" or "x, y".
    Empty lines and lines starting with # are ignored.
    """
end

function read_points(path::AbstractString)
    points = Tuple{Float64, Float64}[]
    for (line_number, line) in enumerate(eachline(path))
        stripped = strip(line)
        (isempty(stripped) || startswith(stripped, "#")) && continue

        matched = match(POINT_RE, stripped)
        matched === nothing && error("Invalid point at $path:$line_number: $line")
        x = parse(Float64, matched.captures[1])
        y = parse(Float64, matched.captures[2])
        isfinite(x) && isfinite(y) || error("Non-finite point at $path:$line_number: $line")
        push!(points, (x, y))
    end

    isempty(points) && error("No points found in $path")
    return points
end

function draw_text(
        ctx::CairoContext,
        x::Float64,
        y::Float64,
        text::AbstractString;
        size::Float64 = 16.0,
        halign::Symbol = :left,
        valign::Symbol = :center,
        rotation::Float64 = 0.0,
    )
    save(ctx)
    translate(ctx, x, y)
    rotation == 0.0 || rotate(ctx, rotation)
    select_font_face(ctx, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(ctx, size)
    ext = text_extents(ctx, text)
    x_offset = halign === :center ? -(ext[3] / 2 + ext[1]) :
        halign === :right ? -(ext[3] + ext[1]) : -ext[1]
    y_offset = valign === :center ? -(ext[4] / 2 + ext[2]) :
        valign === :top ? -ext[2] : -(ext[4] + ext[2])
    move_to(ctx, x_offset, y_offset)
    show_text(ctx, text)
    restore(ctx)
    return nothing
end

function x_axis(path::AbstractString)
    name = splitext(basename(path))[1]
    if startswith(name, "avg_domain_size_vs_dom_red")
        return DOMAIN_SIZE_X_LIMITS, DOMAIN_SIZE_X_TICKS
    end
    return DENSITY_X_LIMITS, DENSITY_X_TICKS
end

function draw_scatter(points, input_path::AbstractString, output_path::AbstractString)
    x_limits, x_ticks = x_axis(input_path)
    xmin, xmax = x_limits
    ymin, ymax = Y_LIMITS

    plot_width = WIDTH - LEFT - RIGHT
    plot_height = HEIGHT - TOP - BOTTOM
    plot_bottom = HEIGHT - BOTTOM
    plot_right = WIDTH - RIGHT
    scale_x(x) = LEFT + (x - xmin) / (xmax - xmin) * plot_width
    scale_y(y) = plot_bottom - (y - ymin) / (ymax - ymin) * plot_height

    surface = CairoRGBSurface(WIDTH, HEIGHT)
    ctx = CairoContext(surface)

    set_source_rgb(ctx, 1.0, 1.0, 1.0)
    paint(ctx)

    set_line_width(ctx, GRID_LINE_WIDTH)
    for (index, (tick, label)) in enumerate(x_ticks)
        x = scale_x(tick)
        set_source_rgb(ctx, 0.9, 0.9, 0.9)
        move_to(ctx, x, TOP)
        line_to(ctx, x, plot_bottom)
        stroke(ctx)
        set_source_rgb(ctx, 0.2, 0.2, 0.2)
        halign = index == 1 ? :left : index == length(x_ticks) ? :right : :center
        draw_text(ctx, x, plot_bottom + X_TICK_LABEL_OFFSET, label; size = TICK_LABEL_SIZE, halign = halign, valign = :top)
    end
    for (tick, label) in Y_TICKS
        y = scale_y(tick)
        set_source_rgb(ctx, 0.9, 0.9, 0.9)
        move_to(ctx, LEFT, y)
        line_to(ctx, plot_right, y)
        stroke(ctx)
        set_source_rgb(ctx, 0.2, 0.2, 0.2)
        draw_text(ctx, LEFT - Y_TICK_LABEL_OFFSET, y, label; size = TICK_LABEL_SIZE, halign = :right)
    end

    set_source_rgb(ctx, 0.0, 0.0, 0.0)
    set_line_width(ctx, AXIS_LINE_WIDTH)
    move_to(ctx, LEFT, TOP)
    line_to(ctx, LEFT, plot_bottom)
    line_to(ctx, plot_right, plot_bottom)
    stroke(ctx)

    set_source_rgb(ctx, 0.0, 0.0, 0.0)
    for (x, y) in points
        arc(ctx, scale_x(x), scale_y(y), POINT_RADIUS, 0.0, 2 * pi)
        fill(ctx)
    end

    mkpath(dirname(output_path))
    write_to_png(surface, output_path)
    return output_path
end

function default_output_path(input_path::AbstractString)
    stem, _ = splitext(input_path)
    return stem * ".png"
end

function main(args::Vector{String} = copy(ARGS))
    if isempty(args) || args[1] in ("-h", "--help")
        println(usage())
        return nothing
    end
    1 <= length(args) <= 2 || error("Expected input.txt and optional output.png")

    input_path = abspath(args[1])
    output_path = length(args) == 2 ? abspath(args[2]) : default_output_path(input_path)
    points = read_points(input_path)
    draw_scatter(points, input_path, output_path)
    println("saved plot: $output_path")
    return output_path
end

if _RUNNING_AS_SCRIPT
    main(copy(ARGS))
end
