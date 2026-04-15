#!/usr/bin/env julia

using Pkg
Pkg.activate(joinpath(@__DIR__, ".."))

using Cairo
using Fontconfig
using Printf


const OUTPUT_PATH = joinpath(@__DIR__, "..", "results", "selected_columns_plot.png")
const X_LABEL = "Total edges"
const Y_LABEL = "Average domain size"
const X_COL_IDX = 2
const X_TICK_STEP = 50.0
const Y_TICK_STEP = 1000.0

# All series in one plot are assumed to use the same x-axis semantics at X_COL_IDX.
const SERIES = [
    (
        path = joinpath(@__DIR__, "..", "results", "random_graph_extra_edges.csv"),
        y_col_idx = 3,
        label = "2-connected before",
        color = "red",
    ),
    (
        path = joinpath(@__DIR__, "..", "results", "random_graph_extra_edges.csv"),
        y_col_idx = 4,
        label = "2-connected after",
        color = "blue",
    ),
    (
        path = joinpath(@__DIR__, "..", "results", "random_laman_extra_edges.csv"),
        y_col_idx = 3,
        label = "Laman before",
        color = "darkgreen",
    ),
    (
        path = joinpath(@__DIR__, "..", "results", "random_laman_extra_edges.csv"),
        y_col_idx = 4,
        label = "Laman after",
        color = "orange",
    ),
]

const CANVAS_WIDTH = 1200
const CANVAS_HEIGHT = 800

const LEFT_MARGIN = 0.12
const RIGHT_MARGIN = 0.06
const TOP_MARGIN = 0.12
const BOTTOM_MARGIN = 0.14

const PLOT_LEFT = LEFT_MARGIN
const PLOT_RIGHT = 1.0 - RIGHT_MARGIN
const PLOT_TOP = TOP_MARGIN
const PLOT_BOTTOM = 1.0 - BOTTOM_MARGIN
const PLOT_WIDTH = PLOT_RIGHT - PLOT_LEFT
const PLOT_HEIGHT = PLOT_BOTTOM - PLOT_TOP


function parse_csv_line(line::String)
    return strip.(split(line, ","))
end


function maybe_parse_float(token::AbstractString)
    parsed = tryparse(Float64, strip(token))
    return parsed
end


function format_tick(value::Float64)
    abs_value = abs(value)
    if abs_value ≥ 1.0e4 || (abs_value > 0 && abs_value < 1.0e-2)
        return @sprintf("%.2e", value)
    elseif abs_value ≥ 100
        return @sprintf("%.0f", value)
    elseif abs_value ≥ 10
        return @sprintf("%.1f", value)
    else
        return @sprintf("%.2f", value)
    end
end


function axis_limits(values::Vector{Float64}, step::Float64)
    isempty(values) && error("Cannot determine axis limits from empty data.")
    vmax = maximum(values)
    step > 0 || error("Axis tick step must be positive.")

    upper = max(step, ceil(vmax / step) * step)
    return 0.0, upper
end


function scale_x(x::Float64, xmin::Float64, xmax::Float64)
    return PLOT_LEFT + (x - xmin) / (xmax - xmin) * PLOT_WIDTH
end


function scale_y(y::Float64, ymin::Float64, ymax::Float64)
    return PLOT_BOTTOM - (y - ymin) / (ymax - ymin) * PLOT_HEIGHT
end


function read_series(spec)
    lines = readlines(spec.path)
    length(lines) ≥ 2 || error("File $(spec.path) must contain at least a config line and a header.")

    header = parse_csv_line(lines[2])
    x_col_name = X_COL_IDX <= length(header) ? header[X_COL_IDX] : "col_$X_COL_IDX"
    y_col_name = spec.y_col_idx <= length(header) ? header[spec.y_col_idx] : "col_$(spec.y_col_idx)"

    points = Tuple{Float64, Float64}[]
    for line in lines[3:end]
        isempty(strip(line)) && continue
        startswith(strip(line), "#") && continue

        fields = parse_csv_line(line)
        max(X_COL_IDX, spec.y_col_idx) <= length(fields) || continue

        x = maybe_parse_float(fields[X_COL_IDX])
        y = maybe_parse_float(fields[spec.y_col_idx])
        x === nothing && continue
        y === nothing && continue
        isfinite(x) || continue
        isfinite(y) || continue

        push!(points, (x, y))
    end

    isempty(points) && error("No finite points found for $(spec.label) in $(spec.path).")
    sort!(points, by = first)

    return (
        label = spec.label,
        color = spec.color,
        path = spec.path,
        x_col_name = x_col_name,
        y_col_name = y_col_name,
        points = points,
    )
end


function color_rgb(name::String)
    colors = Dict(
        "black" => (0.0, 0.0, 0.0),
        "red" => (0.85, 0.15, 0.15),
        "blue" => (0.2, 0.35, 0.85),
        "darkgreen" => (0.1, 0.5, 0.2),
        "orange" => (0.9, 0.55, 0.1),
        "purple" => (0.55, 0.3, 0.75),
        "gray" => (0.4, 0.4, 0.4),
    )
    haskey(colors, name) || error("Unsupported color name: $name")
    return colors[name]
end


to_px_x(x::Float64) = x * CANVAS_WIDTH
to_px_y(y::Float64) = y * CANVAS_HEIGHT


function set_source_rgb_name(ctx::CairoContext, color_name::String)
    r, g, b = color_rgb(color_name)
    set_source_rgb(ctx, r, g, b)
    return nothing
end


function draw_text(
        ctx::CairoContext, x::Float64, y::Float64, value::AbstractString;
        size::Float64 = 20.0, halign::Symbol = :left, valign::Symbol = :center,
        rotation::Float64 = 0.0
    )
    save(ctx)
    translate(ctx, to_px_x(x), to_px_y(y))
    rotation == 0.0 || rotate(ctx, rotation)
    select_font_face(ctx, "Sans", Cairo.FONT_SLANT_NORMAL, Cairo.FONT_WEIGHT_NORMAL)
    set_font_size(ctx, size)
    extents = text_extents(ctx, value)
    x_bearing = extents[1]
    y_bearing = extents[2]
    width = extents[3]
    height = extents[4]

    x_offset = halign === :center ? -(width / 2 + x_bearing) :
        halign === :right ? -(width + x_bearing) : -x_bearing
    y_offset = valign === :center ? -(height / 2 + y_bearing) :
        valign === :top ? -y_bearing : -(height + y_bearing)

    move_to(ctx, x_offset, y_offset)
    show_text(ctx, value)
    restore(ctx)
    return nothing
end


function draw_line(
        ctx::CairoContext, points::Vector{Tuple{Float64, Float64}}, color_name::String;
        width::Float64 = 2.0
    )
    isempty(points) && return nothing

    set_source_rgb_name(ctx, color_name)
    set_line_width(ctx, width)
    move_to(ctx, to_px_x(points[1][1]), to_px_y(points[1][2]))
    for point in points[2:end]
        line_to(ctx, to_px_x(point[1]), to_px_y(point[2]))
    end
    stroke(ctx)
    return nothing
end


function draw_plot(ctx::CairoContext, series_data)
    all_x = Float64[]
    all_y = Float64[]
    for series in series_data
        append!(all_x, first.(series.points))
        append!(all_y, last.(series.points))
    end

    xmin, xmax = axis_limits(all_x, X_TICK_STEP)
    ymin, ymax = axis_limits(all_y, Y_TICK_STEP)

    set_source_rgb(ctx, 1.0, 1.0, 1.0)
    paint(ctx)

    set_source_rgb(ctx, 0.0, 0.0, 0.0)
    set_line_width(ctx, 2.0)
    move_to(ctx, to_px_x(PLOT_LEFT), to_px_y(PLOT_TOP))
    line_to(ctx, to_px_x(PLOT_LEFT), to_px_y(PLOT_BOTTOM))
    stroke(ctx)
    move_to(ctx, to_px_x(PLOT_LEFT), to_px_y(PLOT_BOTTOM))
    line_to(ctx, to_px_x(PLOT_RIGHT), to_px_y(PLOT_BOTTOM))
    stroke(ctx)

    x_ticks = collect(xmin:X_TICK_STEP:xmax)
    for tick in x_ticks
        x_pos = scale_x(tick, xmin, xmax)
        move_to(ctx, to_px_x(x_pos), to_px_y(PLOT_BOTTOM))
        line_to(ctx, to_px_x(x_pos), to_px_y(PLOT_BOTTOM + 0.01))
        stroke(ctx)
        draw_text(ctx, x_pos, PLOT_BOTTOM + 0.03, format_tick(tick); size = 16.0, halign = :center, valign = :top)
    end

    y_ticks = collect(ymin:Y_TICK_STEP:ymax)
    for tick in y_ticks
        y_pos = scale_y(tick, ymin, ymax)
        move_to(ctx, to_px_x(PLOT_LEFT - 0.01), to_px_y(y_pos))
        line_to(ctx, to_px_x(PLOT_LEFT), to_px_y(y_pos))
        stroke(ctx)
        draw_text(ctx, PLOT_LEFT - 0.015, y_pos, format_tick(tick); size = 16.0, halign = :right, valign = :center)
    end

    for series in series_data
        scaled_points = [(scale_x(x, xmin, xmax), scale_y(y, ymin, ymax)) for (x, y) in series.points]
        length(scaled_points) == 1 && push!(scaled_points, scaled_points[1])
        draw_line(ctx, scaled_points, series.color; width = 3.0)
    end

    draw_text(ctx, (PLOT_LEFT + PLOT_RIGHT) / 2, 0.965, X_LABEL; size = 20.0, halign = :center, valign = :center)
    draw_text(ctx, 0.03, (PLOT_TOP + PLOT_BOTTOM) / 2, Y_LABEL; size = 20.0, halign = :center, valign = :center, rotation = -pi / 2)

    legend_x = PLOT_RIGHT - 0.2
    legend_y0 = PLOT_TOP + 0.03
    for (idx, series) in enumerate(series_data)
        y = legend_y0 + 0.04 * (idx - 1)
        draw_line(ctx, [(legend_x, y), (legend_x + 0.03, y)], series.color; width = 4.0)
        draw_text(ctx, legend_x + 0.04, y, series.label; size = 16.0, halign = :left, valign = :center)
    end

    draw_text(ctx, PLOT_LEFT, 0.085, "x: $(series_data[1].x_col_name) (column $X_COL_IDX)"; size = 14.0, halign = :left, valign = :center)
    return nothing
end


function main()
    series_data = [read_series(spec) for spec in SERIES]
    mkpath(dirname(OUTPUT_PATH))
    surface = CairoRGBSurface(CANVAS_WIDTH, CANVAS_HEIGHT)
    ctx = CairoContext(surface)
    draw_plot(ctx, series_data)
    write_to_png(surface, OUTPUT_PATH)
    println("saved plot: $OUTPUT_PATH")
end


main()
