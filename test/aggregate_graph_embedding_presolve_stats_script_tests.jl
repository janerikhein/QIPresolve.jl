using Test
import CSV

include(joinpath(@__DIR__, "..", "scripts", "aggregate_graph_embedding_presolve_stats.jl"))
const AggregatePresolveStatsScript = Main.AggregateGraphEmbeddingPresolveStats

const EXPECTED_AGGREGATE_HEADER = [
    "type",
    "n",
    "R",
    "num_anchors",
    "alpha",
    "edge_density",
    "infeas_strategy",
    "infeas_base",
    "box_scale",
    "num_instances",
    "log_domain_sum_orig",
    "log_domain_sum_ps",
    "relative_log_domain_reduction_percent",
    "num_quadratic_constraints_orig",
    "num_quadratic_constraints_ps",
    "num_linear_constraints_orig",
    "num_linear_constraints_ps",
    "parity_presolve_time",
    "residue_presolve_time",
    "num_parity_substitutions",
    "num_pattern_substitutions",
    "num_variables_fixed",
    "num_constraints_tightened",
    "num_constraints_eq_tightened",
    "avg_relative_bound_tightening_percent",
    "infeasibility_detected_count",
    "infeasibility_source_elimination_count",
    "infeasibility_source_propagation_count",
    "infeasibility_source_residue_count",
]

const AGGREGATE_INPUT_COLUMNS = [
    "instance_name";
    EXPECTED_AGGREGATE_HEADER[1:9];
    AggregatePresolveStatsScript.AVERAGE_COLUMNS;
    ["infeasibility_detected", "infeasibility_source"];
]

function aggregate_input_row(
        instance_name::AbstractString;
        type_label::AbstractString = "2_connected",
        n::Int = 20,
        R::Int = 100,
        num_anchors::Int = 0,
        alpha::Float64 = 0.0,
        edge_density = 0.2,
        infeas_strategy = missing,
        infeas_base = missing,
        box_scale = missing,
        metric_offset::Float64 = 0.0,
        metric_overrides = Dict{String, Any}(),
        detected::Bool = false,
        source::AbstractString = "none",
    )
    values = Dict{String, Any}(
        "instance_name" => instance_name,
        "type" => type_label,
        "n" => n,
        "R" => R,
        "num_anchors" => num_anchors,
        "alpha" => alpha,
        "edge_density" => edge_density,
        "infeas_strategy" => infeas_strategy,
        "infeas_base" => infeas_base,
        "box_scale" => box_scale,
        "infeasibility_detected" => detected,
        "infeasibility_source" => source,
    )
    for (metric_index, column) in enumerate(AggregatePresolveStatsScript.AVERAGE_COLUMNS)
        default_value = if column == "log_domain_sum_orig"
            metric_offset + 20.0
        elseif column == "log_domain_sum_ps"
            metric_offset + 10.0
        elseif column == "avg_relative_bound_tightening"
            0.1 + metric_offset / 100.0
        else
            metric_offset + metric_index
        end
        values[column] = get(metric_overrides, column, default_value)
    end
    names = Tuple(Symbol.(AGGREGATE_INPUT_COLUMNS))
    return NamedTuple{names}(Tuple(values[column] for column in AGGREGATE_INPUT_COLUMNS))
end

function aggregate_config(input::AbstractString, output::AbstractString)
    return AggregatePresolveStatsScript.CliConfig(
        input = String(input),
        output = String(output),
    )
end

@testset "aggregate graph embedding presolve stats script" begin
    mktempdir() do dir
        input_csv = joinpath(dir, "instances_with_presolve_stats.csv")
        output_csv = joinpath(dir, "nested", "aggregated.csv")

        CSV.write(
            input_csv,
            [
                aggregate_input_row("group_one_1"; metric_offset = 0.0),
                aggregate_input_row(
                    "group_two_1";
                    type_label = "infeasible",
                    n = 30,
                    num_anchors = 3,
                    edge_density = missing,
                    infeas_strategy = "bounding_box",
                    infeas_base = "globally_rigid",
                    box_scale = 0.756,
                    alpha = 0.126,
                    metric_offset = 10.126,
                ),
                aggregate_input_row(
                    "group_one_2";
                    metric_offset = 2.0,
                    detected = true,
                    source = "elimination",
                ),
                aggregate_input_row(
                    "group_one_3";
                    metric_offset = 4.0,
                    detected = true,
                    source = "propagation",
                ),
                aggregate_input_row(
                    "group_one_4";
                    metric_offset = 6.0,
                    detected = true,
                    source = "residue",
                ),
            ],
        )

        result_path = AggregatePresolveStatsScript.main([
            "--input=$(input_csv)",
            "--output",
            output_csv,
        ])

        @test result_path == output_csv
        @test isfile(output_csv)

        output = CSV.File(output_csv)
        @test String.(propertynames(output)) == EXPECTED_AGGREGATE_HEADER
        rows = collect(output)
        @test length(rows) == 2

        first_group = rows[1]
        @test first_group.type == "2_connected"
        @test first_group.n == 20
        @test first_group.edge_density == 0.2
        @test ismissing(first_group.infeas_strategy)
        @test ismissing(first_group.infeas_base)
        @test ismissing(first_group.box_scale)
        @test first_group.num_instances == 4
        for (metric_index, column) in enumerate(AggregatePresolveStatsScript.AVERAGE_COLUMNS)
            column == "avg_relative_bound_tightening" && continue
            expected_value = if column == "log_domain_sum_orig"
                23.0
            elseif column == "log_domain_sum_ps"
                13.0
            else
                metric_index + 3.0
            end
            @test getproperty(first_group, Symbol(column)) == expected_value
        end
        @test first_group.relative_log_domain_reduction_percent ==
            round(100.0 * (23.0 - 13.0) / 23.0; digits = 3)
        @test first_group.avg_relative_bound_tightening_percent == 13.0
        @test first_group.infeasibility_detected_count == 3
        @test first_group.infeasibility_source_elimination_count == 1
        @test first_group.infeasibility_source_propagation_count == 1
        @test first_group.infeasibility_source_residue_count == 1

        second_group = rows[2]
        @test second_group.type == "infeasible"
        @test second_group.n == 30
        @test second_group.num_anchors == 3
        @test second_group.alpha == 0.126
        @test ismissing(second_group.edge_density)
        @test second_group.infeas_strategy == "bounding_box"
        @test second_group.infeas_base == "globally_rigid"
        @test second_group.box_scale == 0.756
        @test second_group.num_instances == 1
        @test second_group.parity_presolve_time == 17.126
        @test second_group.avg_relative_bound_tightening_percent == 20.126
        @test second_group.relative_log_domain_reduction_percent == round(
            100.0 * (30.126 - 20.126) / 30.126;
            digits = 3,
        )
        @test second_group.infeasibility_detected_count == 0
        @test second_group.infeasibility_source_elimination_count == 0
        @test second_group.infeasibility_source_propagation_count == 0
        @test second_group.infeasibility_source_residue_count == 0
    end
end

@testset "aggregate graph embedding presolve stats validation" begin
    mktempdir() do dir
        output_csv = joinpath(dir, "output.csv")

        missing_column_csv = joinpath(dir, "missing_column.csv")
        write(missing_column_csv, "instance_name\nmissing_columns\n")
        @test_throws ErrorException AggregatePresolveStatsScript.aggregate_stats(
            aggregate_config(missing_column_csv, output_csv)
        )

        unknown_source_csv = joinpath(dir, "unknown_source.csv")
        CSV.write(
            unknown_source_csv,
            [aggregate_input_row(
                "unknown_source";
                detected = true,
                source = "unexpected",
            )],
        )
        @test_throws ErrorException AggregatePresolveStatsScript.aggregate_stats(
            aggregate_config(unknown_source_csv, output_csv)
        )

        inconsistent_csv = joinpath(dir, "inconsistent.csv")
        CSV.write(
            inconsistent_csv,
            [aggregate_input_row(
                "inconsistent";
                detected = false,
                source = "propagation",
            )],
        )
        @test_throws ErrorException AggregatePresolveStatsScript.aggregate_stats(
            aggregate_config(inconsistent_csv, output_csv)
        )

        missing_metric_csv = joinpath(dir, "missing_metric.csv")
        CSV.write(
            missing_metric_csv,
            [aggregate_input_row(
                "missing_metric";
                metric_overrides = Dict("log_domain_sum_orig" => missing),
            )],
        )
        @test_throws ErrorException AggregatePresolveStatsScript.aggregate_stats(
            aggregate_config(missing_metric_csv, output_csv)
        )

        nonnumeric_metric_csv = joinpath(dir, "nonnumeric_metric.csv")
        CSV.write(
            nonnumeric_metric_csv,
            [aggregate_input_row(
                "nonnumeric_metric";
                metric_overrides = Dict("log_domain_sum_orig" => "not-a-number"),
            )],
        )
        @test_throws ErrorException AggregatePresolveStatsScript.aggregate_stats(
            aggregate_config(nonnumeric_metric_csv, output_csv)
        )

        zero_original_csv = joinpath(dir, "zero_original.csv")
        CSV.write(
            zero_original_csv,
            [aggregate_input_row(
                "zero_original";
                metric_overrides = Dict("log_domain_sum_orig" => 0.0),
            )],
        )
        @test_throws ErrorException AggregatePresolveStatsScript.aggregate_stats(
            aggregate_config(zero_original_csv, output_csv)
        )
    end
end
