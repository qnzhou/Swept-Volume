import sweep3d
import lagrange
import argparse
import pathlib

def parse_args():
    parser = argparse.ArgumentParser(description="Compute generalized sweep from configuration files.")
    parser.add_argument("function_file", type=str, help="Path to the space-time function file.")
    parser.add_argument("config_file", type=str, help="Path to the configuration file.")
    parser.add_argument("-o", "--output-dir", type=str, default="output", help="Directory to save output files.")
    return parser.parse_args()

def main():
    args = parse_args()
    r = sweep3d.generalized_sweep_from_config(args.function_file, args.config_file)

    result_dir = pathlib.Path(args.output_dir)
    envelope_file = result_dir / "envelope.msh"
    arrangement_file = result_dir / "arrangement.msh"
    sweep_surface_file = result_dir / "sweep_surface.msh"

    lagrange.io.save_mesh(envelope_file, r.envelope)
    lagrange.io.save_mesh(arrangement_file, r.arrangement)
    lagrange.io.save_mesh(sweep_surface_file, r.sweep_surface)

if __name__ == "__main__":
    main()
