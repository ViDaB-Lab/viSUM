#!/usr/bin/env python3

"""Run the Deep6 models with deterministic model-file selection."""

import argparse
import os
import sys
import warnings
from pathlib import Path


MODEL_LENGTHS = ("250", "500", "1000", "1500")
SCORE_HEADER = ("name", "length", "duplo", "euk", "mono", "pro", "ribo", "vari")


def parse_args():
    parser = argparse.ArgumentParser(
        description="Run Deep6 while selecting exactly one model per length class."
    )
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--minimum-length", required=True, type=int)
    parser.add_argument("--deep6-installation", required=True, type=Path)
    parser.add_argument("--models", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    return parser.parse_args()


def select_model(model_directory, model_length):
    matches = sorted(model_directory.glob("model_{}_*.h5".format(model_length)))
    if len(matches) != 1:
        raise ValueError(
            "Expected exactly one model_{}_*.h5 file in '{}'; found {}.".format(
                model_length, model_directory, len(matches)
            )
        )
    if matches[0].stat().st_size == 0:
        raise ValueError("Deep6 model is empty: '{}'".format(matches[0]))
    return matches[0]


def predict_batch(deep6_functions, names, forward_sequences, reverse_sequences):
    if not names:
        return
    records = zip(names, forward_sequences, reverse_sequences)
    list(map(deep6_functions.predx, records))


def main():
    args = parse_args()

    if args.minimum_length < 250:
        raise SystemExit("ERROR: Deep6 requires --minimum-length of at least 250.")

    master_directory = args.deep6_installation / "Master"
    functions_file = master_directory / "deep6_functions.py"
    if not functions_file.is_file():
        raise SystemExit(
            "ERROR: Deep6 functions file is missing: '{}'".format(functions_file)
        )

    sys.path.insert(0, str(master_directory))

    try:
        from Bio import SeqIO
        from tensorflow.keras.models import load_model
        import deep6_functions as deep6_functions
    except ImportError as error:
        raise SystemExit("ERROR: Deep6 dependency import failed: {}".format(error))

    warnings.filterwarnings("ignore", "Error in loading the saved optimizer")

    model_dictionary = {}
    for model_length in MODEL_LENGTHS:
        model_file = select_model(args.models, model_length)
        print("Loading Deep6 {}-nt model: {}".format(model_length, model_file.name))
        model_dictionary[model_length] = load_model(str(model_file), compile=False)

    deep6_functions.mdict = model_dictionary
    deep6_functions.outfile = str(args.output)

    with args.output.open("w", encoding="utf-8", newline="\n") as output_handle:
        output_handle.write("\t".join(SCORE_HEADER) + "\n")

    forward_sequences = []
    reverse_sequences = []
    sequence_names = []

    for record in SeqIO.parse(str(args.input), "fasta"):
        sequence = str(record.seq)
        if len(sequence) < args.minimum_length:
            continue

        forward_sequences.append(deep6_functions.encodex(sequence))
        reverse_sequences.append(
            deep6_functions.encodex(str(record.seq.reverse_complement()))
        )
        sequence_names.append(str(record.id))

        if len(sequence_names) == 100:
            predict_batch(
                deep6_functions,
                sequence_names,
                forward_sequences,
                reverse_sequences,
            )
            forward_sequences = []
            reverse_sequences = []
            sequence_names = []

    predict_batch(
        deep6_functions,
        sequence_names,
        forward_sequences,
        reverse_sequences,
    )

    print("Deep6 predictions completed: {}".format(args.output))


if __name__ == "__main__":
    # Prediction is CPU-only in viSUM; set before TensorFlow is imported.
    os.environ.setdefault("CUDA_VISIBLE_DEVICES", "-1")
    main()
