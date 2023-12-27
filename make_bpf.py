#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

__author__ = "Juan Manuel Prieto"

"""
    Make BPF files
"""


def process_input_file(input_file):
    with open(input_file, 'r') as f:
        return f.readlines()

def main():
    parser = argparse.ArgumentParser(description='Process interaction sites.')
    parser.add_argument('input_file', help='Input PDB file')
    parser.add_argument('-p', '--positions', required=True, help='Positions to extract, e.g., 4,5,7')

    args = parser.parse_args()

    try:
        interest = [int(pos) for pos in args.positions.split(',')]
    except ValueError:
        print("Invalid input for positions. Please enter comma-separated integers.")
        exit(1)

    v_set = -2.00
    r_value = 1.20

    output_file = 'ligand_bias.bpf'

    with open(output_file, 'w') as out:
        out.write('x    y    z     Vset    r     type\n')

        for line in process_input_file(args.input_file):
            if int(line[23:27]) in interest:
                out.write('{}     {}     {}     {}\n'.format(line[30:55].strip(), v_set, r_value, line[17:20].lower()))

if __name__ == '__main__':
    main()
