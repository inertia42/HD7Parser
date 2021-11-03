'''
Author: inertia
Date: 2021.10.18
Version: 0.1
Requirement:
    - toml
    - numpy
    - matplotlib

    The purpose of this script is to automatically analyse outputs of HD7 code which includes parameters, eigenvalues and mode structure.
'''

import re
import os
import copy

import numpy as np
import collections


def parameter_reader(filename):
    '''Reads and Formats the parameters from different cases.

    Args:
        filename: the name and dir of outs file.

    Returns:
        A list in which every element is a dict.
    '''

    f = open(filename, 'r')
    para_dict = {}
    para_list = []
    # save the text of previous line
    previous_line = ''

    for line in f.readlines():
        # for the lines about predefined parameters
        if "=" in line:
            # package a function for different patterns
            # but I can write a more general pattern to substitute for this function.
            def figure_reader(pattern, line):
                index_data = re.findall(pattern, line)  # find the value
                # remove those value for find index
                line = re.sub(pattern, '', line)
                # get the index name
                indexes = re.findall(r'(\S+?)\s*?[=$]', line)
                return index_data, indexes

            index_data, indexes = figure_reader(r'(0\.\d+?E[+-]\d+)', line)
            if len(index_data) != len(indexes):
                # for some special cases with a unuseful line, for example HD7-em-tok
                index_data, indexes = figure_reader(r'(-?\d\.\d+)', line)
                if len(index_data) != len(indexes):
                    raise ValueError

            for i, index in enumerate(indexes):
                para_dict[index] = float(index_data[i])
            continue
        else:
            eigen_data = re.findall(r'(-?0\.\d+?E[+-]\d+)', line)
            if len(eigen_data) == 0:
                # store the text of previous line until getting the variable name
                previous_line = line
                continue
            elif previous_line:
                indexes = re.findall(r'\S+', previous_line)
                var_name = indexes[0]
                omr_index = indexes.index("omr")
                omi_index = indexes.index("omi")
                previous_line = ''
            para_dict[var_name] = float(eigen_data[0])
            para_dict['omr'] = float(eigen_data[omr_index])
            para_dict['omi'] = float(eigen_data[omi_index])
            para_list.append(copy.deepcopy(para_dict))
    f.close()
    return para_list


def mode_reader(filename):
    '''Reads and Formats the mode structure of different parameters.

    Args:
        filename: the name and dir of the outp/outa/oute file.

    Return:

    '''
    with open(filename, 'r') as f:
        fulltext = f.read()
    cases = re.split(r'\n{2,}', fulltext)
    mode_list = []
    for case in cases:
        mode_value = re.findall(r'(-?0\.\d+?E[+-]\d+)', case)
        if ('=' not in case) & (len(mode_value) != 0):
            x = []
            real = []
            imag = []
            lines = re.split(r'\n', case)
            for line in lines:
                mode_value = re.findall(r'(-?0\.\d+?E[+-]\d+)', line)
                if len(mode_value) != 3:
                    # raise ValueError
                    continue
                x.append(float(mode_value[0]))
                real.append(float(mode_value[1]))
                imag.append(float(mode_value[2]))
                if len(x) == 0:
                    raise ValueError
            # mode = {"x": x, "real": real, "imag": imag}
            mode = [x, real, imag]
            mode_list.append(mode)
    return mode_list


def ask_path():
    '''
        This function is used to get the path of the data directory with tkinter.
    '''
    from tkinter import Tk
    from tkinter.filedialog import askdirectory

    Tk().withdraw()
    directory = askdirectory()
    return directory


if __name__ == "__main__":
    # para = parameter_reader("./outs-1.dat")
    # print(para)
    pass
    # mode = mode_reader("./outp.dat")
    # print(mode)
