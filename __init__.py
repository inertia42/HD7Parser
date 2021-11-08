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

import toml


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


def ask_path(title=''):
    '''
        This function is used to get the path of the data directory with tkinter.
    '''
    from tkinter import Tk
    from tkinter.filedialog import askdirectory

    Tk().withdraw()
    if title:
        directory = askdirectory(title=title)
    else:
        directory = askdirectory()
    return directory


def get_all_cases(path):
    '''
        This function is used to get cases of a directory or subdirectories of the directory.

        Args:
            path: the path of the directory of input data.
    '''
    cases = []
    file_name = 'outs-1.dat'
    for root, dirs, files in os.walk(path):
        if file_name in files:
            cases = parameter_reader(os.path.join(path, file_name))
            break
        for dir in dirs:
            file_path = os.path.join(path, dir, file_name)
            if os.path.exists(file_path):
                cases.append(parameter_reader(file_path))
        cases = [para for case in cases for para in case]
        break
    return cases


def datef_generator(dir: dict, anchor: dict, target, step, num, addition=''):
    '''
        This function is used to generate datef file from one data point.

        Args:
            dir: a dictionary of the input path and save path.
            anchor: the start point of the scan.
            target: the end point of the scan, it is a list or a value.
            step: the step of scan.
            num: the corresponding number of the scan variable.
            addition: the additional information of the scan.

        Raises:
            ValueError: if more than one case is met.
    '''
    if dir['ask'] == 1:
        dir['input'] = ask_path(title='Please select the input directory')
        dir['save'] = ask_path(title='Please select the save directory')
    elif dir['ask'] == 0:
        if 'input' not in dir.keys():
            raise Exception("Please input the input path.")
        if 'save' not in dir.keys():
            raise Exception("Please input the save path.")

    if isinstance(target, int) or isinstance(target, float):
        target = [target]

    anchor_name = anchor['name']
    var_name = anchor['variable']
    if isinstance(anchor['value'], int) or isinstance(anchor['value'], float):
        anchor_values = [anchor['value']]
    else:
        anchor_values = anchor['value']
    anchor_length = len(anchor_values)

    cases = get_all_cases(dir['input'])
    for anchor_value in anchor_values:
        anchor_cases = [
            case for case in cases if case[anchor_name] == anchor_value
        ]
        if len(anchor_cases) != 1:
            raise ValueError(f'The number of cases is {len(anchor_cases)} not 1')
        case = anchor_cases[0]

        for target_value in target:
            if target_value < case[var_name]:
                step_value = -step
            else:
                step_value = step
            ip1 = int((target_value - case[var_name]) / step_value + 1)
            dataf = ''' $inp
 omr={omr},omi={omi},
 num1={num1},dp1={dp1},ip1={ip1},
 p1={etai},p2={beta},p3={shat},p4={etae},p5={q},p6=0.125,
 p7={rnr},p8={tau},p9=0.00,p10={aky},p11=0.,
 p12=1.00,p13=60,p14=6.00,p15=6.0,p16={rbr},p17={epsil},p18=0.2,p19={sw},ifplt=1,
 deco=0.3
 $'''.format(num1=num, dp1=step_value, ip1=ip1, **case)
            if anchor_length == 1:
                data_dir = '{name}_{value}_to_{target}_step_{step}{addition}'.format(
                    name=var_name,
                    value=case[var_name],
                    target=target_value,
                    step=step,
                    addition=addition)
            elif anchor_length > 1:
                data_dir = '{anchor}_{anchor_value}_{name}_{value}_to_{target}_step_{step}{addition}'.format(
                    anchor=anchor_name,
                    anchor_value=anchor_value,
                    name=var_name,
                    value=case[var_name],
                    target=target_value,
                    step=step,
                    addition=addition)
            data_dir = os.path.join(dir['save'], data_dir)
            try:
                os.mkdir(data_dir)
            except FileExistsError:
                pass
            with open(os.path.join(data_dir, 'datef.dat'), 'w') as f:
                f.write(dataf)
    return


def toml_parser(toml_file):
    '''
        This function is used to parse the toml file.

        Args:
            toml_file: the path of the toml file.
    '''
    with open(toml_file, 'r') as f:
        toml_dict = toml.load(f)
    datef_generator(**toml_dict)
    return


if __name__ == "__main__":
    pass
