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
from math import sqrt

import toml
import json
import numpy as np


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


def jl_data_reader(path, name="out.json"):
    '''
        This function is used to read the data from hd7.jl.

        Args:
            path: the path of the directory of input data.
            name: the path of the json file.
    '''
    filename = os.path.join(path, name)
    with open(filename, 'r') as f:
        data = json.load(f)
    for case in data:
        th = []
        shat = case["shat"]
        aky = case["aky"]
        for k in case["k"]:
            th.append(k/(shat*aky))
        case["phi"] = (th, case["phi_r"], case["phi_i"])
        case["apara"] = (th, case["apara_r"], case["apara_i"])
    return data


def single_parameter_extractor(cases, para_name, lower=None, upper=None):
    '''
        extract three lists from cases.
        Args:
            cases: a list of dicts.
            para_name: the name of the parameter.
        Returns:
            three lists of parameter values , omr and omi.
    '''
    para_list = []
    omr_list = []
    omi_list = []

    # check the lower and upper bound
    if lower:
        cases = [case for case in cases if case[para_name] >= lower]
    if upper:
        cases = [case for case in cases if case[para_name] <= upper]

    # sort the cases by parameter value
    cases = sorted(cases, key=lambda x: x[para_name])
    for case in cases:
        # exclude the duplicated cases
        if case[para_name] in para_list:
            continue
        para_list.append(case[para_name])
        omr_list.append(case["omr"])
        omi_list.append(case["omi"])
    return para_list, omr_list, omi_list


def two_parameter_extractor(cases, var_name='', para_name='', lower=None, upper=None):
    '''
        extract two parameter lists from cases.
        Args:
            cases: a list of dicts.
            para_name: the name of the parameter, which is equal to y.
            var_name: the name of the variable, which is equal to x.
            lower: the lower bound of the parameter.
            upper: the upper bound of the parameter.
        Returns:
            two lists of parameter values and variable parameter values.
    '''
    para_list = []
    var_list = []

    # check the lower and upper bound
    if lower:
        cases = [case for case in cases if case[var_name] >= lower]
    if upper:
        cases = [case for case in cases if case[var_name] <= upper]

    # sort the cases by parameter value
    cases = sorted(cases, key=lambda x: x[var_name])
    for case in cases:
        # exclude the duplicated cases
        if case[var_name] in para_list:
            continue
        para_list.append(case[para_name])
        var_list.append(case[var_name])
    return var_list, para_list


def cal_diffusion_coeff(cases):
    '''
        calculate the diffusion coefficient from the case.
        The method is eq.18 in "Multiple ion temperature gradient driven modes in transport barriers", doi: https://doi.org/10.1088/1741-4326/aa5d02
    '''
    for case in cases:
        shat = case["shat"]
        aky = case["aky"]
        omi = case["omi"]
        rnr = case["rnr"]
        phi = case["phi"]
        q = case["q"]
        beta = case["beta"]
        tau = case["tau"]
        etai = case["etai"]
        alpha = (q**2)*beta*(1+tau)*(1+etai)/rnr
        apara = case["apara"]
        th = np.array(phi[0])
        phi2 = np.array(phi[1])**2 + np.array(phi[2])**2
        tphi = (th**2)*phi2
        sqrt_t = np.sqrt(np.trapz(tphi, x=th)/np.trapz(phi2, x=th))
        alpha
        dcoeff = omi/(aky*rnr*(shat*sqrt_t-alpha*np.sin(sqrt_t))**2)
        case["phi_dcoeff"] = dcoeff
        apara2 = np.array(apara[1])**2 + np.array(apara[2])**2
        tapara = (th**2)*apara2
        sqrt_t = np.sqrt(np.trapz(tapara, x=th)/np.trapz(apara2, x=th))
        dcoeff = omi/(aky*rnr*(shat*sqrt_t-alpha*np.sin(sqrt_t))**2)
        case["apara_dcoeff"] = dcoeff
    return cases


def cal_psi_with_alpha(case):
    '''
    calculate the mode distribution of psi with alpha.
    '''
    q = case['q']
    tau = case['tau']
    aky = case['aky']
    dth = case['dth']*0.447/aky
    omr = case['omr']
    omi = case['omi']
    rnr = case['rnr']
    epsil = case['epsil']
    alpha = sqrt(1+(epsil/q)**2)
    # mode_data = case['apara']
    th, a_real, a_imag = case['apara']
    # th = mode_data['x']
    # a_real = mode_data['real']
    # a_imag = mode_data['imag']

    psi_real = []
    psi_imag = []
    for i, _ in enumerate(th):
        real_integral = sum(a_real[i:])*2-a_real[i]-a_real[-1]
        imag_integral = sum(a_imag[i:])*2-a_imag[i]-a_imag[-1]
        # psi_real.append(q*omr*sqrt(tau)*dk*imag_integral/(2*rnr*shat)+q*omi*sqrt(tau)*dk*real_integral/(2*rnr*shat))
        psi_real.append(q*sqrt(tau)*dth*alpha*aky*(imag_integral*omr+real_integral*omi)/(4*rnr))
        psi_imag.append(q*sqrt(tau)*dth*alpha*aky*(-real_integral*omr+imag_integral*omi)/(4*rnr))
    # case['psi'] = {'x': th, 'real': psi_real, 'imag': psi_imag}
    case['psi'] = [th, psi_real, psi_imag]
    return case


def merge_data(dir):
    file_para = os.path.join(dir, "outs-1.dat")
    file_phi = os.path.join(dir, "outp-1.dat")
    file_apara = os.path.join(dir, "outa-1.dat")
    para_list = parameter_reader(file_para)
    phi_list = mode_reader(file_phi)
    apara_list = mode_reader(file_apara)
    cases = []
    for i, case in enumerate(para_list):
        try:
            case['phi'] = phi_list[i]
            case['apara'] = apara_list[i]
        except IndexError:
            print('The number of cases is not equal to the number of modes.')
            break
        cases.append(case)
    return cases


def difference(case):
    th, phi_real, phi_imag = case['phi']
    _, psi_real, psi_imag = case['psi']

    diff_real = np.array(phi_real)-np.array(psi_real)
    diff_imag = np.array(phi_imag)-np.array(psi_imag)
    case['diff'] = [th, diff_real, diff_imag]
    # return [th, diff_real, diff_imag]
    return case


def mode_plot(case, axes):
    # ax1, ax2, ax3 = axes
    ax1, ax2 = axes
    # th, real, imag = case[name]
    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.5, 3))
    th, phi_real, phi_imag = case['phi']
    _, a_real, a_imag = case['apara']
    _, d_real, d_imag = case['diff']
    phi_real = np.array(phi_real)
    phi_imag = np.array(phi_imag)
    d_real = np.array(d_real)
    d_imag = np.array(d_imag)
    abs_phi = np.sqrt(phi_real**2 + phi_imag**2)
    abs_diff = np.sqrt(d_real**2 + d_imag**2)

    th_max = max(th)

    ax1.plot(th, phi_real, label='real')
    ax1.plot(th, phi_imag, label='imaginary')
    ax1.legend(prop={'size': 8})
    ax1.set_xlim(left=0, right=th_max)
    ax1.set_xlabel(r'$\theta$')
    ax1.set_ylabel(r"$\hat{\phi}$")
    ax2.plot(th, a_real, label='real')
    ax2.plot(th, a_imag, label='imaginary')
    ax2.legend(prop={'size': 8})
    ax2.set_xlim(left=0, right=th_max)
    ax2.set_xlabel(r'$\theta$')
    ax2.set_ylabel(r"$\hat{A}_\parallel$")
    # ax3.plot(th, abs_phi, label=r'$|\hat{\phi}|$')
    # ax3.plot(th, abs_diff, label=r'$|\hat{\phi}-\hat{\psi}|$')
    # ax3.legend(prop={'size': 8})
    # ax3.set_xlim(left=0, right=th_max)
    # ax3.set_xlabel(r'$\theta$')
    # ax3.set_ylabel(r"$|\hat{\phi}-\hat{\psi}|&|\hat{\phi}|$")
    ax2.set_title("{:+.5f}+{:+.5f}i".format(case['omr'], case['omi']))
    return


def mode_plot_with_psi(case, axes):
    ax1, ax2, ax3 = axes
    # ax1, ax2 = axes
    # th, real, imag = case[name]
    # fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=(13.5, 3))
    th, phi_real, phi_imag = case['phi']
    _, a_real, a_imag = case['apara']
    _, d_real, d_imag = case['diff']
    phi_real = np.array(phi_real)
    phi_imag = np.array(phi_imag)
    d_real = np.array(d_real)
    d_imag = np.array(d_imag)
    abs_phi = np.sqrt(phi_real**2 + phi_imag**2)
    abs_diff = np.sqrt(d_real**2 + d_imag**2)

    th_max = max(th)

    ax1.plot(th, phi_real, label='real')
    ax1.plot(th, phi_imag, label='imaginary')
    ax1.legend(prop={'size': 8})
    ax1.set_xlim(left=0, right=th_max)
    ax1.set_xlabel(r'$\theta$')
    ax1.set_ylabel(r"$\hat{\phi}$")
    ax2.plot(th, a_real, label='real')
    ax2.plot(th, a_imag, label='imaginary')
    ax2.legend(prop={'size': 8})
    ax2.set_xlim(left=0, right=th_max)
    ax2.set_xlabel(r'$\theta$')
    ax2.set_ylabel(r"$\hat{A}_\parallel$")
    ax3.plot(th, abs_phi, label=r'$|\hat{\phi}|$')
    ax3.plot(th, abs_diff, label=r'$|\hat{\phi}-\hat{\psi}|$')
    ax3.legend(prop={'size': 8})
    ax3.set_xlim(left=0, right=th_max)
    ax3.set_xlabel(r'$\theta$')
    ax3.set_ylabel(r"$|\hat{\phi}-\hat{\psi}|$"+r"\&"+r"$|\hat{\phi}|$")
    ax2.set_title("{:+.5f}+{:+.5f}i".format(case['omr'], case['omi']))
    return


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
