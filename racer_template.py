#!/usr/bin/env python3
import vectors
from math import sqrt, pi
import time
import argparse

_DEBUG = False

def load_inputs(input_path):

    inputs = []
    n_lm, n_ext, n_sg_ext, n_inputs = None, None, None, None
    with open(input_path,'r') as f:
        for i_l,l in enumerate(f.readlines()):
            if i_l == 0:
                n_lm, n_ext, n_sg_ext, n_inputs = [ int(i_str) for i_str in l.strip().split(' ') ]
                continue
            v_is = l.strip().replace('E','e').replace('D','e').replace('d','e').split(' ')
            inputs.append(
                (
                    [vectors.Vector([float(v_i) for v_i in v_is[i_chunk*4:(i_chunk+1)*4]]) for i_chunk in range((len(v_is)-1)//4)],
                    float(v_is[-1])
                )
            )
    return [n_lm, n_ext, n_sg_ext, inputs]

def numerator(em,spatial_dots):
    sd = spatial_dots
%(numerator)s

def evaluate_ltd_cut(i_ltd_cut, lm, extm, sgextm, em, em_osE, spatial_dots):
    sd = spatial_dots
%(evaluate_ltd_cut)s

def set_em_energies(sglm, sgextm, em):

%(set_em_energies)s

def set_em(sglm, sgextm, em):

%(set_em)s


def race(n_lm, n_ext, n_sg_ext, input_data):

    # Naming conventions (all 4D vectors):
    # lm : loop momenta of the amplitude
    # sglm : loop momenta of the supergraph (input from rust)
    # extm : external momenta from the amplitude (but not containing external momenta of the supergraph)
    # sgextm: external momenta from the supergraph (each appears one, not twice for each side of the forward scattering graph)
    # em: all edge momenta of the supergraph
    # spatial_dots: all spatial dots between edge momenta

    # SG loop momenta
    sglm = input_data[:n_lm+n_ext]

    # Amplitude loop momenta
    lm = input_data[:n_lm]
    if _DEBUG:
        print("Amplitude loop momenta:")
        print("-----------------------")
        for i_lm, m in enumerate(lm):
            print('{:<5s} : {:s}'.format('k{:d}'.format(i_lm),','.join('{:>+30.16e}'.format(m_i) for m_i in m)))

    # Amplitude external momenta
    extm = input_data[n_lm:n_lm+n_ext]
    if _DEBUG:
        print("Amplitude external momenta:")
        print("---------------------------")
        for i_m, m in enumerate(extm):
            print('{:<5s} : {:s}'.format('pl{:d}'.format(i_m),','.join('{:>+30.16e}'.format(m_i) for m_i in m)))
    # Supergraph externals
    sgextm = input_data[n_lm+n_ext:]
    if _DEBUG:
        print("Supergraph external momenta:")
        print("-----------------------")
        for i_m, m in enumerate(sgextm):
            print('{:<5s} : {:s}'.format('p{:d}'.format(i_m),','.join('{:>+30.16e}'.format(m_i) for m_i in m)))
    # Amplitude edges momenta
    em = [None,]*%(n_amplitude_edges)s
    set_em(sglm, sgextm, em)
    if _DEBUG:
        mom_id_to_name_map = %(mom_id_to_name_map)s
        print("All edges momenta:")
        print("------------------")
        for i_em, m in enumerate(em):
            print('{:<5s} : {:s}'.format(mom_id_to_name_map[i_em],','.join('{:>+30.16e}'.format(m_i) for m_i in m)))
    em_osE = [None,]*%(n_amplitude_edges)s
    for i_e, m in enumerate(em):
        em_osE[i_e] = sqrt(m[1]*m[1]+m[2]*m[2]+m[3]*m[3])
    spatial_dots = []
    for i, m_i in enumerate(em):
        spatial_dots.append( [ m_i[1]*m_j[1]+m_i[2]*m_j[2]+m_i[3]*m_j[3] for j, m_j in enumerate(em) ] )

    tot = 0.
    for i_ltd_cut in range(%(n_ltd_cuts)d):
        res = evaluate_ltd_cut(i_ltd_cut, lm, extm, sgextm, em, em_osE, spatial_dots)
        if _DEBUG: print("Result from LTD cut #{:d}: {:.16e}".format(i_ltd_cut, res))
        tot += res
    return tot*(%(normalisation)s)

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Python racer.""")

    parser.add_argument('--n_runs', '-n_runs', dest='n_runs', type=int, default=1,
        help='Specify number of runs to average over (default: %%(default)s).')
    parser.add_argument('--debug', '-d', dest='debug', action='store_true', default=False,
        help='Enable debug mode (default: %%(default)s).')
    args = parser.parse_args()

    _DEBUG = args.debug

    n_lm, n_ext, n_sg_ext, inputs  = load_inputs('./%(input_file_name)s')
    
    if _DEBUG:
        n_runs = 1
    else:
        n_runs = args.n_runs

    all_res = []
    all_targets = []
    t_start = time.time()
    for i_r in range(n_runs):
        for i, (input_data, target) in enumerate(inputs):
            res = race(n_lm, n_ext, n_sg_ext, input_data)
            all_res.append(res)
            all_targets.append(target)
            if _DEBUG:
                header="Sample #{:d}:".format(i)
                print("-"*len(header))
                print(header)
                print("-"*len(header))
                print("Result   : {:.16e}".format(res))
                print("Target   : {:.16e}".format(target))
                print("Rel diff : {:.16e}".format(abs((res-target)/res)))
                print("Ratio    : {:.16e}".format(res/target))
    tot_time = time.time()-t_start
    print('{:.16e}'.format(float(tot_time/(n_runs*len(inputs)))*1000000.))
    max_diff = -1.
    max_index = None
    for i, (res, target) in enumerate(zip(all_res,all_targets)):
        if abs((res-target)/target)>max_diff:
            max_diff = abs((res-target)/target)
            max_index = i
    print('{:.5e}'.format(max_diff*100.))
    print("Max difference of {:.3e}%% for input #{:d} with res={:.16e} and target_res={:.16e}".format(
        max_diff*100., max_index, all_res[max_index], all_targets[max_index]
    ))