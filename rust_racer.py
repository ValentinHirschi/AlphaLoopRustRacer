#!/usr/bin/env python3
from vectors import Vector, LorentzVector
import sys
import os
from pprint import pprint, pformat
pjoin = os.path.join

import time
import yaml
import random
import argparse
import numpy as np
import itertools
import subprocess

root_path = os.path.dirname(os.path.realpath( __file__ ))

class bcolors:
    HEADER = '\033[95m'
    BLUE = '\033[94m'
    GREEN = '\033[92m'
    YELLOW = '\033[93m'
    WARNING = YELLOW
    FAIL = '\033[91m'
    ENDC = '\033[0m'
    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'
    PURPLE = '\033[95m'
    CYAN = '\033[96m'
    DARKCYAN = '\033[36m'
    RED = '\033[91m'
    END = ENDC

class RacerError(Exception):
    pass

# Below only needs to be specified when generating input data
_AL_PATH = os.path.abspath(pjoin(root_path,os.pardir))
_MG_PATH = os.path.abspath(pjoin(root_path,os.pardir,os.pardir,os.pardir))

_PARAMETERS = {
    'ge' : 0.30795376724436885,
    'gs' : 1.2177157847767197,
    'h1' : 0.48380377501264738928 # h-function evaluated at 1, i.e. h(1)=1/(E^((t^2 + 1)/t)*N[2*BesselK[1, 2], 20]) /. {t -> 1}
}

_INCLUDE_NUMERATOR = True

import logging
logger = logging.getLogger("RustRacer")
logger.setLevel(logging.INFO)
FORMAT = '%(asctime)s %(message)s'
logging.basicConfig(format=FORMAT)

class Racer(object):
    
    def __init__(self, recycle_inputs, n_inputs, application):

        self._APPLICATION = application
        self.set_process()
        self.recycle_inputs = recycle_inputs
        self._N_INPUTS = n_inputs
        self.sg_info = yaml.load(open(pjoin(self._PROCESS_PATH,'Rust_inputs','%s.yaml'%self._SG_ID),'r'), Loader=yaml.Loader)
        self.hyperparams = yaml.load(open(pjoin(root_path,'racer_hyperparameters.yaml'),'r'), Loader=yaml.Loader)
        self.externals = [ LorentzVector(v) for v in self.hyperparams['CrossSection']['incoming_momenta'] ]
        self.Q0 = sum(v[0] for v in self.externals)
        self.cut_id = None
        for i_cc, cc_info in enumerate(self.sg_info['cutkosky_cuts']):
            if set(c['name'] for c in cc_info['cuts']) == self._CUT_TO_RACE:
                self.cut_id = i_cc
                break
        
        self.lmb = self.sg_info['loop_momentum_basis']
        self.n_loops = len(self.lmb)
        self.edge_signatures = self.sg_info['edge_signatures']
        self.sorted_edges = sorted(list(self.edge_signatures.keys()))
        self.sorted_edges_map = {k : i_k for i_k, k in enumerate(self.sorted_edges)}
        if self.cut_id is None:
            raise RacerError("Could not find Cutkosky cut (%s)."%(','.join(sorted(list(self._CUT_TO_RACE)))))
        

        self.cut_info = self.sg_info['cutkosky_cuts'][self.cut_id]
        self.amp_info = None
        if len(self.cut_info['diagram_sets'])!=1:
            raise RacerError("Can only race cuts with exactly one diagram set.")
        if len(self.cut_info['diagram_sets'][0]['diagram_info'])!=2:
            raise RacerError("Can only race cuts with exactly two diagrams in the diagram set.")
        for diag_info in self.cut_info['diagram_sets'][0]['diagram_info']:
            if not diag_info['conjugate_deformation']:
                self.amp_info = diag_info['graph']
                break
        if self.amp_info is None:
            raise RacerError("Could not find amplitude with no complex conjugation in diagram info.")

        logger.info("Will race cut #%d (%s) of %s%s%s over %s%d%s inputs"%(
            self.cut_id, ', '.join(sorted(list(self._CUT_TO_RACE))),bcolors.GREEN,self._APPLICATION,bcolors.END,bcolors.BLUE,self._N_INPUTS,bcolors.END)
        )

    def set_process(self):
        if self._APPLICATION == 'a_ddx_NLO_SG_QG0':

            self._SG_ID = 'SG_QG0'
            self._PROCESS_PATH = pjoin(_MG_PATH,'LU_a_ddx_NLO_SG_QG0_bare')
            self._CUT_TO_RACE = {'pq4','pq3'}
            if _INCLUDE_NUMERATOR:
                self._RAW_NUMERATOR = """ -32*(em[pq1][0]*em[pq4][0] - sd[pq1][pq4])*(em[pq2][0]*em[pq3][0] - sd[pq2][pq3]) """
                self._NORMALISATION = """ %(h1)s*%(ge)s**2*%(gs)s**2/((2.0e0*pi)**5*9.0e0)"""
            else:
                self._RAW_NUMERATOR = """ 1 """
                self._NORMALISATION = """ 1 """
            #############################################################################################

        elif self._APPLICATION == 'a_ddx_NNLO_SG_QG3':

            self._SG_ID = 'SG_QG3'
            self._PROCESS_PATH = pjoin(_MG_PATH,'LU_a_ddx_NNLO_SG_QG3_bare')
            self._CUT_TO_RACE = {'pq4','pq3'}
            if _INCLUDE_NUMERATOR:
                self._RAW_NUMERATOR = """ 32*((em[pq1][0]*em[pq4][0] - sd[pq1][pq4])*(em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq3][0]*em[pq6][0] - sd[pq3][pq6]) - 
        (em[pq1][0]*em[pq4][0] - sd[pq1][pq4])*(em[pq2][0]*em[pq6][0] - sd[pq2][pq6])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7]) - 
        (em[pq1][0]*em[pq7][0] - sd[pq1][pq7])*((em[pq2][0]*em[pq6][0] - sd[pq2][pq6])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4]) - 
            (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*(em[pq3][0]*em[pq6][0] - sd[pq3][pq6]) + (em[pq2][0]*em[pq3][0] - sd[pq2][pq3])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6])) + 
        (em[pq1][0]*em[pq3][0] - sd[pq1][pq3])*(em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) - 
        (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) + 
        (em[pq1][0]*em[pq6][0] - sd[pq1][pq6])*(-((em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4])) + 
            (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7]) + (em[pq2][0]*em[pq3][0] - sd[pq2][pq3])*(em[pq4][0]*em[pq7][0] - sd[pq4][pq7])) + 
        (em[pq1][0]*em[pq3][0] - sd[pq1][pq3])*(em[pq2][0]*em[pq6][0] - sd[pq2][pq6])*(em[pq4][0]*em[pq7][0] - sd[pq4][pq7]) - 
        (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq3][0]*em[pq6][0] - sd[pq3][pq6])*(em[pq4][0]*em[pq7][0] - sd[pq4][pq7]) - 
        ((em[pq1][0]*em[pq4][0] - sd[pq1][pq4])*(em[pq2][0]*em[pq3][0] - sd[pq2][pq3]) + (em[pq1][0]*em[pq3][0] - sd[pq1][pq3])*(em[pq2][0]*em[pq4][0] - sd[pq2][pq4]) - 
            (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4]))*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7])) """
                self._NORMALISATION = """ %(h1)s*%(ge)s**2*%(gs)s**4/((2*pi)**5*(9.0e0))/(-(3.0e0/4.0e0)*(2.0e0*pi)**3)"""
            else:
                self._RAW_NUMERATOR = """ 1 """
                self._NORMALISATION = """ 1 """

        elif self._APPLICATION == 'a_ddx_N3LO_SG_QG108':

            self._SG_ID = "SG_QG108"
            self._PROCESS_PATH = pjoin(_MG_PATH,'LU_a_ddx_N3LO_SG_QG108_bare')
            self._CUT_TO_RACE = {'pq4','pq3'}
            if _INCLUDE_NUMERATOR:
                self._RAW_NUMERATOR = """-128*((em[pq1][0]*em[pq6][0] - sd[pq1][pq6])*((em[pq10][0]*em[pq9][0] - sd[pq10][pq9])*((em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4]) - 
                    (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7])) + (em[pq10][0]*em[pq7][0] - sd[pq10][pq7])*
                    (-((em[pq2][0]*em[pq9][0] - sd[pq2][pq9])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4])) + (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*
                    (em[pq3][0]*em[pq9][0] - sd[pq3][pq9])) + (em[pq10][0]*em[pq4][0] - sd[pq10][pq4])*
                    ((em[pq2][0]*em[pq9][0] - sd[pq2][pq9])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7]) - (em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq3][0]*em[pq9][0] - sd[pq3][pq9]))) + 
                (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq10][0]*em[pq9][0] - sd[pq10][pq9])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) - 
                (em[pq1][0]*em[pq10][0] - sd[pq1][pq10])*(em[pq2][0]*em[pq9][0] - sd[pq2][pq9])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) - 
                (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq10][0]*em[pq7][0] - sd[pq10][pq7])*(em[pq3][0]*em[pq9][0] - sd[pq3][pq9])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) + 
                (em[pq1][0]*em[pq10][0] - sd[pq1][pq10])*(em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq3][0]*em[pq9][0] - sd[pq3][pq9])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) - 
                (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq10][0]*em[pq9][0] - sd[pq10][pq9])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4])*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7]) + 
                (em[pq1][0]*em[pq10][0] - sd[pq1][pq10])*(em[pq2][0]*em[pq9][0] - sd[pq2][pq9])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4])*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7]) + 
                (em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*(em[pq10][0]*em[pq4][0] - sd[pq10][pq4])*(em[pq3][0]*em[pq9][0] - sd[pq3][pq9])*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7]) - 
                (em[pq1][0]*em[pq10][0] - sd[pq1][pq10])*(em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*(em[pq3][0]*em[pq9][0] - sd[pq3][pq9])*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7]) + 
                (em[pq1][0]*em[pq3][0] - sd[pq1][pq3])*((em[pq10][0]*em[pq9][0] - sd[pq10][pq9])*(-((em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6])) + 
                    (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7])) + (em[pq10][0]*em[pq7][0] - sd[pq10][pq7])*
                    ((em[pq2][0]*em[pq9][0] - sd[pq2][pq9])*(em[pq4][0]*em[pq6][0] - sd[pq4][pq6]) - (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*(em[pq6][0]*em[pq9][0] - sd[pq6][pq9])) + 
                    (em[pq10][0]*em[pq4][0] - sd[pq10][pq4])*(-((em[pq2][0]*em[pq9][0] - sd[pq2][pq9])*(em[pq6][0]*em[pq7][0] - sd[pq6][pq7])) + 
                    (em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq6][0]*em[pq9][0] - sd[pq6][pq9]))) + 
                ((em[pq1][0]*em[pq2][0] - sd[pq1][pq2])*((em[pq10][0]*em[pq7][0] - sd[pq10][pq7])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4]) - 
                    (em[pq10][0]*em[pq4][0] - sd[pq10][pq4])*(em[pq3][0]*em[pq7][0] - sd[pq3][pq7])) + (em[pq1][0]*em[pq10][0] - sd[pq1][pq10])*
                    (-((em[pq2][0]*em[pq7][0] - sd[pq2][pq7])*(em[pq3][0]*em[pq4][0] - sd[pq3][pq4])) + (em[pq2][0]*em[pq4][0] - sd[pq2][pq4])*
                    (em[pq3][0]*em[pq7][0] - sd[pq3][pq7])))*(em[pq6][0]*em[pq9][0] - sd[pq6][pq9]) - 
                ((em[pq10][0]*em[pq6][0] - sd[pq10][pq6])*(em[pq2][0]*em[pq3][0] - sd[pq2][pq3]) + (em[pq10][0]*em[pq3][0] - sd[pq10][pq3])*(em[pq2][0]*em[pq6][0] - sd[pq2][pq6]) - 
                    (em[pq10][0]*em[pq2][0] - sd[pq10][pq2])*(em[pq3][0]*em[pq6][0] - sd[pq3][pq6]))*((em[pq1][0]*em[pq9][0] - sd[pq1][pq9])*(em[pq4][0]*em[pq7][0] - sd[pq4][pq7]) - 
                    (em[pq1][0]*em[pq7][0] - sd[pq1][pq7])*(em[pq4][0]*em[pq9][0] - sd[pq4][pq9]) - (em[pq1][0]*em[pq4][0] - sd[pq1][pq4])*(em[pq7][0]*em[pq9][0] - sd[pq7][pq9])))"""
                self._NORMALISATION = """ %(h1)s*%(ge)s**2*%(gs)s**6/((2.0e0*pi)**5*(9.0e0))/(36.0e0*pi**6) """
            else:
                self._RAW_NUMERATOR = """ 1 """
                self._NORMALISATION = """ 1 """
            #############################################################################################

        else:
            raise RacerError("Application '%s' not reckognised."%self._APPLICATION)

    def generate_python_racer(self):

        logger.info("Generating Python implementation...")

        repl_dict = {}
        repl_dict['input_file_name'] = 'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS)

        repl_dict['n_amplitude_edges'] = len(self.sorted_edges)

        repl_dict['set_em'] = [
            'em[%d]=(%s%s)'%(i_e,''.join(
                '%ssglm[%d]'%('+' if sig > 0 else '-', i_sig) for i_sig, sig in enumerate(self.edge_signatures[e][0]) if sig!=0
            ),''.join(
                '%ssgextm[%d]'%('+' if sig > 0 else '-', (i_sig%(len(self.edge_signatures[e][1])//2))) for i_sig, sig in enumerate(self.edge_signatures[e][1]) if sig!=0
            ))
            for i_e, e in enumerate(self.sorted_edges)
        ]
        repl_dict['set_em_energies'] = [
            'em[%d][0]=(%s%s)'%(i_e,''.join(
                '%ssglm[%d][0]'%('+' if sig > 0 else '-', i_sig) for i_sig, sig in enumerate(self.edge_signatures[e][0]) if sig!=0
            ),''.join(
                '%ssgextm[%d][0]'%('+' if sig > 0 else '-', (i_sig%(len(self.edge_signatures[e][1])//2))) for i_sig, sig in enumerate(self.edge_signatures[e][1]) if sig!=0
            ))
            for i_e, e in enumerate(self.sorted_edges)
        ]

        repl_dict['mom_id_to_name_map'] = str({v:k for k,v in self.sorted_edges_map.items()})

        repl_dict['evaluate_ltd_cut'] = []
        all_ltd_cuts = []
        for i_cs, cs in enumerate(self.amp_info['ltd_cut_structure']):
            loop_lines_cut = [ [ (i_ll, i_prop, s, i_cs, prop['name']) for i_prop, prop in enumerate(self.amp_info['loop_lines'][i_ll]['propagators']) ] 
                                for i_ll, s in [(i, s_i) for i, s_i in enumerate(cs) if s_i!=0]]
            for cuts in itertools.product(*loop_lines_cut):
                all_ltd_cuts.append(cuts)
        repl_dict['n_ltd_cuts'] = len(all_ltd_cuts)

        for i_ltd_cut, ltd_cuts in enumerate(all_ltd_cuts):
            repl_dict['evaluate_ltd_cut'].append("# LTD cut #%d, cut structure #%d: %s"%(i_ltd_cut, ltd_cuts[0][3], ','.join(p[-1] for p in ltd_cuts)))
            repl_dict['evaluate_ltd_cut'].append('%sif i_ltd_cut==%d:'%('el' if i_ltd_cut!=0 else '',i_ltd_cut))
            repl_dict['evaluate_ltd_cut'].append('    # Set on-shell energies of cut momenta')
            for (i_ll, i_prop, cut_sign, i_cs, edge_name) in ltd_cuts:
                # Caching is better here
                #repl_dict['evaluate_ltd_cut'].append('    em[%(en)d][0]=sqrt(em[%(en)d][1]*em[%(en)d][1]+em[%(en)d][2]*em[%(en)d][2]+em[%(en)d][3]*em[%(en)d][3])'%{'en':self.sorted_edges_map[edge_name]})
                repl_dict['evaluate_ltd_cut'].append('    em[%(en)d][0]=%(cs)sem_osE[%(en)d]'%{'cs': '+' if cut_sign>0 else '-','en':self.sorted_edges_map[edge_name]})
            repl_dict['evaluate_ltd_cut'].append('    # Set on-shell energies in LMB')
            # construct the cut basis to LTD loop momentum basis mapping
            mat = [self.amp_info['loop_lines'][i_ll]['signature'] for (i_ll, i_prop, cut_sign, i_cs, name) in ltd_cuts]
            transfo = np.linalg.inv(np.array(mat))
            if any(any(abs(sgn) not in [-1,0,1] for sgn in l_transfo) for l_transfo in transfo):
                raise RacerError("Unary basis transformation matrix has elements not in -1, 0, 1: %s"%str(transfo))

            shifts = []
            for (i_ll, i_prop, cut_sign, i_cs, edge_name) in ltd_cuts:
                # Note: change of sign intentional!
                shift = ''.join('%sextm[%d][0]'%('-' if sgn>0 else '+', i_lm) for i_lm, sgn in enumerate(self.amp_info['loop_lines'][i_ll]['propagators'][i_prop]['parametric_shift'][0]) if sgn!=0)
                shift += ''.join('%ssgextm[%d][0]'%('-' if sgn>0 else '+', i_lm) for i_lm, sgn in enumerate(self.amp_info['loop_lines'][i_ll]['propagators'][i_prop]['parametric_shift'][1]) if sgn!=0)
                if shift!='':
                    shifts.append(shift)
                else:
                    shifts.append(None)
            for i_lm, row in enumerate(transfo):                    
                repl_dict['evaluate_ltd_cut'].append('    lm[%d][0]=%s'%(
                    i_lm, ''.join('%s%s'%('+' if sgn>0 else '-', 'em[%d][0]'%self.sorted_edges_map[ltd_cuts[i_col][-1]] if shifts[i_col] is None else 
                                    '(em[%d][0]%s)'%(self.sorted_edges_map[ltd_cuts[i_col][-1]], shifts[i_col])) for i_col, sgn in enumerate(row) if sgn!=0)
                ))
            repl_dict['evaluate_ltd_cut'].append('    # Set all_edge energies')
            repl_dict['evaluate_ltd_cut'].append('    sglm = lm+[%s]'%(','.join('em[%d]'%self.sorted_edges_map[e_name] for e_name in self.lmb[self.amp_info['n_loops']:])))
            repl_dict['evaluate_ltd_cut'].append('    set_em_energies(sglm, sgextm, em)')
            repl_dict['evaluate_ltd_cut'].append('    # Compute denominator')
            
            ltd_cut_signs = {c[-1]:c[2] for c in  ltd_cuts}
            repl_dict['evaluate_ltd_cut'].append('    denom = %s'%('*'.join( 
                ('(em[%(i)d][0]*em[%(i)d][0]-sd[%(i)d][%(i)d])'%{'i':self.sorted_edges_map[prop['name']]}  if prop['name'] not in ltd_cut_signs else 
                                            '(%s2*em[%d][0])'%('+' if ltd_cut_signs[prop['name']]>0 else '-',self.sorted_edges_map[prop['name']]))
                for ll in self.amp_info['loop_lines'] for prop in ll['propagators']
            )))
            repl_dict['evaluate_ltd_cut'].append('    # Compute numerator')
            repl_dict['evaluate_ltd_cut'].append('    num = numerator(em,sd)')
            repl_dict['evaluate_ltd_cut'].append('    return num/denom')

        repl_dict['evaluate_ltd_cut'].append('else:')
        repl_dict['evaluate_ltd_cut'].append('    print("Cut {:d} not found.".format(i_ltd_cut))')

        # numerator
        numerator = self._RAW_NUMERATOR
        for i_edge, e_name in enumerate(self.sorted_edges):
            numerator = numerator.replace('[%s]'%e_name,'[%d]'%self.sorted_edges_map[e_name])
        repl_dict['numerator'] = ['return (%s)'%(numerator),]

        python_parameter = {k:('%.16e'%v) for k,v in _PARAMETERS.items()}
        repl_dict['normalisation'] = self._NORMALISATION%python_parameter

        template = open(pjoin(root_path, 'racer_template.py'),'r').read()
        for k in ['set_em','set_em_energies','evaluate_ltd_cut','numerator']:
            repl_dict[k] = '\n'.join('    %s'%l for l in repl_dict[k])
        with open(pjoin(root_path,'racer.py'),'w') as f:
            f.write(template%repl_dict)

    def get_t_scaling(self, moms):

        scaling_solutions = list(self.rust_worker.get_scaling(moms,self.cut_id))
        scaling, scaling_jacobian = scaling_solutions.pop(0)
        while scaling < 0.0:
            if len(scaling_solutions)==0:
                break
            scaling, scaling_jacobian = scaling_solutions.pop(0)

        return scaling, scaling_jacobian

    def generate_fortran_racer(self):

        logger.info("Generating Fortran implementation...")

        repl_dict = {}

        repl_dict['input_file_name'] = 'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS)
        repl_dict['n_lm'] = self.amp_info['n_loops']
        repl_dict['n_sglm'] = self.n_loops
        repl_dict['n_extm'] = self.n_loops-self.amp_info['n_loops']
        repl_dict['n_sgextm'] = len(self.externals)
        repl_dict['n_edges'] =  len(self.sorted_edges)
        repl_dict['n_inputs'] = self._N_INPUTS

        repl_dict['warmup_code'] = []
        repl_dict['warmup_code'].append("! Compute edge momenta")
        repl_dict['warmup_code'].extend([
            'em(%d,:)=(%s%s)'%(i_e,''.join(
                '%ssglm(%d,:)'%('+' if sig > 0 else '-', i_sig) for i_sig, sig in enumerate(self.edge_signatures[e][0]) if sig!=0
            ),''.join(
                '%ssgextm(%d,:)'%('+' if sig > 0 else '-', (i_sig%(len(self.edge_signatures[e][1])//2))) for i_sig, sig in enumerate(self.edge_signatures[e][1]) if sig!=0
            ))
            for i_e, e in enumerate(self.sorted_edges)
        ])
        repl_dict['warmup_code'].append("! Compute on-shell energies")
        
        repl_dict['warmup_code'].append('do I=0,%d-1'%len(self.sorted_edges))
        repl_dict['warmup_code'].append('    do J=0,%d-1'%len(self.sorted_edges))
        repl_dict['warmup_code'].append('        sd(I,J) = em(I,1)*em(J,1)+em(I,2)*em(J,2)+em(I,3)*em(J,3)')
        repl_dict['warmup_code'].append('    enddo')
        repl_dict['warmup_code'].append('enddo')
        repl_dict['warmup_code'].extend([
            'em_osE(%d)=dsqrt(sd(%d,%d))'%(i_e,i_e,i_e)
            for i_e, e in enumerate(self.sorted_edges)
        ])

        repl_dict['evaluate_ltd_cut'] = []

        all_ltd_cuts = []
        for i_cs, cs in enumerate(self.amp_info['ltd_cut_structure']):
            loop_lines_cut = [ [ (i_ll, i_prop, s, i_cs, prop['name']) for i_prop, prop in enumerate(self.amp_info['loop_lines'][i_ll]['propagators']) ] 
                                for i_ll, s in [(i, s_i) for i, s_i in enumerate(cs) if s_i!=0]]
            for cuts in itertools.product(*loop_lines_cut):
                all_ltd_cuts.append(cuts)

        for i_ltd_cut, ltd_cuts in enumerate(all_ltd_cuts):
            repl_dict['evaluate_ltd_cut'].append('')
            repl_dict['evaluate_ltd_cut'].append("! LTD cut #%d, cut structure #%d: %s"%(i_ltd_cut, ltd_cuts[0][3], ','.join(p[-1] for p in ltd_cuts)))
            repl_dict['evaluate_ltd_cut'].append('! Set on-shell energies of cut momenta')
            for (i_ll, i_prop, cut_sign, i_cs, edge_name) in ltd_cuts:
                # Caching is better here
                #repl_dict['evaluate_ltd_cut'].append('    em[%(en)d][0]=sqrt(em[%(en)d][1]*em[%(en)d][1]+em[%(en)d][2]*em[%(en)d][2]+em[%(en)d][3]*em[%(en)d][3])'%{'en':self.sorted_edges_map[edge_name]})
                repl_dict['evaluate_ltd_cut'].append('em(%(en)d,0)=%(cs)sem_osE(%(en)d)'%{'cs': '+' if cut_sign>0 else '-','en':self.sorted_edges_map[edge_name]})
            repl_dict['evaluate_ltd_cut'].append('! Set on-shell energies in LMB')
            # construct the cut basis to LTD loop momentum basis mapping
            mat = [self.amp_info['loop_lines'][i_ll]['signature'] for (i_ll, i_prop, cut_sign, i_cs, name) in ltd_cuts]
            transfo = np.linalg.inv(np.array(mat))
            if any(any(abs(sgn) not in [-1,0,1] for sgn in l_transfo) for l_transfo in transfo):
                raise RacerError("Unary basis transformation matrix has elements not in -1, 0, 1: %s"%str(transfo))

            shifts = []
            for (i_ll, i_prop, cut_sign, i_cs, edge_name) in ltd_cuts:
                # Note: change of sign intentional!
                shift = ''.join('%sextm(%d,0)'%('-' if sgn>0 else '+', i_lm) for i_lm, sgn in enumerate(self.amp_info['loop_lines'][i_ll]['propagators'][i_prop]['parametric_shift'][0]) if sgn!=0)
                shift += ''.join('%ssgextm(%d,0)'%('-' if sgn>0 else '+', i_lm) for i_lm, sgn in enumerate(self.amp_info['loop_lines'][i_ll]['propagators'][i_prop]['parametric_shift'][1]) if sgn!=0)
                if shift!='':
                    shifts.append(shift)
                else:
                    shifts.append(None)
            for i_lm, row in enumerate(transfo):                    
                repl_dict['evaluate_ltd_cut'].append('lm(%d,0)=%s'%(
                    i_lm, ''.join('%s%s'%('+' if sgn>0 else '-', 'em(%d,0)'%self.sorted_edges_map[ltd_cuts[i_col][-1]] if shifts[i_col] is None else 
                                    '(em(%d,0)%s)'%(self.sorted_edges_map[ltd_cuts[i_col][-1]], shifts[i_col])) for i_col, sgn in enumerate(row) if sgn!=0)
                ))
            repl_dict['evaluate_ltd_cut'].append('! Set all_edge energies')
            repl_dict['evaluate_ltd_cut'].extend([
                'em(%d,0)=(%s%s)'%(i_e,''.join(
                    '%s%s(%d,0)'%('+' if sig > 0 else '-', 'sglm' if i_sig>=self.amp_info['n_loops'] else 'lm', i_sig) for i_sig, sig in enumerate(self.edge_signatures[e][0]) if sig!=0
                ),''.join(
                    '%ssgextm(%d,0)'%('+' if sig > 0 else '-', (i_sig%(len(self.edge_signatures[e][1])//2))) for i_sig, sig in enumerate(self.edge_signatures[e][1]) if sig!=0
                ))
                for i_e, e in enumerate(self.sorted_edges)
            ])
            repl_dict['evaluate_ltd_cut'].append('! Compute denominator')
            
            ltd_cut_signs = {c[-1]:c[2] for c in  ltd_cuts}
            repl_dict['evaluate_ltd_cut'].append('denom = %s'%('*'.join( 
                ('(em(%(i)d,0)*em(%(i)d,0)-sd(%(i)d,%(i)d))'%{'i':self.sorted_edges_map[prop['name']]}  if prop['name'] not in ltd_cut_signs else 
                                            '(%s2*em(%d,0))'%('+' if ltd_cut_signs[prop['name']]>0 else '-',self.sorted_edges_map[prop['name']]))
                for ll in self.amp_info['loop_lines'] for prop in ll['propagators']
            )))
            repl_dict['evaluate_ltd_cut'].append('! Compute numerator')
            # numerator
            numerator = self._RAW_NUMERATOR
            for i_edge, e_name in enumerate(self.sorted_edges):
                numerator = numerator.replace('[%s]'%e_name,'[%d]'%self.sorted_edges_map[e_name])
            for i_edge in range(len(self.sorted_edges)):
                for j_edge in range(len(self.sorted_edges)):
                    numerator = numerator.replace('sd[%d][%d]'%(i_edge,j_edge),'sd(%d,%d)'%tuple(sorted([i_edge,j_edge])))
                    numerator = numerator.replace('em[%d][%d]'%(i_edge,j_edge),'em(%d,%d)'%(i_edge,j_edge))

            numerator = ''.join(nl.strip() for nl in numerator.split('\n'))
            repl_dict['evaluate_ltd_cut'].append('num = %s'%numerator)
            repl_dict['evaluate_ltd_cut'].append('IF (DBG) THEN')
            repl_dict['evaluate_ltd_cut'].append('    write(*,*) "Result for LTD cut #%d, cut structure #%d, %s :", num/denom'%(i_ltd_cut, ltd_cuts[0][3], ','.join(p[-1] for p in ltd_cuts)))
            repl_dict['evaluate_ltd_cut'].append('ENDIF')
            repl_dict['evaluate_ltd_cut'].append('res = res + num/denom')

        repl_dict['wrapup_code'] = []
        fortran_parameter = {k:('%.16e'%v).replace('e','D') for k,v in _PARAMETERS.items()}
        repl_dict['wrapup_code'].append('res = res*(%s)'%(self._NORMALISATION%fortran_parameter).replace('e','D'))

        template = open(pjoin(root_path, 'racer_template.f'),'r').read()
        for k in ['warmup_code','evaluate_ltd_cut','wrapup_code']:
            repl_dict[k] = '\n'.join('%s       %s'%( (' ',l) if not l.startswith('!') else ('!',l[1:]) ) for l in repl_dict[k])
        with open(pjoin(root_path,'racer.f'),'w') as f:
            f.write(template%repl_dict)

    def generate_racer(self):

        # Sanity check
        if self.amp_info['loop_momentum_map'] != [ 
            [ [ 0 if j!=i else 1 for j in range(self.n_loops)], [0,]*(self.n_loops-self.amp_info['n_loops']+len(self.externals))]
            for i in range(self.amp_info['n_loops'])]:
                raise RacerError("The supergraph must be generated with an LMB that directly matches the cut momentum basis.")

        for i_cut, cut in enumerate(self.cut_info['cuts'][:-1]):
            if not all(s==0 for i_s, s in enumerate(cut['signature'][0]) if i_s != i_cut+self.amp_info['n_loops']) or cut['signature'][0][i_cut+self.amp_info['n_loops']]!=1:
                raise RacerError("Unsupported signature for cut edge %s : %s"%(cut['name'],str(cut['signature'])))
            if not all(s==0 for i_s, s in enumerate(cut['signature'][1])):
                raise RacerError("Unsupported signature for cut edge %s : %s"%(cut['name'],str(cut['signature'])))
            if abs(cut['particle_id']) not in range(1,7) and abs(cut['particle_id']) not in range(21,23):
                raise RacerError("Cut propagators must be massless.")
        
        # Generate input data for race
        if not os.path.exists(pjoin(root_path,'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS))) or not self.recycle_inputs:

            random.seed(1)

            if sys.path[0]!=_AL_PATH:
                sys.path.insert(0,_AL_PATH)
            import ltd
            #print(
            #    pjoin(self._PROCESS_PATH,'Rust_inputs','%s.yaml'%self._SG_ID),
            #    pjoin(root_path,'racer_hyperparameters.yaml')                
            #)
            self.rust_worker = ltd.CrossSection(
                pjoin(self._PROCESS_PATH,'Rust_inputs','%s.yaml'%self._SG_ID),
                pjoin(root_path,'racer_hyperparameters.yaml')
            )

            input_data = ['%d %d %d %d'%(self.amp_info['n_loops'], self.n_loops-self.amp_info['n_loops'], len(self.externals), self._N_INPUTS),]
            logger.info("Generating %d input sample points into file '%s'..."%(self._N_INPUTS,'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS)))
            for i_inputs in range(self._N_INPUTS):

                # inputs = [
                #     Vector([0.1,0.2,0.3]), # k == pq2
                #     Vector([0.4,0.5,0.6]), # l == pq7
                #     Vector([0.7,0.8,0.9]), # m == pq10
                #     Vector([1.0,1.1,1.2]), # n == pq4
                # ]
                inputs = [ Vector([random.random() for _ in range(3)]) for i_loop in range(self.n_loops) ]

                scaling, scaling_jacobian = self.get_t_scaling(inputs)
                # res = self.rust_worker.evaluate_cut(inputs,self.cut_id,scaling,scaling_jacobian)
                # res_f128 = self.rust_worker.evaluate_cut_f128(inputs,self.cut_id,scaling,scaling_jacobian)
                # print(scaling, scaling_jacobian,res,res_f128)
                
                rescaled_inputs = [ scaling*v for v in inputs]
                scaling, scaling_jacobian = self.get_t_scaling(rescaled_inputs)
                res = self.rust_worker.evaluate_cut(rescaled_inputs,self.cut_id,scaling,scaling_jacobian)
                res_f128 = self.rust_worker.evaluate_cut_f128(rescaled_inputs,self.cut_id,scaling,scaling_jacobian)
                # print(inputs,self.cut_id,scaling,scaling_jacobian)
                # print(res,res_f128)
                # print(scaling, scaling_jacobian,res,res_f128)

                rescaled_lms = [ LorentzVector([0.]+list(v)) for v in rescaled_inputs[:self.amp_info['n_loops']] ]
                rescaled_exts = [ LorentzVector([0.]+list(v)) for v in rescaled_inputs[self.amp_info['n_loops']:] ]
                rescaled_exts.append( sum( float(self.cut_info['cuts'][i_v]['sign']*self.cut_info['cuts'][-1]['sign'])*v for i_v, v in enumerate(rescaled_exts)) )

                # Set cut particles on-shell
                for i_v, v in enumerate(rescaled_exts):
                    v.set_square(0., negative=(self.cut_info['cuts'][i_v]['sign']<0))
                
                recomputed_Q0 = sum(self.cut_info['cuts'][i_v]['sign']*v[0] for i_v, v in enumerate(rescaled_exts))
                if abs(self.Q0 - recomputed_Q0)/self.Q0 > 1.0e-14:
                    raise RacerError("Error when computing cut energies: %.16e != %.16e."%(self.Q0, recomputed_Q0))

                # print("Amplitude loop momenta:")
                # print("-----------------------")
                # for i_lm, m in enumerate(rescaled_lms):
                #     print('{:<5s} : {:s} : {:>+30.16e}'.format('k{:d}'.format(i_lm),','.join('{:>+30.16e}'.format(m_i) for m_i in m),m.square()))

                # print("Amplitude external momenta:")
                # print("---------------------------")
                # for i_m, m in enumerate(rescaled_exts):
                #     print('{:<5s} : {:s} : {:>+30.16e}'.format('pl{:d}'.format(i_m),','.join('{:>+30.16e}'.format(m_i) for m_i in m),m.square()))

                input_data.append(
                    ('%s %s %s %.16e'%(
                        ' '.join(' '.join('%.16e'%v_i for v_i in v) for v in rescaled_lms),
                        ' '.join(' '.join('%.16e'%v_i for v_i in v) for v in rescaled_exts[:-1]),
                        ' '.join(' '.join('%.16e'%v_i for v_i in v) for v in self.externals),
                        res_f128[0]
                    )).replace('e','D')
                )

            with open(pjoin(root_path,'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS)),'w') as f:
                f.write('\n'.join(input_data))

        else:
            with open(pjoin(root_path,'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS)),'r') as f:
                new_n_inputs=int(f.readline().strip().split(' ')[-1])
                if new_n_inputs != self._N_INPUTS:
                    logger.warning("Recycling input file '%s' which contains a different number of inputs (%d) than requested."%('inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS),new_n_inputs))
                    self._N_INPUTS = new_n_inputs
            logger.info("Recycling existing %d inputs found in %s..."%(self._N_INPUTS,'inputs_%s_%d.txt'%(self._APPLICATION,self._N_INPUTS)))
        self.generate_python_racer()

        self.generate_fortran_racer()

    def race(self):

        n_runs = 1

        logger.info("Profiling Python implementation...")
        python_result = subprocess.check_output(' '.join([sys.executable,'racer.py','--n_runs','%d'%n_runs]),cwd=root_path,shell=True).decode(encoding='utf8').split('\n')
        logger.info("%-20s : %s%.6f \u00B5s%s / eval"%("Python timing (LTD)",bcolors.GREEN,float(python_result[0]),bcolors.END))
        logger.info("%-20s : %.2e%%"%("Python max diff.",float(python_result[1])))

        if os.path.exists(pjoin(self._PROCESS_PATH,'Rust_inputs','%s.yaml'%self._SG_ID)):

            if sys.path[0]!=_AL_PATH:
                sys.path.insert(0,_AL_PATH)
            import ltd
            self.rust_worker = ltd.CrossSection(
                pjoin(self._PROCESS_PATH,'Rust_inputs','%s.yaml'%self._SG_ID),
                pjoin(root_path,'racer_hyperparameters.yaml')
            )
            
            random.seed(1)
            rust_inputs = []
            for i_inputs in range(self._N_INPUTS):
                inputs = [ Vector([random.random() for _ in range(3)]) for i_loop in range(self.n_loops) ]
                scaling, scaling_jacobian = self.get_t_scaling(inputs)            
                rescaled_inputs = [ scaling*v for v in inputs]
                scaling, scaling_jacobian = self.get_t_scaling(rescaled_inputs)
                rust_inputs.append([rescaled_inputs,self.cut_id,scaling,scaling_jacobian])
        
            logger.info("Profiling Rust implementation...")
            t_start = time.time()
            for i_r in range(n_runs):
                for rust_input in rust_inputs:
                    self.rust_worker.evaluate_cut(*rust_input)
            rust_timing = (time.time()-t_start)/float(n_runs*len(rust_inputs))*1000000.
            logger.info("%-20s : %s%.6f \u00B5s%s / eval"%("Rust timing (cLTD)",bcolors.GREEN,rust_timing,bcolors.END))
        else:
            logger.info("Skipping profiling of rust because input file '%s' is not found."%pjoin(self._PROCESS_PATH,'Rust_inputs','%s.yaml'%self._SG_ID))

        logger.info("Compiling Fortran implementation...")
        subprocess.call(' '.join(['make','clean']),cwd=root_path,shell=True,stdout=subprocess.DEVNULL)
        subprocess.call(' '.join(['make','OPTLEVEL=3']),cwd=root_path,shell=True,stdout=subprocess.DEVNULL)
        if not os.path.exists(pjoin(root_path,'racer')):
            raise RacerError("Could not successfully compile Fortran racer.")
        logger.info("Profiling Fortran implementation...")
        fortran_result = subprocess.check_output(' '.join(['./racer',]),cwd=root_path,shell=True).decode(encoding='utf8').split('\n')
        logger.info("%-20s : %s%.6f \u00B5s%s / eval"%("Fortran timing (LTD)",bcolors.GREEN,float(fortran_result[0]),bcolors.END))
        logger.info("%-20s : %.2e%%"%("Fortran max diff.",float(fortran_result[1])))

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description="""Python racer generator""")

    parser.add_argument('--n_inputs', '-n_inputs', dest='n_inputs', type=int, default=1000,
        help='Specify number of n_inputs to perform timing over (default: %%(default)s).')
    parser.add_argument('--application', '-a', dest='application', type=str, default="a_ddx_NLO_SG_QG0",
        choices=['a_ddx_NLO_SG_QG0','a_ddx_NNLO_SG_QG3','a_ddx_N3LO_SG_QG108'],
        help='Supergraph cut contribution to race (default: %%(default)s).')
    parser.add_argument('--debug', '-d', dest='debug', action='store_true', default=False,
        help='Enable debug mode (default: %%(default)s).')
    parser.add_argument('--no_recycle_inputs', '-nr', dest='no_recycle_inputs', action='store_true', default=False,
        help='Force input generation even if inputs_<application>_<n_points>.txt is found (default: %%(default)s).')
    parser.add_argument('--generate_only', '-go', dest='generate_only', action='store_true', default=False,
        help='Only generate codes and not race (default: %%(default)s).')
    args = parser.parse_args()

    if args.debug:
        logger.setLevel(logging.DEBUG)

    racer = Racer((not args.no_recycle_inputs),args.n_inputs,args.application)
    racer.generate_racer()
    if not args.generate_only:
        racer.race()