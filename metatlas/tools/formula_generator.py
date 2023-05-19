from __future__ import absolute_import
from __future__ import print_function
from copy import deepcopy
from six.moves import range
def get_formulae(mass,tol=5,charge=0,cho_only=False,tol_type='ppm',max_tests=1e7,
    min_h=0,max_h=200,
    min_c=0,max_c=50,
    min_n=0,max_n=6,
    min_o=0,max_o=20,
    min_p=0,max_p=4,
    min_s=0,max_s=4,
    min_na=0,max_na=0,
    min_k=0,max_k=0,
    min_cl=0,max_cl=0):
    """
    performs brute force chemical formula generation.  enumerate all possible matches within bounds.
    Considers only [C,H,N,O,P,S,Na,K,Cl]

    returns list of tuples (error,formula).

    use it like this:
    from metatlas.helpers import formula_generator as formulae
    out = formulae.get_formulae(634.13226,1,min_c=20,min_h=15)

    (4.457399995771993e-06, 'H26C32O14')
    (9.535900062473956e-06, 'H29C34N5P1S3')
    (3.076309997140925e-05, 'H151C24N6O3P2')
    (5.4717500006518094e-05, 'H35C33N1O2P3S2')
    ....



    """

    elements = [
    {'num':0,'symbol':"H",'mass': 1.0078250321,'min':min_h,'max':max_h,'guess':0},
    {'num':1,'symbol':"C",'mass': 12.000000000,'min':min_c,'max':max_c,'guess':0},
    {'num':2,'symbol':"N",'mass': 14.003074005,'min':min_n,'max':max_n,'guess':0},
    {'num':3,'symbol':"O",'mass': 15.994914622,'min':min_o,'max':max_o,'guess':0},
    {'num':4,'symbol':"P",'mass': 30.97376151,'min':min_p,'max':max_p,'guess':0},
    {'num':5,'symbol':"S",'mass': 31.97207069,'min':min_s,'max':max_s,'guess':0},
    {'num':6,'symbol':"Na",'mass':22.989770,'min':min_na,'max':max_na,'guess':0},
    {'num':7,'symbol':"K",'mass': 38.963708,'min':min_k,'max':max_k,'guess':0},
    {'num':8,'symbol':"Cl",'mass':34.968852682,'min':min_cl,'max':max_cl,'guess':0}
    ]
    
    electron=0.0005485799094

    mass=mass+charge*electron #neutralize the molecule if charge is provided
    if tol_type=='ppm':
        tol=tol*mass/1e6
    if cho_only==True:
        hits = do_calculations_cho(mass,tol,elements,max_tests)
    else:
        hits = do_calculations(mass,tol,elements,max_tests)
    
    hits = sorted(hits,key=lambda x:x[0])

    formulae = [] #store all formulae
    for hit in hits:
        formula = [] #store list of elements and stoichiometry
        for element in hit[-1]:
            if element['guess'] != 0:
                formula.append('%s%d'%(element['symbol'],element['guess']))
        formulae.append((hit[:-1],''.join(formula)))
    
    return formulae

def calc_dbe(nC,nH):
    '''
    Calculate double bond equivalent given number of Carbon
    and number of Hydrogen atoms.
    
    From:
    
    W.C. Hockaday, A.M. Grannas, S. Kim, P.G. Hatcher, Org. Geochem. 37 (2006)
    501.
    
    Normalizeding dbe to nC is a measure of aromaticity
    
    condensed tannics are dbe/c around 0.6

    fused ring structures are dbe/c > 0.7
    '''
    dbe = 0.5 * ((nC*2.0)+2-nH)
    return dbe

def calc_mass(elements):
    """
    ?
    """
    sum = 0.0
    for el in elements:
        sum += el['mass'] * el['guess']
    return sum


def do_calculations_cho(mass,tol,elements,max_tests):
    """
    ?
    """
    limit_low = mass - tol
    limit_high = mass + tol
    test_counter = 0
    hits = []
    for n3 in range(elements[3]['min'],elements[3]['max']+1):
        elements[3]['guess'] = n3
        for n1 in range(elements[1]['min'],elements[1]['max']+1):
            elements[1]['guess'] = n1
            for n0 in range(elements[0]['min'],elements[0]['max']+1):
                elements[0]['guess'] = n0
                test_counter += 1
                if test_counter > max_tests:
                    print('ERROR test limit exceeded')
                    return
                theoretical_mass = calc_mass(elements)

                if ((theoretical_mass >= limit_low) & (theoretical_mass <= limit_high)):
                    hits.append((theoretical_mass,mass,abs(theoretical_mass-mass),deepcopy(elements)))

                if theoretical_mass > limit_high: # n0
                    break
        #         # if (theoretical_mass > limit_high) & (n0==elements[0]['min']):
        #             # break
        #     if (theoretical_mass > limit_high) & (n1==elements[1]['min']):
        #         break
        # if (theoretical_mass > limit_high) & (n3==elements[3]['min']):
        #     break
    return hits


def do_calculations(mass,tol,elements,max_tests):
    """
    ?
    """
    limit_low = mass - tol
    limit_high = mass + tol
    test_counter = 0
    hits = []
    for n8 in range(elements[8]['min'],elements[8]['max']+1):
        elements[8]['guess'] = n8
        for n7 in range(elements[7]['min'],elements[7]['max']+1):
            elements[7]['guess'] = n7
            for n6 in range(elements[6]['min'],elements[6]['max']+1):
                elements[6]['guess'] = n6
                for n5 in range(elements[5]['min'],elements[5]['max']+1):
                    elements[5]['guess'] = n5
                    for n4 in range(elements[4]['min'],elements[4]['max']+1):
                        elements[4]['guess'] = n4
                        for n3 in range(elements[3]['min'],elements[3]['max']+1):
                            elements[3]['guess'] = n3
                            for n2 in range(elements[2]['min'],elements[2]['max']+1):
                                elements[2]['guess'] = n2
                                for n1 in range(elements[1]['min'],elements[1]['max']+1):
                                    elements[1]['guess'] = n1
                                    for n0 in range(elements[0]['min'],elements[0]['max']+1):
                                        elements[0]['guess'] = n0
                                        test_counter += 1
                                        if test_counter > max_tests:
                                            print('ERROR test limit exceeded')
                                            return

                                        theoretical_mass = calc_mass(elements)

                                        if ((theoretical_mass >= limit_low) & (theoretical_mass <= limit_high)):
                                            hits.append((theoretical_mass,mass,abs(theoretical_mass-mass),deepcopy(elements)))

                                        if theoretical_mass > limit_high: # n0
                                            break
                                    if (theoretical_mass > limit_high) & (n0==elements[0]['min']):
                                        break
                                if (theoretical_mass > limit_high) & (n1==elements[1]['min']):
                                    break
                            if (theoretical_mass > limit_high) & (n2==elements[2]['min']):
                                break
                        if (theoretical_mass > limit_high) & (n3==elements[3]['min']):
                            break
                    if (theoretical_mass > limit_high) & (n4==elements[4]['min']):
                        break
                if (theoretical_mass > limit_high) & (n5==elements[5]['min']):
                    break
            if (theoretical_mass > limit_high) & (n6==elements[6]['min']):
                break
        if (theoretical_mass > limit_high) & (n7==elements[7]['min']):
            break
    return hits