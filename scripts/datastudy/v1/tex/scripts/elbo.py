import matplotlib.pyplot as plt
import numpy as np
from sklearn.linear_model import LinearRegression
import sys
import re

def compute_residuals(k, Ks, dics):
    dics = (dics - np.mean(dics)) / np.std(dics)
    K1, K2 = Ks[:k, None], Ks[k-1:, None]
    dic1, dic2 = dics[:k, None], dics[k-1:, None]
    lm1 = LinearRegression().fit(K1, dic1)
    lm2 = LinearRegression().fit(K2, dic2)
    r1 = sum((lm1.predict(K1) - dic1) ** 2)[0]
    r2 = sum((lm2.predict(K2) - dic2) ** 2)[0]
    return r1 + r2

def find_elbo(Ks, dics):
    assert len(Ks) == len(dics)
    knots = np.arange(1, len(Ks))
    resids = [compute_residuals(k, Ks, dics) for k in knots]
    best_K = Ks[knots[np.argmin(resids)] - 1]
    return best_K

def genpath(simname, K, prefix=''):
    return f'{prefix}results/K={K}_{simname}/info.txt'

def read_dic(simname, K, **kwargs):
    path = genpath(simname, K, **kwargs)
    with open(path, 'r') as f:
        contents = f.read()
        dic = re.findall(r'(?<=DIC:\s)\d+\.\d+', contents)[0]
    return float(dic)


def get_best_K(simname, Ks=np.arange(2, 10), prefix=''):
    dics = np.array([read_dic(simname, K, prefix=prefix) for K in Ks])
    best_K = find_elbo(Ks, dics)
    return best_K


if __name__ == '__main__':
    if len(sys.argv[0]) > 0:
        prefix = ''
    else:
        prefix = '../'

    models = ['true', 'false']
    markers = ['CD3z', 'CD56', 'CD57', 'EOMES', 'Granzyme_A', 'LAG3',
               'Perforin', 'Siglec7']
    simnames = [f'marker={marker}_skewtmix={model}' for marker in markers for model in models]

    with open(f'{prefix}results/img/best-K.txt', 'w') as F:
        for simname in simnames:
            best_K = get_best_K(simname, prefix=prefix)
            marker = re.findall(r'(?<=marker=)\w+(?=_skewtmix)', simname)[0]
            model = re.findall(r'(?<=skewtmix=)\w+', simname)[0]
            F.write(f'marker: {marker} | skewtmix: {model} | K: {best_K}\n')

