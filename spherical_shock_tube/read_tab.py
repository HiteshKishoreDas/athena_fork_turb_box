import numpy as np
import matplotlib.pyplot as plt
import sys


N = sys.argv[1]
fn = f'Sod.block0.out1.{N.zfill(5)}.tab'
# fn = 'Sod.hst'

hdr = []
with open(fn, 'r') as fp:
    while True:
        l = fp.readline()

        if l[0] != '#':
            print(f'l[0]: {l}')
            raise Exception("No header in %s found!" %(fn))
        if '# i' in l:
            hdr = [i for i in l[1:].split(" ") if i!='' and i!='\n'] #[i.split(" ")[1].strip() for i in l[1:].split("[") if ']' in i]
            break

prim = np.loadtxt(fn, dtype={'names' : hdr, 'formats' : len(hdr) * (float,)})


N = sys.argv[1]
fn = f'Sod.block0.out2.{N.zfill(5)}.tab'
# fn = 'Sod.hst'

hdr = []
with open(fn, 'r') as fp:
    while True:
        l = fp.readline()

        if l[0] != '#':
            print(f'l[0]: {l}')
            raise Exception("No header in %s found!" %(fn))
        if '# i' in l:
            hdr = [i for i in l[1:].split(" ") if i!='' and i!='\n'] #[i.split(" ")[1].strip() for i in l[1:].split("[") if ']' in i]
            break

cons = np.loadtxt(fn, dtype={'names' : hdr, 'formats' : len(hdr) * (float,)})

plt.figure()

plt.plot(prim['x1v'], prim['press'], label='Pressure')
plt.plot(prim['x1v'], prim['rho'], label = 'Density')
plt.plot(prim['x1v'], prim['press']/prim['rho'], label='Temperature')

plt.plot(cons['x1v'], cons['Etot'], label='Total Energy')

plt.legend()
plt.yscale('log')
plt.show()

print(np.average(prim['press']/cons['Etot']))   