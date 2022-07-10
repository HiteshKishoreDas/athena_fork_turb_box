from re import I
import numpy as np
import matplotlib.pyplot as plt
import sys


def read_tab(fn):

    hdr = []
    try:
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

    except:
        print(f"{fn} doesn't exist!")
        prim = -1

    return prim

N = sys.argv[1]

if int(N)>=0:
    fn = f'Sod.block0.out1.{N.zfill(5)}.tab'
    prim = read_tab(fn)

    fn = f'Sod.block0.out2.{N.zfill(5)}.tab'
    cons = read_tab(fn)


    plt.figure()

    plt.plot(prim['x1v'], prim['rho'], label = 'Density')
    plt.plot(prim['x1v'], prim['press'], label='Pressure')
    # plt.plot(prim['x1v'], prim['press']/prim['rho'], label='Temperature')

    plt.yscale('log')
    plt.legend()
    plt.show()

else:

    for n in range(np.abs(int(N))):

        print(f"Reading File {n}!")

        fn = f'Sod.block0.out1.{str(n).zfill(5)}.tab'
        prim = read_tab(fn)

        fn = f'Sod.block0.out2.{str(n).zfill(5)}.tab'
        cons = read_tab(fn)

        if prim==-1 or cons==-1:
            print("Last file reached!")
            break

        
        plt.figure()

        plt.plot(prim['x1v'], prim['rho'], label = 'Density')
        plt.plot(prim['x1v'], prim['press'], label='Pressure')
        # plt.plot(prim['x1v'], prim['press']/prim['rho'], label='Temperature')

        plt.ylim(5e-1, 2e3)
    
        plt.yscale('log')
        plt.legend()
        plt.savefig(f'Plots/spherical_{str(n).zfill(5)}.png')
        # plt.savefig(f'Plots/cartesian_{str(n).zfill(5)}.png')
        
        print(f"Plotted File {n}!")