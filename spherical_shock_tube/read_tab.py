from re import I
import numpy as np
import matplotlib.pyplot as plt
import sys


gamma = 5/3

def read_tab(fn):

    time = 0.0

    hdr = []
    try:
        with open(fn, 'r') as fp:
            while True:
                l = fp.readline()

                if l[0] != '#':
                    print(f'l[0]: {l}')
                    raise Exception("No header in %s found!" %(fn))

                if '# Athena++ data' in l:
                    for st in l[1:].split(" "):
                        if 'time' in st:
                            time = st.split("=")[1]
                if '# i' in l:
                    hdr = [i for i in l[1:].split(" ") if i!='' and i!='\n'] #[i.split(" ")[1].strip() for i in l[1:].split("[") if ']' in i]
                    break

        prim = np.loadtxt(fn, dtype={'names' : hdr, 'formats' : len(hdr) * (float,)})

    except:
        print(f"{fn} doesn't exist!")
        prim = -1

    return prim, float(time)

if __name__=="__main__":
    N = sys.argv[1]

    if int(N)>=0:
        fn = f'Sod.block0.out1.{N.zfill(5)}.tab'
        prim, time = read_tab(fn)

        fn = f'Sod.block0.out2.{N.zfill(5)}.tab'
        cons, time = read_tab(fn)


        plt.figure()

        plt.plot(prim['x1v'], prim['rho'], label = 'Density')
        plt.plot(prim['x1v'], prim['press'], label='Pressure')
        plt.plot(prim['x1v'], prim['press']/prim['rho'], label='Temperature')

        plt.yscale('log')
        plt.legend()
        plt.show()

    else:

        for n in range(np.abs(int(N))+1):

            print(f"Reading File {n}!")

            fn = f'Sod.block0.out1.{str(n).zfill(5)}.tab'
            prim = read_tab(fn)

            fn = f'Sod.block0.out2.{str(n).zfill(5)}.tab'
            cons = read_tab(fn)

            if prim==-1 or cons==-1:
                print("Last file reached!")
                break

            cs = np.sqrt(gamma * prim['press']/prim['rho'])

            plt.figure()

            # plt.plot(prim['x1v'], prim['rho'], label = 'Density')
            # plt.plot(prim['x1v'], prim['press'], label='Pressure')
            # plt.plot(prim['x1v'], prim['press']/prim['rho'], label='Temperature')
            plt.plot(prim['x1v'], prim['vel1'], label='vel_x')

            # plt.ylim(-5, 5)
            plt.ylim(-1, 1)

            # plt.axhline(-1, linestyle='dashed')
            # plt.axhline( 1, linestyle='dashed')

            # plt.yscale('log')
            plt.legend()
            # plt.savefig(f'Plots/spherical_{str(n).zfill(5)}.png')
            # plt.savefig(f'Plots/Mach/spherical_mach_{str(n).zfill(5)}.png')
            plt.savefig(f'Plots/vel/cartesian_vel_{str(n).zfill(5)}.png')

            print(f"Plotted File {n}!")