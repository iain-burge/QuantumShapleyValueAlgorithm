
import numpy as np
import matplotlib.pyplot as plt

from typing import Callable

def binomialFunct(n: int, m: int) -> Callable[[float], float]:
    return lambda p: (p**m)*((1-p)**(n-m))

def plotFunction(funct: Callable[[float], float], ax: ...,
        fRange: tuple[float, float] = (0,1), resolution: int = 100, 
    ) -> None:

    x = np.linspace(fRange[0], fRange[1], resolution)
    y = funct(x)

    ax.plot(x,y)

def getSamples(funct: Callable[[float], float], shift: bool = True,
        fRange: tuple[float, float] = (0,1), resolution: int = 100
    ) -> Callable[[float], float]:

    x = np.linspace(
        fRange[0] + 1/(2*resolution), fRange[1] - 1/(2*resolution), resolution) 
    return funct(x)

def plotApproxArea(hFunct: Callable[[float], float], ax: int, 
        samples: list[float] = None, widths: list[float] = None, 
        fRange: tuple[float, float] = (0,1), resolution: int = 100
    ) -> None:

    #Samples
    if samples is None:
        samples = np.linspace(fRange[0], fRange[1], resolution)
    
    #Bar widths
    if widths is None:
        widths    = abs(np.roll(samples, -1) - np.roll(samples, 1))/2
        widths[0] = widths[-1] = (1-np.sum(widths[1:-1]))/2

    #Bar placement
    barLefts = np.zeros(samples.shape)
    current  = 0
    for i, width in enumerate(widths[:-1]):
        barLefts[i] = current
        current    += width
    barLefts[-1] = current

    #Plotting bars / sampled points
    ax.plot(samples, hFunct(samples), 'bo')
    ax.bar(barLefts, hFunct(samples), width=widths,
        alpha=0.2, align='edge', edgecolor='b')
    ax.set(xlim=(0,1))

def initAxBinomialArea(
        funct: Callable[[float], float], ax: ..., areaRes: int = 100,
        functionRes: int = 1000
    ) -> None:

    widths = getSamples(lambda x: np.sin(np.pi*x), resolution=areaRes)
    widths /= np.sum(widths)

    samples = getSamples(lambda x: np.sin(np.pi*x/2)**2, resolution=areaRes)

    estArea  = np.sum([funct(sample)*width for sample, width in zip(samples,widths)])
    realArea = np.sum([getSamples(funct,resolution=functionRes)])/functionRes
    print(f"{estArea, realArea = }")
    print(f"delta = {abs(1-(estArea / realArea))}")

    plotFunction(funct, ax, resolution=functionRes)
    plotApproxArea(
        funct, 
        ax, 
        samples=samples,
        widths=widths
    )

def main():
    # plt.rcParams['text.usetex'] = True

    initEll = 0

    numRows = 2
    numCols = 3

    n = 4
    m = 1

    funct = binomialFunct(n,m)
    fig, ax = plt.subplots(numRows,numCols)

    # plt.suptitle(f"Binomial Area Approximation, {(n, m) = }")
    for i in range(numRows*numCols):
        curRow = i//numCols
        curCol = i% numCols

        initAxBinomialArea(
            funct, ax[curRow,curCol], areaRes=(2**(i+initEll))
        )

        ax[curRow, curCol].set_title(fr"$\ell$ = {i+initEll}")

        if curRow == numRows-1:
            ax[curRow, curCol].set_xlabel(fr"$x$")
        else:
            ax[curRow, curCol].get_xaxis().set_visible(False)

        if curCol == 0:
            ax[curRow, curCol].set_ylabel(r"$x^n(1-x)^{n-m}$")
        else:
            ax[curRow, curCol].get_yaxis().set_visible(False)


    plt.show()
    # plt.savefig(
    #     './plot.pdf', 
    #     bbox_inches="tight", 
    #     orientation="landscape",
    # )

if __name__ == "__main__":
    main()
