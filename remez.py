import numpy as np
import matplotlib.pyplot as plt

class RemezContext:

    def __init__(self,num_taps):
        if num_taps % 2 == 0:
            raise ValueError("This implementation requires odd num_taps")

        self.numTaps = num_taps
        self.r = (num_taps - 1) // 2

        self.grid = []
        self.D = []
        self.W = []

        self.gridSize = 0

        self.ext = []

        self.x = []
        self.y = []
        self.ad = []
        self.E = []

        self.coef = []
        self.delta = 0.0


def create_dense_grid(bands,desired,weight,gridDensity,ctx):

    ctx.grid.clear()
    ctx.D.clear()
    ctx.W.clear()

    delta = 0.5 / (gridDensity * ctx.r)

    for band in range(len(weight)):

        low = bands[2 * band]
        high = bands[2 * band + 1]

        d = desired[2 * band]

        f = low

        while f < high:

            ctx.grid.append(f)
            ctx.D.append(d)
            ctx.W.append(weight[band])

            f += delta

        ctx.grid.append(high)
        ctx.D.append(d)
        ctx.W.append(weight[band])

    ctx.gridSize = len(ctx.grid)


def initial_guess(ctx):
    n = ctx.r + 2
    ctx.ext = []
    for i in range(n):
        index = int(
            round(i * (ctx.gridSize - 1) / (n - 1)))
        ctx.ext.append(index)


def calc_params(ctx):

    n = ctx.r + 2
    M = ctx.r

    A = np.zeros((n,M + 2))
    b = np.zeros(n)

    for i in range(n):

        index = ctx.ext[i]
        f = ctx.grid[index]

        for k in range(M + 1):
            A[i,k] = np.cos(
                np.pi * f * k
            )

        A[i,M + 1] = (-1.0) ** i / ctx.W[index]

        b[i] = ctx.D[index]

    solution = np.linalg.solve(A,b)

    ctx.coef = solution[:M + 1]
    ctx.delta = solution[M + 1]


def computeA(ctx,freq):

    value = 0.0

    for k in range(ctx.r + 1):
        value += (ctx.coef[k] *np.cos(np.pi * freq * k))
    return value


def calcError(ctx):

    ctx.E = np.zeros(ctx.gridSize)
    for i in range(ctx.gridSize):
        A = computeA(ctx,ctx.grid[i])
        ctx.E[i] = ctx.W[i] * (ctx.D[i] - A)


def find_candidates(ctx):

    candidates = []

    for i in range(1,ctx.gridSize - 1):

        if ctx.D[i] != ctx.D[i - 1]:
            candidates.append(i)
            continue

        if ctx.D[i] != ctx.D[i + 1]:
            candidates.append(i)
            continue

        left = abs(ctx.E[i - 1])
        current = abs(ctx.E[i])
        right = abs(ctx.E[i + 1])

        if current >= left and current >= right:
            candidates.append(i)

    candidates.append(0)
    candidates.append(ctx.gridSize - 1)

    return sorted(set(candidates))


def select_extrema(ctx,candidates):

    required = ctx.r + 2

    candidates = sorted(candidates,key=lambda i: ctx.grid[i])

    extrema = []

    for index in candidates:

        if len(extrema) == 0:
            extrema.append(index)
            continue

        previous = extrema[-1]
        sign_current = np.sign(ctx.E[index])
        sign_previous = np.sign(ctx.E[previous])
        if sign_current == 0:
            continue
        if sign_current != sign_previous:
            extrema.append(index)
        elif abs(ctx.E[index]) > abs(ctx.E[previous]):
            extrema[-1] = index
    if len(extrema) > required:
        while len(extrema) > required:
            smallest = min(range(len(extrema)),key=lambda i:abs(ctx.E[extrema[i]]) )
            extrema.pop(smallest)

    if len(extrema) < required:

        candidates = sorted(candidates,key=lambda i: abs(ctx.E[i]),reverse=True)
        candidates = sorted(candidates,key=lambda i: ctx.grid[i])

        extrema = []

        for index in candidates:

            if len(extrema) == 0:
                extrema.append(index)
                continue

            previous = extrema[-1]

            if np.sign(ctx.E[index]) == np.sign(ctx.E[previous]):
                continue

            extrema.append(index)

            if len(extrema) == required:
                break

    if len(extrema) != required:

        raise RuntimeError(f"Not enough alternating extrema: "f"{len(extrema)} found, "f"{required} required")
    return extrema


def search(ctx):

    candidates = find_candidates(ctx)
    ctx.ext = select_extrema(ctx,candidates)


def is_done(ctx):

    required = ctx.r + 2

    if len(ctx.ext) != required:
        return False

    errors = np.array([abs(ctx.E[i])for i in ctx.ext])

    emax = np.max(errors)
    emin = np.min(errors)

    if emax == 0:
        return True

    ripple_ok = ((emax - emin) / emax< 1e-6)
    extrema_ok = True

    for i in ctx.ext:

        if i > 0 and i < ctx.gridSize - 1:

            if abs(ctx.E[i]) < abs(ctx.E[i - 1]):
                extrema_ok = False

            if abs(ctx.E[i]) < abs(ctx.E[i + 1]):
                extrema_ok = False

    return ripple_ok and extrema_ok


def coef_to_fir(a,num_taps):

    M = (num_taps - 1) // 2

    h = np.zeros(num_taps)

    h[M] = a[0]

    for k in range(1,M + 1):

        value = a[k] / 2.0

        h[M - k] = value
        h[M + k] = value

    return h


def freq_sample(ctx):

    a = np.array(ctx.coef)

    h = coef_to_fir(a,ctx.numTaps)

    return h,ctx


def plot_remez_iteration(ctx,iteration):

    
    plt.plot( ctx.grid,ctx.E,label="Weighted error")
    plt.scatter([ctx.grid[i] for i in ctx.ext],[ctx.E[i] for i in ctx.ext],label="Extrema")
    plt.xlabel("Normalized frequency")
    plt.ylabel("Weighted error E(f)")
    plt.title(f"Remez error — iteration {iteration}")
    plt.legend()
    plt.grid()


def remez(num_taps,bands,desired,weight,grid_density=16):

    ctx = RemezContext(num_taps)

    create_dense_grid(
        bands,
        desired,
        weight,
        grid_density,
        ctx
    )

    initial_guess(ctx)
    plt.figure()
    for iteration in range(40):

        calc_params(ctx)
        calcError(ctx)
        plot_remez_iteration(ctx, iteration)
        if is_done(ctx):
            break

        search(ctx)

    calc_params(ctx)
    calcError(ctx)
    

    h,ctx = freq_sample(ctx)
    print("\nFinal extrema:")

    return h,ctx