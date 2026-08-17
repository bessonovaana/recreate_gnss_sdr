import numpy as np
import matplotlib.pyplot as plt

class RemezContext:

    def __init__(self, num_taps):

        self.numTaps = num_taps

        # число неизвестных косинусных коэффициентов
        self.r = num_taps // 2

        self.grid = []
        self.D = []
        self.W = []

        self.gridSize = 0

        self.ext = []

        self.x = []
        self.y = []
        self.ad = []
        self.E = []


def create_dense_grid(bands,
                      desired,
                      weight,
                      gridDensity,
                      ctx):

    ctx.grid.clear()
    ctx.D.clear()
    ctx.W.clear()

    delta = 0.5 / (gridDensity * ctx.r)

    for band in range(len(weight)):

        low = bands[2 * band]
        high = bands[2 * band + 1]

        f = low

        while f <= high + 1e-12:

            ctx.grid.append(f)

            ctx.D.append(desired[2 * band])

            ctx.W.append(weight[band])

            f += delta

    ctx.gridSize = len(ctx.grid)


def initial_guess(ctx):

    n = ctx.r + 2

    ctx.ext = [0] * n

    for i in range(n):

        ctx.ext[i] = (i * (ctx.gridSize - 1)) // (ctx.r + 1)
    plt.plot(ctx.grid, np.zeros(ctx.gridSize))
    plt.scatter([ctx.grid[i] for i in ctx.ext],np.zeros(len(ctx.ext)),label='Initial extrema')

    plt.xlabel("Normalized frequency")
    plt.yticks([])

    plt.title("Initial guess of extremal frequencies")

    plt.grid(True)

def calc_params(ctx):

    n = ctx.r + 2

    ctx.x = [0.0] * n
    ctx.y = [0.0] * n
    ctx.ad = [0.0] * n

    for i in range(n):

        ctx.x[i] = np.cos(2 * np.pi * ctx.grid[ctx.ext[i]])

    for i in range(n):

        prod = 1.0

        for k in range(n):

            if i == k:
                continue

            prod *= (ctx.x[i] - ctx.x[k])

        ctx.ad[i] = 1.0 / prod

    numerator = 0.0
    denominator = 0.0

    sign = 1

    for i in range(n):

        numerator += ctx.ad[i] * ctx.D[ctx.ext[i]]

        denominator += sign * ctx.ad[i] / ctx.W[ctx.ext[i]]

        sign = -sign

    delta = numerator / denominator

    sign = 1

    for i in range(n):

        ctx.y[i] = ctx.D[ctx.ext[i]] - sign * delta / ctx.W[ctx.ext[i]]

        sign = -sign


def computeA(ctx, freq):

    xc = np.cos(2 * np.pi * freq)

    numer = 0.0
    denom = 0.0

    eps = 1e-12

    for i in range(ctx.r + 2):

        diff = xc - ctx.x[i]

        if abs(diff) < eps:

            return ctx.y[i]

        c = ctx.ad[i] / diff

        numer += c * ctx.y[i]

        denom += c

    return numer / denom


def calcError(ctx):

    ctx.E = [0.0] * ctx.gridSize

    for i in range(ctx.gridSize):

        A = computeA(ctx, ctx.grid[i])

        ctx.E[i] = ctx.W[i] * (ctx.D[i] - A)


def search(ctx):

    extrema = [0]

    for i in range(1, ctx.gridSize - 1):

        if abs(ctx.E[i]) >= abs(ctx.E[i - 1]) and \
           abs(ctx.E[i]) >= abs(ctx.E[i + 1]):

            extrema.append(i)

    extrema.append(ctx.gridSize - 1)

    alternating = [extrema[0]]

    sign = np.sign(ctx.E[extrema[0]])

    for idx in extrema[1:]:

        s = np.sign(ctx.E[idx])

        if s != sign:

            alternating.append(idx)

            sign = s

        else:

            if abs(ctx.E[idx]) > abs(ctx.E[alternating[-1]]):

                alternating[-1] = idx

    while len(alternating) > ctx.r + 2:

        smallest = min(
            range(len(alternating)),
            key=lambda i: abs(ctx.E[alternating[i]])
        )

        alternating.pop(smallest)

    if len(alternating) < ctx.r + 2:

        raise RuntimeError("Not enough extrema")

    ctx.ext = alternating
    
    plt.plot(ctx.grid, ctx.E, label = "Ошибка")
    plt.scatter([ctx.grid[i] for i in ctx.ext], [ctx.E[i] for i in ctx.ext])
    plt.xlabel("Normalized frequency")
    plt.ylabel("Weighted error E(f)")
    plt.legend()

def is_done(ctx):

    err = [abs(ctx.E[i]) for i in ctx.ext]

    emax = max(err)
    emin = min(err)

    return (emax - emin) / emax < 1e-4



def get_coef(ctx):
    M=ctx.r
    K=M+2

    freqs = np.array([ctx.grid[i] for i in ctx.ext])
    y=np.array(ctx.y)
    X = np.zeros((K, M + 1))

    for i in range(K):

        for k in range(M + 1):

            X[i, k] = np.cos(2 * np.pi * freqs[i] * k)
    a, *_ = np.linalg.lstsq(X, y, rcond=None)

    return a

def coef_to_fir(a, num_taps):
    M = (num_taps - 1) // 2

    h = np.zeros(num_taps)
    for k in range(1, M + 1):

        value = a[k] / 2.0

        h[M - k] = value
        h[M + k] = value

    return h

def remez(num_taps,
          bands,
          desired,
          weight,
          grid_density=16):

    ctx = RemezContext(num_taps)

    create_dense_grid(
        bands,
        desired,
        weight,
        grid_density,
        ctx
    )

    initial_guess(ctx)

    for _ in range(40):

        calc_params(ctx)

        calcError(ctx)

        search(ctx)

        if is_done(ctx):
            break

    calc_params(ctx)

    return freq_sample(ctx)





class ParksMaclellan:

    def __init__(
        self,
        num_taps,
        bands,
        desired,
        weight,
        grid_density=16
    ):

        self.num_taps = num_taps
        self.bands = bands
        self.desired = desired
        self.weight = weight
        self.grid_density = grid_density

    def design(self):

        h, ctx = remez(
            self.num_taps,
            self.bands,
            self.desired,
            self.weight,
            self.grid_density
        )

        return h, ctx
    



number_of_taps = 5

bands = [
    0.0,
    0.44,
    0.55,
    1.0
]

desired = [
    1.0,
    1.0,
    0.0,
    0.0
]

weight = [
    1.0,
    1.0
]

h = remez(
    number_of_taps,
    bands,
    desired,
    weight,
    grid_density=16
)

print(h)

plt.show()