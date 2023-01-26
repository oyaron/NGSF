import statistics
import numpy as np
import scipy.signal as mf


def linear_error(spec_object):

    flux = spec_object[:, 1]
    lam = spec_object[:, 0]

    # For how many points do we make the lines
    num = 10

    if len(flux) % num != 0:
        c = len(flux) % num
        flux = flux[:-c]
        lams = lam[:-c]

    else:
        lams = lam
        c = 0

    flux_new = flux.reshape((-1, num))
    lam_new = lams.reshape((-1, num))
    m = []
    b = []
    sigma = []
    r = []

    for n in range(len(lam_new)):

        a = np.polyfit(lam_new[n], flux_new[n], 1)
        m.append(a[0])
        b.append(a[1])
        y = m[n] * lam_new[n] + b[n]

        r = flux_new - y
    for i in r:
        s = statistics.stdev(i)
        sigma.append(s)

    # Here we make the error be the same size as the original lambda and then
    # take the transpose

    error = list(np.repeat(sigma, num))
    ls = [error[-1]] * c
    error = error + ls

    error = np.asarray(error)

    return np.array([lam, error]).T


def savitzky_golay(spec):

    x = spec[:, 0]
    y = spec[:, 1] / spec[:, 1].mean()

    # Find residuals from smooth line
    smooth = mf.savgol_filter(
        y, 31, 3, deriv=0, delta=1.0, axis=-1, mode="nearest", cval=0.0
    )
    resid = y - smooth

    # Calculate the variance
    def moving_average(a, n=3):
        ret = np.cumsum(a, dtype=float)
        ret[n:] = ret[n:] - ret[:-n]
        return ret[n - 1:] / n

    mov_var = moving_average(resid**2, n=100)
    mov_var = np.concatenate((mov_var, [mov_var[-1]] * (resid.size -
                                                        mov_var.size)))
    err_std = np.sqrt(mov_var)

    for i in range(0, len(err_std)):
        if err_std[i] == 0:
            err_std[i] = 1e-50

    return np.array([x, err_std]).T
