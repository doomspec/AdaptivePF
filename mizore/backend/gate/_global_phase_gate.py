from qulacs.gate import DenseMatrix
from numpy import exp


def GlobalPhaseGate(phase):
    return DenseMatrix(0, [[exp(1j * phase), 0], [0, exp(1j * phase)]])


if __name__ == '__main__':
    from math import pi

    print(GlobalPhaseGate(2 * pi))
