class DihdForceConst(object):
    def __init__(self, value, dihd):
        self.forceconst = value
        self.dihd = dihd
        self.repr = dihd.repr

    def __repr__(self):
        return repr(self.forceconst)

    __str__ = __repr__

    def _call__(self, value):
        self.value = value
