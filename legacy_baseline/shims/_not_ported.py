"""Helpers for import stubs of out-of-scope (image-analysis) modules.

Used only by the legacy baseline container: the stubs sit ahead of the original source on
PYTHONPATH so the *unmodified* original downstream modules can be imported without cv2,
scikit-image, image data or the missing SkeletonAnalysis/dataValidation modules.
"""
MESSAGE = "image analysis not yet ported"


def stub(name):
    def _f(*args, **kwargs):
        raise NotImplementedError("%s (%s)" % (MESSAGE, name))
    _f.__name__ = name.split('.')[-1]
    return _f


class StubClass(object):
    def __init__(self, *args, **kwargs):
        raise NotImplementedError("%s (%s)" % (MESSAGE, type(self).__name__))
