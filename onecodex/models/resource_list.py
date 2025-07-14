from onecodex.exceptions import OneCodexException
from onecodex.models.base import OneCodexBase

DEFAULT_PAGE_SIZE = 200


# TODO: Check these changes *very* carefully
class ResourceList(object):
    """Wrapper around lists of OneCodexBase API objects.

    Parameters
    ----------
    _res_list : `list`
        A list of OneCodexBase model objects.
    oc_model : `OneCodexBase`
        A class which inherits from `OneCodexBase` (e.g., Samples, Classifications, Projects).

    Notes
    -----
    In OneCodexBase, when attributes are lists (e.g., `Samples.tags`), actions performed on the
    returned lists are not passed through to the underlying object's list. This class passes
    those actions through, and will generally act like a list.
    """

    def __init__(self, _res_list, oc_model, **kwargs):
        if not issubclass(oc_model, OneCodexBase):
            raise ValueError(
                "Expected object of type '{}', got '{}'".format(
                    OneCodexBase.__name__, oc_model.__name__
                )
            )

        # turn potion Resource objects into OneCodex objects
        self._oc_model = oc_model
        self._kwargs = kwargs
        self._res_list = _res_list

    def _check_valid_resource(self, other, check_for_dupes=True):
        if not isinstance(other, list):
            other = [other]

        other_ids = []
        for o in other:
            if not isinstance(o, self._oc_model):
                raise ValueError(
                    "Expected object of type '{}', got '{}'".format(
                        self._oc_model.__name__, type(o).__name__
                    )
                )

            other_ids.append(o.id)

        if check_for_dupes:
            # duplicates are not allowed
            self_ids = [s.id for s in self._res_list]

            if len(set(self_ids + other_ids)) != len(self_ids + other_ids):
                raise OneCodexException(
                    "{} cannot contain duplicate objects".format(self.__class__.__name__)
                )

    def __eq__(self, other):
        if not isinstance(other, self.__class__):
            return False

        return all(id(x) == id(y) for x, y in zip(self._res_list, other._res_list))

        # two ResourceLists are equal if they refer to the same underlying Resource
        return id(self._res_list) == id(other._res_list)

    def __contains__(self, other):
        return other.__hash__() in [x.__hash__() for x in self._res_list]

    @property
    def __repr__(self):
        return self._res_list.__repr__

    @property
    def __len__(self):
        return self._res_list.__len__

    def __getitem__(self, x):
        wrapped = self._res_list[x]
        if isinstance(wrapped, list):
            return self.__class__(self._res_list[x], self._oc_model, **self._kwargs)
        else:
            return wrapped

    def __setitem__(self, k, v):
        self._check_valid_resource(v)
        self._res_list[k] = v

    def __delitem__(self, x):
        del self._res_list[x]

    @property
    def __iter__(self):
        return self._res_list.__iter__

    @property
    def __reversed__(self):
        return self._res_list.__reversed__

    def __add__(self, other):
        if not isinstance(other, self.__class__):
            raise TypeError(
                'can only concatenate {} (not "{}") to {}'.format(
                    self.__class__.__name__, type(other), self.__class__.__name__
                )
            )
        new_obj = self.copy()
        new_obj.extend(other._res_list)
        return new_obj

    def append(self, x):
        self._check_valid_resource(x)
        self._res_list.append(x)

    def clear(self):
        self._res_list.clear()

    def copy(self):
        new_obj = self.__class__(self._res_list[:], self._oc_model, **self._kwargs)
        return new_obj

    def count(self, x):
        # assume that ResourceList objects are identical if they share the same underlying resource
        self._check_valid_resource(x, check_for_dupes=False)
        n = 0
        for res_obj in self._res_list:
            if res_obj == x:
                n += 1
        return n

    def extend(self, iterable):
        self._check_valid_resource(iterable)
        self._res_list.extend([x for x in iterable])

    def index(self, x):
        # assume that ResourceList objects are identical if they share the same underlying resource
        self._check_valid_resource(x, check_for_dupes=False)
        for res_obj_idx, res_obj in enumerate(self._res_list):
            if res_obj == x:
                return res_obj_idx
        raise ValueError("{} is not in list".format(x))

    def insert(self, idx, x):
        self._check_valid_resource(x)
        self._res_list.insert(idx, x)

    def pop(self):
        return self._res_list.pop()

    def remove(self, x):
        del self._res_list[self.index(x)]
