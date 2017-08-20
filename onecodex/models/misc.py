import warnings

from onecodex.models import OneCodexBase
from onecodex.models.helpers import truncate_string


class Tags(OneCodexBase):
    _resource_path = '/api/v1/tags'

    def __repr__(self):
        return '<{} {}: "{}">'.format(self.__class__.__name__, self.id,
                                      truncate_string(self.name, 24))


class Users(OneCodexBase):
    _resource_path = '/api/v1/users'


class Projects(OneCodexBase):
    _resource_path = '/api/v1/projects'

    @classmethod
    def search_public(cls, *filters, **keyword_filters):
        warnings.warn('Now supported via `.where(..., public=True)`', DeprecationWarning)
        keyword_filters['public'] = True
        keyword_filters['limit'] = 1000
        return cls.where(*filters, **keyword_filters)


class Jobs(OneCodexBase):
    _resource_path = '/api/v1/jobs'
