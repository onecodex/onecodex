# flake8: noqa
import collections.abc
from pprint import pformat

from .utils import escape


class PaginatedList(collections.abc.Sequence):
    def __init__(self, binding, params):
        self._pages = {}
        self._per_page = per_page = params.pop('per_page', 20)
        self._binding = binding
        self._total_count = 0
        self._request_params = params
        self.fetch_page(1, per_page)

    def __getitem__(self, item):
        if isinstance(item, slice):
            return [self.__getitem__(index) for index in range(*item.indices(self._total_count))]

        if item < 0 or item >= self._total_count:
            raise IndexError()

        page, offset = item // self._per_page + 1, item % self._per_page
        if page not in self._pages:
            self.fetch_page(page, self._per_page)
        return self._pages[page][offset]

    def __len__(self):
        return self._total_count

    def fetch_page(self, page, per_page):
        params = dict(page=page, per_page=per_page)
        params.update(self._request_params)
        response, response_data = self._binding.make_request(None, params)

        try:
            self._total_count = int(response.headers['X-Total-Count'])
        except KeyError:
            self._total_count = len(response_data)

        self._pages[page] = response_data

    def _repr_html_(self):
        if len(self) <= 10:
            items = [escape(pformat(item)) for item in self[:]]
        else:
            items = [escape(pformat(item)) for item in self[:5]] + \
                    ['<strong>...</strong>'] + \
                    [escape(pformat(item)) for item in self[-5:]]

        return '''<table>
        <thead>
            <tr>
                <th><code>PaginatedList({params})</code>: <em>{count}</em> items</th>
            </tr>
        </thead>
        <tbody>{items}</tbody>
        </table>'''.format(
            params=', '.join(['{}.{}'.format(self._binding.owner.__name__, self._binding.link.rel)] +
                             ['{}={}'.format(k, repr(v)) for k, v in self._request_params.items()]),
            count=self._total_count,
            items='\n'.join('<tr><td><code>{}</code></td></tr>'.format(item) for item in items))

    def __repr__(self):
        return 'PaginatedList({params})'.format(params=', '.join(
            ['{}.{}'.format(self._binding.owner.__name__, self._binding.link.rel)] +
            ['{}={}'.format(k, repr(v)) for k, v in self._request_params.items()]), )
