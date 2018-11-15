import re

__author__ = 'lyschoening'


def camel_case(s):
    return s[0].lower() + s.title().replace('-', '').replace('_', '')[1:] if s else s


def upper_camel_case(s):
    return s.title().replace('-', '').replace('_', '') if s else s


def snake_case(s):
    return s[0].lower() + re.sub('([A-Z])', r'_\1', s[1:]).lower()


def escape(html):
    """Returns the given HTML with ampersands, quotes and carets encoded."""
    return html \
        .replace('&', '&amp;') \
        .replace('<', '&lt;') \
        .replace('>', '&gt;') \
        .replace('"', '&quot;') \
        .replace("'", '&#39;')
