
=============
Potion client
=============


.. image:: https://img.shields.io/travis/biosustain/potion-client/new-potion-client.svg?style=flat-square
    :target: https://travis-ci.org/biosustain/potion-client

.. image:: https://img.shields.io/coveralls/biosustain/potion-client/new-potion-client.svg?style=flat-square
    :target: https://coveralls.io/r/biosustain/potion-client

.. image:: https://img.shields.io/pypi/v/Potion-Client.svg?style=flat-square
    :target: https://pypi.python.org/pypi/Potion-Client

.. image:: https://img.shields.io/pypi/l/Potion-Client.svg?style=flat-square
    :target: https://pypi.python.org/pypi/Potion-Client

.. image:: https://badges.gitter.im/Join%20Chat.svg
   :alt: Join the chat at https://gitter.im/biosustain/potion
   :target: https://gitter.im/biosustain/potion?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge

Description
===========

This is a Python client for APIs written in `Flask-Potion <https://github.com/biosustain/potion>`_ (a powerful Flask extension for self-documenting JSON APIs).

The package uses `Requests <https://github.com/kennethreitz/requests>`_ to provide a super-simple interface to Potion APIs that
works with all common authentication methods. It generates classes for each of the resources in the API and automatically handles pagination
and resolving and serializing references. It also has some basic `IPython Notebook <http://ipython.org/notebook.html>`_ support.

Example
=======

.. code-block:: python

    from potion_client import Client
    from potion_client.auth import HTTPBearerAuth
    from potion_client.exceptions import ItemNotFound

    client = Client('http://localhost/api', auth=HTTPBearerAuth('79054025255fb1a26e4bc422aef54eb4'))

    u123 = client.User(123)

    chomp = client.Animal()
    chomp.owner = u123
    chomp.name = "Chomp"
    chomp.species = "hamster"
    chomp.save()

    pets = client.Animal.instances(where={"owner": u123}, sort={"created_at": True})

    print("{} has {} pet(s)".format(u123.first_name, len(pets))

    for pet in pets:
        if pet is not chomp:
            pet.add_friend(chomp)
            print("{} is now friends with Chomp".format(pet.name)))

    try:
        foo = client.User.first(where={"username": "foo"})
    except ItemNotFound:
        print("User 'foo' does not exist!")
    else:
        chomp.update(owner=foo)
        print("Chomp has been sold to {}".format(foo.name))

    chomp.destroy()
    print("RIP, Chomp. You lived a happy life.")

Installation
============

To install ``potion-client``, run:

::

    pip install potion-client




Authors
=======

Potion-client was written by `João Cardoso <https://github.com/joaocardoso>`_ and `Lars Schöning <https://github.com/lyschoening>`_.