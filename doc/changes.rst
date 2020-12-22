=====================
prospect's Change Log
=====================

1.0.1 (unreleased)
------------------

* No changes yet.

1.0.0 (2020-12-22)
------------------

*This is a major refactor to allow prospect to be ``pip``-installable,
along with other standard Python package features.*  Some API changes should
be expected.  See PR `#54`_ for details.

.. _`#54`: https://github.com/desihub/prospect/pull/54

0.3.0 (2020-10-27)
------------------

* Planned final reference tag before package refactor.
* Allow prospect to be loaded in an environment without the DESI software stack (PR `#50`_).
* Allow specutils objects to be plotted (PR `#43`_).

.. _`#50`: https://github.com/desihub/prospect/pull/50
.. _`#43`: https://github.com/desihub/prospect/pull/43

0.2.2 (2020-07-04)
------------------

* Command line scripts (particularly for static HTML display) and package documentation (PR `#48`_).

.. _`#48`: https://github.com/desihub/prospect/pull/48

0.2.1 (2020-06-12)
------------------

* Fix some data handling issues related to Andes release (`dbcde2f`_).

.. _`dbcde2f`: https://github.com/desihub/prospect/commit/dbcde2f0be2b13e96138a9fbac036f083e2f7b24)

0.2.0 (2020-06-12)
------------------

Summary of new features from this branch (PR `#42`_):

- display "second model" which can be the Nth best fit, or standard templates
- table with redrock results shows N best fits
- ivar-weighting when smoothing
- imaging cross-hair
- line list
- several bug fixes + code restructure
- new widgets (*e.g.* "standard VI comments") + widget layout

.. _`#42`: https://github.com/desihub/prospect/pull/42

0.1.1 (2020-04-07)
------------------

* Static HTML pages (PR `#39`_).

.. _`#39`: https://github.com/desihub/prospect/pull/39

0.1.0 (2020-04-01)
------------------

* Initial reference tag.
