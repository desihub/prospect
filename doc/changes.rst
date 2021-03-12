=====================
prospect's Change Log
=====================

1.1.1 (unreleased)
------------------

* Make specutils_ imports optional, for DESI-only users.
* Improved metadata internal handling and display in VI pages, including FIRST/LAST/NUM_EXPID/TILEID/NIGHT/FIBER, MORPHTYPE, redrock version, and support all current DESI targeting masks (PR `#21`_, `#51`_ and `#55`_).
* List of "major" spectral lines updated (PR `#69`_).

.. _specutils: https://specutils.readthedocs.io
.. _`#21`: https://github.com/desihub/prospect/issues/21
.. _`#51`: https://github.com/desihub/prospect/issues/51
.. _`#55`: https://github.com/desihub/prospect/issues/55
.. _`#69`: https://github.com/desihub/prospect/issues/69

1.1.0 (2021-02-10)
------------------

* Restructure :mod:`prospect.viewer` into classes (PR `#67`_).

.. _`#67`: https://github.com/desihub/prospect/pull/67

1.0.2 (2021-01-27)
------------------

* Optional countdown widget (PR `#58`_).
* Updated requests to legacysurvey (PR `#65`_)
* Merge of DESI and specutils code in :mod:`prospect.viewer` (PR `#66`_).

.. _`#58`: https://github.com/desihub/prospect/pull/58
.. _`#65`: https://github.com/desihub/prospect/pull/65
.. _`#66`: https://github.com/desihub/prospect/pull/66

1.0.1 (2021-01-07)
------------------

* Make desiutil_ imports optional in ``setup.py`` (PR `#56`_).
* Adapt scripts to blanc ("deep" directories) (8552fd8e_).
* Support ``SV1_BGS_TARGET`` (7a5ca41f_).

.. _desiutil: https://github.com/desihub/desiutil
.. _`#56`: https://github.com/desihub/prospect/pull/56
.. _8552fd8e: https://github.com/desihub/prospect/commit/8552fd8ec1801d322e9df3b468ed319109410763
.. _7a5ca41f: https://github.com/desihub/prospect/commit/7a5ca41f41d1e7475c579b256b1e9fdccafe530f

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
* Allow specutils_ objects to be plotted (PR `#43`_).

.. _`#50`: https://github.com/desihub/prospect/pull/50
.. _`#43`: https://github.com/desihub/prospect/pull/43
.. _specutils: https://specutils.readthedocs.io

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
