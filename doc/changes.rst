=====================
prospect's Change Log
=====================

1.3.1 (unreleased)
------------------

* Don't use mask information when marking bad pixels in SDSS data (PR `#94`_).

.. _`#94`: https://github.com/desihub/prospect/pull/94

1.3.0 (2023-09-06)
------------------

* Added/renamed options to be run with :command:`prospect_pages` (PR `#90`_):
  - ``--no_imaging``
  - ``--no_noise``
  - ``--no_thumb_tab``
  - ``--no_vi_widgets``
  - ``--no_coaddcam``
* Renamed options in prospect_pages:
  - ``--nspecperfile`` is renamed ``--nspec_per_page``
  - ``--no-clean_fiberstatus`` is renamed ``--no_clean_fiberstatus``

.. _`#90`: https://github.com/desihub/prospect/pull/90

1.2.5 (2023-06-14)
------------------

* Update spectrum service notebook to `NOIRLab SPARCL`_; ensure continued support
  of both SDSS and BOSS spectra (PR `#87`_).
* Handle high-proper motion objects: a second cross-hair is shown (Issue `#84`_).
* Adapt code and examples to recent DESI releases, up to iron (Issue `#82`_).
* Minor bug fixes.

.. _`NOIRLab SPARCL`: https://astrosparcl.datalab.noirlab.edu/
.. _`#87`: https://github.com/desihub/prospect/pull/87
.. _`#84`: https://github.com/desihub/prospect/issues/84
.. _`#82`: https://github.com/desihub/prospect/issues/82

1.2.4 (2023-01-09)
------------------

* Workaround masked ``SUBTYPE`` when loading templates (PR `#81`_).
* Update API documentation and test infrastructure (PR `#79`_).

.. _`#81`: https://github.com/desihub/prospect/pull/81
.. _`#79`: https://github.com/desihub/prospect/pull/79

1.2.3 (2022-09-15)
------------------

* "Standard templates" now provided from standalone files, independently of Redrock (Issue `#68`_).
* Handling of different Redrock templates (Issue `#77`_).
* Solved oversmoothing issue (Issue `#60`_).
* New (convenience) options in ``prospect_pages.py``.

.. _`#68`: https://github.com/desihub/prospect/issues/68
.. _`#77`: https://github.com/desihub/prospect/issues/77
.. _`#60`: https://github.com/desihub/prospect/issues/60

1.2.2 (2022-05-24)
------------------

* Support for astropy 5 (redrock files I/O).

1.2.1 (2022-03-02)
------------------

* Bug fixes and support for astropy 5 (PR `#75`_).

.. _`#75`: https://github.com/desihub/prospect/pull/75

1.2.0 (2021-07-29)
------------------

* Tuned widget layout, move legend location (Issues `#61`_ and `#63`_), added "Redrock_deltachi2" in VI outputs.
* Various bug fixes, in particular handling z_input / rest-frame issues (Issue `#44`_)
* Code cleaning (PR `#73`_).
* Now a single script to create static html pages for DESI: ``prospect_pages.py``.
* Two notebooks updated, and a script ``example_prospect_pages``, to be used both as doc and tests.
* Compatible with files and directory trees from andes to everest.

.. _`#73`: https://github.com/desihub/prospect/pull/73
.. _`#61`: https://github.com/desihub/prospect/issues/61
.. _`#63`: https://github.com/desihub/prospect/issues/63
.. _`#44`: https://github.com/desihub/prospect/issues/44

1.1.1 (2021-03-31)
------------------

* Make specutils_ imports optional, for DESI-only users.
* Improved metadata internal handling and display in VI pages,
  including FIRST/LAST/NUM_EXPID/TILEID/NIGHT/FIBER, MORPHTYPE,
  redrock version, and support all current DESI targeting masks (Issues `#21`_, `#51`_ and `#55`_).
* List of "major" spectral lines updated (Issue `#69`_).

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
