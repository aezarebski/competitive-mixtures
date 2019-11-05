# Changelog for competitive-mixtures

## 2019-10-24

* Clean up the dumping code so that it does not rely on hard-coded values for the number of days the experiment ran for.

## 2019-10-23

* Use `packrat` to ensure the correct versions of packages are used.
* Fix a bug where the data needed by Stan was not being found because it was in the global environment.

## 2019-09-01

* Include some safety checks before reading and writing files to assist in debugging.
* Include a reference to the publication that used this code.
