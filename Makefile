test: pytest

nosetest:
	nosetests --verbosity=2 --with-doctest obarasaika tests

pytest:
	pytest -v --doctest-modules obarasaika tests
