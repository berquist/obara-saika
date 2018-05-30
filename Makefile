test: pytest-cov

nosetest:
	nosetests --verbosity=2 --with-doctest obarasaika tests

pytest:
	pytest -v --doctest-modules obarasaika tests

pytest-cov:
	pytest -v --doctest-modules --cov=obarasaika obarasaika tests
