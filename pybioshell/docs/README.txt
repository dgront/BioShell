To compile the documentation:
  1. start in bioshell4/pybioshell
  2. build the library with Maturin; make sure the default python can import the bioshell library
      maturin develop
  3. remove the old docs:
      rm -rf docs/_build/
  4. run Sphinx:
      sphinx-autobuild docs docs/_build/html --open-browser