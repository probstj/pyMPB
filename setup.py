from setuptools import setup

setup(name='pympb',
      version='1.0',
      description='Python package for running simulations with MIT Photonic-Bands (MPB)',
      url='http://github.com/probstj/pympb',
      author='Jürgen Probst',
      author_email='juergen.probst@gmail.com',
      license='GPLv3',
      packages=['pympb'],
      install_requires=[
          'numpy', 
          'matplotlib',
      ],
      zip_safe=False)
