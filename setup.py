import setuptools

setuptools.setup(
    name='scampy',
    version='1.0.0',
    packages=['scampy'],
    license=['MIT'],
    install_requires=[
          'numpy',
          'sympy',
          'scipy'
    ]
)