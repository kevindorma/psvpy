from distutils.core import setup
setup(
  name = 'psvpy',         # How you named your package folder (MyLib)
  package = 'psvpy',
  version = '0.5',      # Start with a small number and increase it with every change you make
  license='MIT',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'PSV calculations',   # Give a short description about your library
  author = 'Kevin Dorma',                   # Type in your name
  author_email = 'kevin@kevindorma.ca',      # Type in your E-Mail
  url = 'https://github.com/kevindorma/psvpy',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/kevindorma/psvpy/archive/v0.5.tar.gz',    # I explain this later on
  keywords = ['PSV', 'engineering', 'safety'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'numpy',
          'pandas',
          'scipy',
      ],
)