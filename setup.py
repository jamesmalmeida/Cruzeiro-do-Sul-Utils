from distutils.core import setup
setup(
  name = 'czds_utils',         # How you named your package folder (MyLib)
  packages = ['czds_utils'],   # Chose the same as "name"
  version = '1.5',      # Start with a small number and increase it with every change you make
  license='GPL-3.0',        # Chose a license from here: https://help.github.com/articles/licensing-a-repository
  description = 'Libraries for the Cruzeiro do Sul Database',   # Give a short description about your library
  author = 'James Moraes de Almeida',                   # Type in your name
  author_email = 'james@almeida.page',      # Type in your E-Mail
  url = 'https://github.com/jamesmalmeida/Cruzeiro-do-Sul-Utils',   # Provide either the link to your github or to your website
  download_url = 'https://github.com/jamesmalmeida/Cruzeiro-do-Sul-Utils/archive/refs/tags/v1.5.tar.gz',
  keywords = ['Cruzeiro do Sul', 'Database'],   # Keywords that define your package best
  install_requires=[            # I get to this in a second
          'regex',
          'pandas',
          'numpy',
          'scipy',
      ],
  classifiers=[
    'Development Status :: 3 - Alpha',      # Chose either "3 - Alpha", "4 - Beta" or "5 - Production/Stable" as the current state of your package

    'Intended Audience :: Developers',      # Define that your audience are developers
    'Topic :: Software Development :: Build Tools',

    "License :: OSI Approved :: GNU General Public License (GPL)",

    'Programming Language :: Python :: 3',      #Specify which pyhton versions that you want to support
    'Programming Language :: Python :: 3.10',
  ],
)
