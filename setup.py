try : 
    from setuptools import setup
except ImportError :
    from distutils.core import setup


config = {
    'description' : 'My Project',
    'author' : 'Subhodeep Moitra',
    'url' : 'https://github.com/smoitra87/...' ,
    'download_url' : 'https://github.com/smoitra87/.../downloads',
    'author_email' : 'subho@cmu.edu',
    'version' : '0.1',
    'install_requires' : ['nose'],
    'packages' : ['NAME'],
    'scripts' : [],
    'name' : 'NAME'
}

setup(**config)

