from distutils.core import setup, Extension



# config = {
#     'description' : 'My Project',
#     'author' : 'Subhodeep Moitra',
#     'url' : 'https://github.com/smoitra87/...' ,
#     'download_url' : 'https://github.com/smoitra87/.../downloads',
#     'author_email' : 'subho@cmu.edu',
#     'version' : '0.1',
#     'install_requires' : ['nose'],
#     'packages' : ['NAME'],
#     'scripts' : [],
#     'name' : 'NAME'
# }

#setup(**config)


setup(name='S-Systems using ALR',
	version='0.1.0',
	url = 'http://github.com/smoitra87/ssystem',
	author='Subhodeep Moitra',
	author_email='smoitra@cs.cmu.edu',
	license='BSD License',
	description='Solves S-system problems using ALR',
	packages=['ssystem'],
	ext_modules=[Extension('ssystem.cparser',
	['ssystem/keywordReader.c'])])


