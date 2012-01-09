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

module1 = Extension('spam',sources=['ssystem/spammodule.c'])

setup(name='spam',
	version=1.0,
	description='Just a spam module',
	ext_modules = [module1])


