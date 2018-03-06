from distutils.core import setup, Extension

setup (
   name = 'mergeseq',
   version = '0.1',
   description = 'Do sequencing merging with C code.',
   ext_modules = [
      Extension('mergeseq', sources = ['mergeseqmodule.c']),
   ]
)

