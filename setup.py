from setuptools import setup
import versioneer
from setuptools import find_packages
from glob import glob
from os.path import basename
from os.path import dirname
from os.path import join
from os.path import splitext


# A list of authors and their email addresses.
authors=("Antonia Mey <antonia.mey@gmail.com, "
         "Jordi Juárez Jiménez <jordi.juarez@gmail.com>")

# Run the setup.
setup(name='StructureSelector',
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass(),
          description='StructureSelector makes it easy to select the correct PdB to Dock against or run MD with',
          author=authors,
          url='https://github.com/michellab/StructureSelector',
          license='MIT',
          packages=find_packages('src'),
          package_dir={'': 'src'},
          py_modules=[splitext(basename(path))[0] for path in glob('src/*.py')],
          include_package_data=True,
          zip_safe=False,
          classifiers=[
            'Development Status :: 2 - Pre-Alpha',
            'Environment :: Console',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Natural Language :: English',
            'Operating System :: MacOS :: MacOS X',
            'Operating System :: POSIX :: Linux',
            'Programming Language :: C',
            'Programming Language :: Cython',
            'Programming Language :: Python :: 3.7',
            'Topic :: Scientific/Engineering :: Bio-Informatics',
            'Topic :: Scientific/Engineering :: Chemistry',
            'Topic :: Scientific/Engineering :: Mathematics',
            'Topic :: Scientific/Engineering :: Physics'
        ],
          keywords=[
            'clustering',
            'RMSD',
            'Diffucsion Map',
            'Protein database'
        ],
          install_requires=[
        'numpy', 'scipy', 'scikit-learn', 'matplotlib'
        ]
)
