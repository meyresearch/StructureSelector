from setuptools import setup
import versioneer


# A list of authors and their email addresses.
authors=("Antonia Mey <antonia.mey@gmail.com, "
         "Jordi Juárez Jiménez <jordi.juarez@gmail.com>")

# Run the setup.
try:
    setup(name='StructureSelector',
          version=versioneer.get_version(),
          cmdclass=versioneer.get_cmdclass(),
          description='StructureSelector makes it easy to select the correct PdB to Dock against or run MD with',
          author=authors,
          url='https://github.com/michellab/StructureSelector',
          license='MIT',
          packages=find_packages(),
          zip_safe=True,
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
        ]
    )

