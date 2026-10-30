# Meson build for the pvcraze *Python package* only.
#
# The C library and the GNU Octave package are built with autotools
# (./configure && make). This file compiles the library sources directly
# into the Python extension module `pvcraze._core`, so `pip install .`
# does not need an installed libpvcraze.
project(
  'pvcraze',
  'c',
  version: files('VERSION'),
  license: 'GPL-3.0-only',
  meson_version: '>=1.1',
  default_options: ['buildtype=release', 'b_ndebug=if-release', 'warning_level=1'],
)

cc = meson.get_compiler('c')
py = import('python').find_installation(pure: false)

m_dep = cc.find_library('m', required: false)
gsl_dep = dependency('gsl')
lapack_dep = dependency('lapack', required: false)
if not lapack_dep.found()
  # OpenBLAS bundles LAPACK (e.g. MSYS2 on Windows ships no lapack.pc)
  lapack_dep = dependency('openblas', required: false)
endif
if not lapack_dep.found()
  lapack_dep = cc.find_library('lapack')
endif
numpy_dep = dependency('numpy')

core_c_args = ['-DNPY_NO_DEPRECATED_API=NPY_1_20_API_VERSION']
# Python's Windows headers only define MS_WIN64 for MSVC; mingw-w64 needs it
# explicitly or Py_ssize_t and friends come out 32-bit.
if host_machine.system() == 'windows' and cc.get_id() == 'gcc' and host_machine.cpu_family() == 'x86_64'
  core_c_args += '-DMS_WIN64'
endif

lib_sources = [
  'emalloc.c', 'curve.c', 'fitmodels.c', 'curiosity.c', 'Error.c',
  'quickselect.c', 'curvemodel.c', 'spline_iv.c', 'splinefit.c', 'KH_iv.c',
  'PC_iv.c', 'diode_iv.c', 'bezier2_iv.c', 'bezier3_iv.c', 'pwl_iv.c',
  'apwl_iv.c', 'neldermead.c', 'diffevol.c', 'furthestsampling.c',
  'lib-pv-craze.c', 'tictoc.c',
]
model_headers = [
  'bezier2_iv.h', 'KH_iv.h', 'spline_iv.h', 'bezier3_iv.h', 'diode_iv.h',
  'PC_iv.h', 'pwl_iv.h', 'apwl_iv.h',
]
other_headers = [
  'emalloc.h', 'curve.h', 'fitmodels.h', 'curiosity.h', 'Error.h',
  'quickselect.h', 'curvemodel.h', 'splinefit.h', 'neldermead.h',
  'diffevol.h', 'furthestsampling.h', 'lib-pv-craze.h',
  'lib-pv-craze-structs.h', 'tictoc.h',
]

# Copy the library sources into the build directory. The compiler searches
# the directory of the including file first, so compiling straight from src/
# would pick up stale generated headers left there by an in-tree autotools
# build instead of the ones generated below.
copied_sources = []
foreach f : lib_sources
  copied_sources += configure_file(input: 'src' / f, output: f, copy: true)
endforeach
model_header_files = []
foreach f : model_headers
  model_header_files += configure_file(input: 'src' / f, output: f, copy: true)
endforeach
foreach f : other_headers
  configure_file(input: 'src' / f, output: f, copy: true)
endforeach

bash = find_program('bash')
srcdir = meson.current_source_dir() / 'src'
model_gen = custom_target(
  'models',
  input: ['src/register_models.sh'] + model_header_files,
  output: ['models.h', 'model_table.h'],
  command: [bash, '@INPUT0@', srcdir] + model_headers,
)
flag_gen = custom_target(
  'ivanalyzeflags',
  input: ['src/gen_ivanalyzeflags.sh', 'src/curve.c'],
  output: ['ivanalyzeflags.h', 'ivanalyzemessages.h'],
  command: [bash, '@INPUT0@', srcdir, 'curve.c'],
)

py.extension_module(
  '_core',
  ['python/_core.c'] + copied_sources,
  model_gen,
  flag_gen,
  dependencies: [numpy_dep, gsl_dep, lapack_dep, m_dep],
  c_args: core_c_args,
  gnu_symbol_visibility: 'hidden',
  install: true,
  subdir: 'pvcraze',
)

py.install_sources(
  [
    'python/pvcraze/__init__.py',
    'python/pvcraze/_analysis.py',
    'python/pvcraze/_arrays.py',
    'python/pvcraze/_esp.py',
    'python/pvcraze/_flags.py',
    'python/pvcraze/_models.py',
    'python/pvcraze/_verbosity.py',
    'python/pvcraze/diode_gt.py',
    'python/pvcraze/iec60891.py',
    'python/pvcraze/_core.pyi',
    'python/pvcraze/py.typed',
  ],
  subdir: 'pvcraze',
)
