dnl                                               -*- Autoconf -*-
dnl Process this file with autoconf to produce a configure script.

AC_CONFIG_MACRO_DIR([m4])
AC_INIT
AM_PATH_PYTHON
AC_CONFIG_SRCDIR([$AEPHEM_NAME])
AC_CONFIG_HEADER([src/config.h])

AEPHEM_NAME=aephem
AC_SUBST(AEPHEM_NAME)
AEPHEM_RELEASE=canopus
AC_SUBST(AEPHEM_RELEASE)
AEPHEM_VERSION=2:0:0
#              | | |
#       +------+ | +---+
#       |        |     |
#    current:revision:age
#       |        |     |
#       |        |     +- increment if interfaces have been added
#       |        |        set to zero if interfaces have been removed
#       |        |        or changed
#       |        +- increment if source code has changed
#       |           set to zero if current is incremented
#       +- increment if interfaces have been added, removed or changed
AC_SUBST(AEPHEM_VERSION)

AC_DEFINE([AE_MAJOR], 2, [ The major version number ])
AC_DEFINE([AE_MINOR], 0, [ The minor version number ])
AC_DEFINE([AE_REVISION], 0, [ The revision number ])
AC_DEFINE([AE_VERSION_SUFFIX], [0],
          [ Other qualifiers on the version number triplet (if any) ])


AM_INIT_AUTOMAKE($AEPHEM_NAME, $AEPHEM_VERSION, no-define)

dnl Checks for programs.
AC_PROG_CC
LT_INIT
AM_SANITY_CHECK

dnl Checks for header files.
AC_HEADER_STDC
AC_CHECK_HEADERS([math.h stdio.h stdlib.h time.h ctype.h string.h], ,
                 AC_MSG_ERROR([*** missing standard header!]))

dnl Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_CONST

dnl Checks for library functions.
AC_FUNC_REALLOC
AC_FUNC_STAT

dnl Bindings
make_pybindings="yes"
AC_ARG_ENABLE(python, AS_HELP_STRING([--disable-python],
                                    [do not build python bindings (aephem)]),
              [case "{$enableval}" in
                 no) make_pybindings="no" ;;
                 *) make_pybindings="yes" ;
               esac
              ])

AC_ARG_WITH([python],
            AS_HELP_STRING([--with-python=PATH],
                           [use python interpreter in PATH [autodetect]]),
            [case "${withval}" in
               no) have_python="no" ;;
               yes) python_path= ; have_python= ;;
               *) python_path="${withval}"; have_python= ;;
             esac],
            [python_path=; have_python=])

AC_ARG_WITH([python-module-dir],
            AS_HELP_STRING([--with-python-module-dir=PATH],
                           [install the python bindings in PATH [autodetect]]),
            [case "${withval}" in
               no) python_mod_path= ;;
               *) python_mod_path="${withval}"
             esac],
            [python_mod_path= ])


AC_MSG_CHECKING([whether to include python bindings])
if test "x$make_pybindings" = "xyes"; then
  AC_PROG_CC_C99
  if test "x$ac_cv_prog_cc_c99" = "xno"; then
    make_pybindings="no";
    AC_MSG_RESULT([no (requires C99 compiler)])
  else
    AC_MSG_RESULT([$make_pybindings])
  fi
fi

if test "x$make_pybindings" = "xyes"; then
  echo
  echo "--- Configuring python bindings ---"

  first_vers=2.3
  final_vers=2.7
  python_list="python python2"
  AC_CHECK_PROGS([SEQ], [seq], [not found])
  if test "x$SEQ" == "xnot found"; then
    AC_CHECK_PROGS([JOT], [jot], [not found])
    if test "x$JOT" != "xnot found"; then
      python_list="$python_list `$JOT -w 'python%.1f' - \
                                 $final_vers $first_vers -0.1`"
    fi
  else
    python_list="$python_list `$SEQ -f 'python%.1f' \
                               $final_vers -0.1 $first_vers`"
  fi

  dnl Try to find a suitably new version of python.
  if test "x$python_path" == "x"; then
    AC_MSG_CHECKING([for python version >= $first_vers])
    PYTHON="not found"
    for p in $python_list; do
      _AS_PATH_WALK([$PATH],
                    [for e in '' $ac_executable_extensions; do
                       if AS_EXECUTABLE_P(["$as_dir/$p$e"]); then
                         AM_PYTHON_CHECK_VERSION(["$as_dir/$p$e"],
                                                 [$first_vers],
                                               [PYTHON="$as_dir/$p$e"; break 3])
                       fi
                     done])
    done
    AC_MSG_RESULT([$PYTHON])
  else
    AC_MSG_CHECKING([whether $python_path >= $first_vers])
    AM_PYTHON_CHECK_VERSION([$python_path], [$first_vers],
                             AC_MSG_RESULT([yes]), PYTHON=$python_path,
                             AC_MSG_RESULT([no]), PYTHON="not found"])
  fi

  if test "x$PYTHON" = "xnot found"; then
    have_python="no"
    PYTHON=
  fi

  AC_SUBST([PYTHON])
fi

if test "x${have_python}" != "xno"; then
  AC_MSG_CHECKING([$PYTHON version])
  PYTHON_VERSION=`$PYTHON -c "import sys; print sys.version[[:3]]"`
  AC_MSG_RESULT([$PYTHON_VERSION])

  dnl Figure out python CPPFLAGS and LIBS.
  if test -x $PYTHON-config; then
    if test -n "$python_path"; then
      python_exec_prefix=`$PYTHON-config --exec-prefix`
      PYTHON_LIBS="-L${python_exec_prefix}/lib "
    else
      PYTHON_LIBS=
    fi
    PYTHON_CPPFLAGS=`$PYTHON-config --includes 2>/dev/null`
    PYTHON_LIBS="${PYTHON_LIBS} `$PYTHON-config --ldflags 2>/dev/null`"
  else
    python_prefix=`$PYTHON -c "import sys; print sys.prefix"`
    python_exec_prefix=`$PYTHON -c "import sys; pring sys.exec_prefix"`
    python_libdir=`$PYTHON -c "from distutils import sysconfig; \
                   print sysconfig.get_config_var('LIBDIR')"`
    python_syslibs=`$PYTHON -c "from distutils import sysconfig; \
                    print sysconfig.get_config_var('SYSLIBS')"`
    python_shlibs=`$PYTHON -c "from distutils import sysconfig; \
                   print sysconfig.get_config_var('SHLIBS')"`
    python_shlibs=`$PYTHON -c "from distutils import sysconfig; \
                   print sysconfig.get_config_var('SHLIBS')"`
    python_modlibs=`$PYTHON -c "from distutils import sysconfig; \
                    print sysconfig.get_config_var('MODLIBS')"`
    PYTHON_CPPFLAGS="-I${python_prefix}/include/python${PYTHON_VERSION} -I${python_exec_prefix}/include/python${PYTHON_VERSION}"
    if test -n "$python_path"; then
      PYTHON_LIBS="-L${python_libdir} "
    else
      PYTHON_LIBS=
    fi
  fi
  PYTHON_LIBS="${PYTHON_LIBS} $python_syslibs $python_shlibs $python_modlibs -lpython${PYTHON_VERSION}"

  AC_MSG_CHECKING([python includes])
  AC_MSG_RESULT([$PYTHON_CPPFLAGS])
  AC_SUBST([PYTHON_CPPFLAGS])
  AC_MSG_CHECKING([python libraries])
  AC_MSG_RESULT([$PYTHON_LIBS])
  AC_SUBST([PYTHON_LIBS])

  dnl Check for C header.
  saved_CPPFLAGS=${CPPFLAGS}
  CPPFLAGS="${CPPFLAGS} ${PYTHON_CPPFLAGS}"
  AC_CHECK_HEADERS(Python.h,,[have_python="no"])
  CPPFLAGS=${saved_CPPFLAGS}
fi

if test "x${have_python}" != "xno"; then
  dnl Calculate the exec prefix.
  pyexec_prefix=$exec_prefix
  test "x$pyexec_prefix" = xNONE && pyexec_prefix=$prefix
  test "x$pyexec_prefix" = xNONE && pyexec_prefix=$ac_default_prefix

  dnl Calculate the extension module directory
  AC_MSG_CHECKING([python extension module directory])
  if test "x${python_mod_path}" = "x"; then
    pythondir=`$PYTHON -c "from distutils import sysconfig; print sysconfig.get_python_lib(1,0,prefix='${pyexec_prefix}')" 2>/dev/null || echo "${pyexec_prefix}/lib/python${PYTHON_VERSION}/site-packages"`
  else
    pythondir=$python_mod_path
  fi

  AC_SUBST([pythondir])
  AC_MSG_RESULT([$pythondir])

  saved_CPPFLAGS=${CPPFLAGS}
  CPPFLAGS="${CPPFLAGS} ${PYTHON_CPPFLAGS}"
  saved_LIBS=$LIBS
  LIBS="$LIBS ${PYTHON_LIBS}"
  
  dnl Try to compile a module
  AC_MSG_CHECKING([if we can compile a simple python extension module])
  AC_TRY_LINK_FUNC([Py_Initialize], [AC_MSG_RESULT([yes])],
                   [AC_MSG_RESULT([no]); have_python="no"])
  CPPFLAGS=$saved_CPPFLAGS
  LIBS=$saved_LIBS
fi

AM_CONDITIONAL(MAKE_PYBINDINGS, [test "x$make_pybindings" = "xyes"])

if test "x$make_pybindings" != "xno"; then
  BINDINGS_BUILT="${BINDINGS_BUILT} python"
else
  BINDINGS_LEFT="${BINDINGS_LEFT} python"
fi


AC_CONFIG_FILES([Makefile src/Makefile demo/Makefile
                 bindings/Makefile bindings/python/Makefile])
AC_OUTPUT
