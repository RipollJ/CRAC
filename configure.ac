dnl Process this file with autoconf to produce a configure script.

dnl Copyright © 2010-2013 -- IRB/INSERM
dnl                          (Institut de Recherches en Biothérapie /
dnl                          Institut National de la Santé et de la Recherche
dnl                          Médicale)
dnl                          LIFL/INRIA
dnl                          (Laboratoire d'Informatique Fondamentale de
dnl                          Lille / Institut National de Recherche en
dnl                          Informatique et Automatique)
dnl                          LIRMM/CNRS
dnl                          (Laboratoire d'Informatique, de Robotique et de
dnl                          Microélectronique de Montpellier /
dnl                          Centre National de la Recherche Scientifique)
dnl                          LITIS
dnl                          (Laboratoire d'Informatique, du Traitement de
dnl                          l'Information et des Systèmes).
dnl
dnl Auteurs/Authors:
dnl                   Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>
dnl                   Mikaël SALSON    <mikael.salson@lifl.fr>
dnl                   Thérèse COMMES   <commes@univ-montp2.fr>
dnl                   Éric RIVALS      <eric.rivals@lirmm.fr>
dnl
dnl Programmeurs
dnl      /Progammers:
dnl                   Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>
dnl                   Mikaël SALSON    <mikael.salson@lifl.fr>
dnl  with additional contribution for the packaging of:
dnl                   Alban MANCHERON  <alban.mancheron@lirmm.fr>
dnl
dnl  Contact:         CRAC list   <crac-bugs@lists.gforge.inria.fr>
dnl -------------------------------------------------------------------------
dnl
dnl Ce fichier fait partie du programme CRAC.
dnl
dnl Ce logiciel est régi  par la licence CeCILL  soumise au droit français et
dnl respectant les principes  de diffusion des logiciels libres.  Vous pouvez
dnl utiliser, modifier et/ou redistribuer ce programme sous les conditions de
dnl la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur
dnl le site "http://www.cecill.info".
dnl
dnl En contrepartie de l'accessibilité au code source et des droits de copie,
dnl de modification et de redistribution accordés par cette licence, il n'est
dnl offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
dnl seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le
dnl titulaire des droits patrimoniaux et les concédants successifs.
dnl
dnl À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques
dnl associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au
dnl développement  et à la reproduction du  logiciel par  l'utilisateur étant
dnl donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à
dnl manipuler et qui le réserve donc à des développeurs et des professionnels
dnl avertis  possédant  des  connaissances  informatiques  approfondies.  Les
dnl utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du
dnl logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la
dnl sécurité de leurs systêmes et ou de leurs données et,  plus généralement,
dnl à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.
dnl
dnl Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez
dnl pris connaissance  de la licence CeCILL,  et que vous en avez accepté les
dnl termes.
dnl
dnl -------------------------------------------------------------------------
dnl
dnl This File is part of the CRAC program.
dnl
dnl This software is governed by the CeCILL license under French law and
dnl abiding by the rules of distribution of free software. You can use,
dnl modify and/ or redistribute the software under the terms of the CeCILL
dnl license as circulated by CEA, CNRS and INRIA at the following URL
dnl "http://www.cecill.info".
dnl
dnl As a counterpart to the access to the source code and rights to copy,
dnl modify and redistribute granted by the license, users are provided only
dnl with a limited warranty and the software's author, the holder of the
dnl economic rights, and the successive licensors have only limited
dnl liability.
dnl
dnl In this respect, the user's attention is drawn to the risks associated
dnl with loading, using, modifying and/or developing or reproducing the
dnl software by the user in light of its specific status of free software,
dnl that may mean that it is complicated to manipulate, and that also
dnl therefore means that it is reserved for developers and experienced
dnl professionals having in-depth computer knowledge. Users are therefore
dnl encouraged to load and test the software's suitability as regards their
dnl requirements in conditions enabling the security of their systems and/or
dnl data to be ensured and, more generally, to use and operate it in the same
dnl conditions as regards security.
dnl
dnl The fact that you are presently reading this means that you have had
dnl knowledge of the CeCILL license and that you accept its terms.

dnl $Id$

dnl Process this file with autoconf to produce a configure script.

AC_PREREQ([2.64])
m4_define([V_MAJOR], [1])
m4_define([V_MINOR], [5])
m4_define([V_MICRO], [0])

AC_INIT(
  [crac],
  [V_MAJOR.V_MINOR.V_MICRO],
  [crac-bugs@lists.gforge.inria.fr],
  [crac],
  [http://crac.gforge.inria.fr/])

PACKAGE_DESCRIPTION="Crac is a tool to analyse RNA-Seq data provided by NGS"
AC_SUBST([PACKAGE_DESCRIPTION])

AC_CONFIG_SRCDIR([src/crac.cpp])
AC_REVISION([$Revision$])
AC_COPYRIGHT([
*******************************************************************************
*                                                                             *
*  Copyright © 2010-2013 -- IRB/INSERM                                        *
*                           (Institut de Recherches en Biothérapie /          *
*                           Institut National de la Santé et de la Recherche  *
*                           Médicale)                                         *
*                           LIFL/INRIA                                        *
*                           (Laboratoire d'Informatique Fondamentale de       *
*                           Lille / Institut National de Recherche en         *
*                           Informatique et Automatique)                      *
*                           LIRMM/CNRS                                        *
*                           (Laboratoire d'Informatique, de Robotique et de   *
*                           Microélectronique de Montpellier /                *
*                           Centre National de la Recherche Scientifique)     *
*                           LITIS                                             *
*                           (Laboratoire d'Informatique, du Traitement de     *
*                           l'Information et des Systèmes).                   *
*                                                                             *
*                                                                             *
*  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *
*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *
*                   Thérèse COMMES   <commesd@univ-montp2.fr>                 *
*                   Éric RIVALS      <eric.rivals@lirmm.fr>                   *
*                                                                             *
*  Programmeurs                                                               *
*      /Progammers: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *
*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *
*  with additional contribution of:                                           *
*                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               *
*                                                                             *
*  Contact:         CRAC list        <crac-bugs@lists.gforge.inria.fr>        *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  Ce fichier fait partie du programme CRAC.                                  *
*                                                                             *
*  Crac est un outil d'analyse de données de RNA-Seq provenant des NGS.       *
*                                                                             *
*  Ce logiciel est régi  par la licence CeCILL  soumise au droit français et  *
*  respectant les principes  de diffusion des logiciels libres.  Vous pouvez  *
*  utiliser, modifier et/ou redistribuer ce programme sous les conditions de  *
*  la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur  *
*  le site "http://www.cecill.info".                                          *
*                                                                             *
*  En contrepartie de l'accessibilité au code source et des droits de copie,  *
*  de modification et de redistribution accordés par cette licence, il n'est  *
*  offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,  *
*  seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le  *
*  titulaire des droits patrimoniaux et les concédants successifs.            *
*                                                                             *
*  À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques  *
*  associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au  *
*  développement  et à la reproduction du  logiciel par  l'utilisateur étant  *
*  donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à  *
*  manipuler et qui le réserve donc à des développeurs et des professionnels  *
*  avertis  possédant  des  connaissances  informatiques  approfondies.  Les  *
*  utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du  *
*  logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la  *
*  sécurité de leurs systêmes et ou de leurs données et,  plus généralement,  *
*  à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.         *
*                                                                             *
*  Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez  *
*  pris connaissance  de la licence CeCILL,  et que vous en avez accepté les  *
*  termes.                                                                    *
*                                                                             *
*  -------------------------------------------------------------------------  *
*                                                                             *
*  This File is part of the CRAC program.                                     *
*                                                                             *
*  Crac is a tool to analyse RNA-Seq data provided by NGS.                    *
*                                                                             *
*  This software is governed by the CeCILL license under French law and       *
*  abiding by the rules of distribution of free software. You can use,        *
*  modify and/ or redistribute the software under the terms of the CeCILL     *
*  license as circulated by CEA, CNRS and INRIA at the following URL          *
*  "http://www.cecill.info".                                                  *
*                                                                             *
*  As a counterpart to the access to the source code and rights to copy,      *
*  modify and redistribute granted by the license, users are provided only    *
*  with a limited warranty and the software's author, the holder of the       *
*  economic rights, and the successive licensors have only limited            *
*  liability.                                                                 *
*                                                                             *
*  In this respect, the user's attention is drawn to the risks associated     *
*  with loading, using, modifying and/or developing or reproducing the        *
*  software by the user in light of its specific status of free software,     *
*  that may mean that it is complicated to manipulate, and that also          *
*  therefore means that it is reserved for developers and experienced         *
*  professionals having in-depth computer knowledge. Users are therefore      *
*  encouraged to load and test the software's suitability as regards their    *
*  requirements in conditions enabling the security of their systems and/or   *
*  data to be ensured and, more generally, to use and operate it in the same  *
*  conditions as regards security.                                            *
*                                                                             *
*  The fact that you are presently reading this means that you have had       *
*  knowledge of the CeCILL license and that you accept its terms.             *
*                                                                             *
*******************************************************************************
])

AS_IF([test x"$CFLAGS" == "x"], [CFLAGS=""])
AS_IF([test x"$CPPFLAGS" == "x"], [CPPFLAGS=""])
AS_IF([test x"$CXXFLAGS" == "x"], [CXXFLAGS=""])
AS_IF([test x"$LDFLAGS" == "x"], [LDFLAGS=""])

AC_CONFIG_AUX_DIR([config])

dnl create a config.h file (Automake will add -DHAVE_CONFIG_H)
AC_CONFIG_HEADERS([config/config.h])

dnl Every other copy of the package version number gets its value from here
AM_INIT_AUTOMAKE()

AC_CANONICAL_HOST
AC_LANG([C++])

AC_CONFIG_MACRO_DIR([config/m4])

dnl Checks for programs.
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_CC
AC_PROG_CXX
AC_PROG_CPP
AC_PROG_MAKE_SET
AC_PROG_LN_S

AC_CHECK_PROG([TEXI2DVI], [texi2dvi], [texi2dvi])
AC_CHECK_PROG([PDFTEX], [pdftex], [pdftex])
AC_CHECK_PROG([KPSEWHICH], [kpsewhich], [kpsewhich])
BUILD_DOC=false
AS_IF([test "x$TEXI2DVI" != "xtexi2dvi" || test "x$PDFTEX" != "xpdftex"],
      [BUILD_DOC=false
       AC_MSG_WARN([The 'texi2dvi' or 'pdftex' program could not be found. This isn't critical, but it means that you won't be able to create the documentation and thus to make releases.])])

AS_IF([test "x$KPSEWHICH" != "xkpsewhich" || test "$(kpsewhich  -show-path=fmt)" == "/nonesuch"],
      [BUILD_DOC=false
       AC_MSG_WARN([The 'texi2dvi' or 'pdftex' is missing a correct latex installation. This isn't critical, but it means that you won't be able to create the documentation and thus to make releases.])])
AM_CONDITIONAL([BUILD_DOC], [$BUILD_DOC])

AC_SUBST([VERSION])

dnl Disable shared library construction by default (must appear before LT_INIT)
AC_DISABLE_SHARED

LT_INIT
LT_PATH_LD

ISODATE=`date +%Y-%m-%d`
AC_SUBST([ISODATE])

ARCH=`uname -m | sed s/i.86/i386/ | sed s/x86_64/amd64/`
AC_SUBST([ARCH])

AC_PATH_PROG([DEB_BUILDER], [dpkg])
AC_PATH_PROG([RPM_BUILDER], [rpmbuild])

# Checks for libraries.
AC_CHECK_LIB([m], [sqrt])
AC_CHECK_LIB([pthread], [pthread_join])

# Checks for header files.
AC_FUNC_ALLOCA
AC_HEADER_MAJOR
AC_CHECK_HEADERS([\
  fcntl.h inttypes.h libintl.h limits.h locale.h \
  memory.h stdint.h stdlib.h string.h strings.h \
  sys/file.h sys/param.h sys/time.h unistd.h utime.h \
])

# Checks for typedefs, structures, and compiler characteristics.
AC_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_CHECK_MEMBERS([struct stat.st_blksize])
AC_TYPE_UINT32_T
AC_TYPE_UINT8_T
AC_CHECK_TYPES([ptrdiff_t])

# Checks for library functions.
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([\
  floor gettimeofday memchr memmove memset \
  mkfifo pow select sqrt strrchr strstr strtol \
])


CPPFLAGS="$CPPFLAGS -Wall -ansi -funroll-loops"

dnl Add option --enable-debug (default to --disable-debug)
AC_ARG_ENABLE([debug],
              [AS_HELP_STRING([--enable-debug],
                              [enable debugging informations])],
              [CPPFLAGS="$CPPFLAGS -g -O0"],
              [CPPFLAGS="$CPPFLAGS -O4 -DNDEBUG"])

dnl Add option --without-included-GkArrays (default --with-included-GkArrays=check)
AC_ARG_WITH([included-GkArrays],
            [AS_HELP_STRING([--without-included-GkArrays],
                            [don't use the GkArrays library included here])],
            [],
            [with_included_GkArrays="check"])

dnl Add option --with-libGkArrays-prefix=PATH (default --without-libGkArrays-prefix)
AC_ARG_WITH([libGkArrays-prefix],
            [AS_HELP_STRING([--with-libGkArrays-prefix=DIR],
                            [search for libGkArrays in DIR/include/libGkArrays and DIR/lib])],
            [LDFLAGS="$LDFLAGS -L$with_libGkArrays_prefix/lib"
             CPPFLAGS="$CPPFLAGS -I$with_libGkArrays_prefix/include/libGkArrays"],
            [with_libGkArrays_prefix=""])

dnl Check if we need included GkArrays library.
AS_IF([test "x$with_included_GkArrays" == "xcheck"],
      [AS_IF([test "x$with_libGkArrays_prefix" != "x"],
             [with_included_GkArrays="no"],
             [with_included_GkArrays="yes"])])

AS_IF([test "x$with_included_GkArrays" == "xyes"],
      [AS_IF([test "x$with_libGkArrays_prefix" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-GkArrays and --with-libGkArrays-prefix options])])])

dnl Add option --without-included-ProgressBar (default --with-included-ProgressBar=check)
AC_ARG_WITH([included-ProgressBar],
            [AS_HELP_STRING([--without-included-ProgressBar],
                            [don't use the ProgressBar library included here])],
            [],
            [with_included_ProgressBar="check"])

dnl Add option --with-libProgressBar-prefix=PATH (default --without-libProgressBar-prefix)
AC_ARG_WITH([libProgressBar-prefix],
            [AS_HELP_STRING([--with-libProgressBar-prefix=DIR],
                            [search for libProgressBar in DIR/include/libProgressBar and DIR/lib])],
            [LDFLAGS="$LDFLAGS -L$with_libProgressBar_prefix/lib"
             CPPFLAGS="$CPPFLAGS -I$with_libProgressBar_prefix/include/libProgressBar"],
            [with_libProgressBar_prefix=""])

dnl Check if we need included ProgressBar library.
PB_OK=0
AS_IF([test "x$with_included_ProgressBar" == "xcheck"],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [with_included_ProgressBar="no"],
             [with_included_ProgressBar="$with_included_GkArrays"])])

AS_IF([test "x$with_included_ProgressBar" == "xyes"],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [AC_MSG_FAILURE([You can't use both --with-included-ProgressBar and --with-libProgressBar-prefix options])],
             [AS_IF([test "x$with_included_GkArrays" != "xyes"],
                    [AC_MSG_FAILURE([You can't use both --with-included-ProgressBar and --without-included-GkArrays options])],
                    [PB_OK=1])
             ])],
      [AS_IF([test "x$with_libProgressBar_prefix" != "x"],
             [PB_OK=1
              LDFLAGS="$LDFLAGS -L$with_libProgressBar_prefix -L$with_libProgressBar_prefix/lib"
              CPPFLAGS="$CPPFLAGS -I$with_libProgressBar_prefix -I$with_libProgressBar_prefix/include"],
             [AC_CHECK_LIB([ProgressBar],
                           [libProgressBarVersion],
                           [PB_OK=1
                            ProgressBar_LIBS=-lProgressBar])])]) dnl Override the default behavior of AC_CHECK_LIB (we don't want -lProgressBar in LIBS).

AS_IF([test "$PB_OK" == "1"], [AC_DEFINE([HAVE_LIBPROGRESSBAR], [1], [Define to 1 in order to Use libProgressBar])])

dnl Set INCLUDED_PROGRESSBAR conditionnal for automake
AM_CONDITIONAL([INCLUDED_PROGRESSBAR], [test "x$with_included_ProgressBar" == "xyes"])

AC_CHECK_LIB([z], [gzread], ,
             [AC_MSG_ERROR([zlib not found, see http://www.zlib.net])])

AS_IF([test "x$with_included_GkArrays" != "xyes"],
      [AC_CHECK_LIB([GkArrays],
                    [libGkArraysVersion],
                    [:],dnl Override the default behavior of AC_CHECK_LIB (we don't want -lGkArrays in LIBS).
                    [AC_MSG_FAILURE([[--without-included-GkArrays was given, but test for libGkArrays failed]])])],
      [AC_CONFIG_SUBDIRS([src/libGkArrays])])

dnl Set INCLUDED_GKARRAYS conditionnal for automake
AM_CONDITIONAL([INCLUDED_GKARRAYS], [test "x$with_included_GkArrays" == xyes])

dnl Forward options given to current configure script to libGkArrays configure script if needed
AS_IF([test "x$enable_shared" == "xyes" ],
      [AS_VAR_APPEND([ac_configure_args], [" --enable-shared"])],
      [AS_VAR_APPEND([ac_configure_args], [" --disable-shared"])])
AS_IF([test "x$enable_static" == "xyes" ],
      [AS_VAR_APPEND([ac_configure_args], [" --enable-static"])],
      [AS_VAR_APPEND([ac_configure_args], [" --disable-static"])])

dnl Forward FLAGS to libGkArrays/libProgressBar configure scripts if needed
AS_VAR_APPEND([ac_configure_args], [" CPPFLAGS=\"$CPPFLAGS\""])
AS_VAR_APPEND([ac_configure_args], [" CFLAGS=\"$CFLAGS\""])
AS_VAR_APPEND([ac_configure_args], [" CXXFLAGS=\"$CXXFLAGS\""])

LDFLAGS="$LDFLAGS $ProgressBar_LIBS"
CPPFLAGS="$CPPFLAGS $ProgressBar_CFLAGS"

AC_CONFIG_FILES([\
  control crac.spec                             \
  Makefile src/Makefile scripts/Makefile        \
  src/libReadsInfo/Makefile src/libSSA/Makefile \
  src/libgzstream/Makefile                      \
  doc/Makefile doc/package.texi                 \
  man/Makefile man/crac.1                       \
  scripts/crac-client                           \
  test_f/Makefile test_f/bug/Makefile           \
  test_f/causes/Makefile test_f/cigar/Makefile  \
], [ \
  test -f scripts/crac-client && chmod 755 scripts/crac-client ; \
])

AC_OUTPUT
