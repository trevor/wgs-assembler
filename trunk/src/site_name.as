#
# $Id: site_name.as,v 1.3 2006-01-12 21:34:56 eliv Exp $
#
# This file sets the site name where the assembler is to be built
# Site-specific settings can be made in c_make.as and other makfiles,
# such as AS_UID/uid_transport.as
#

# Default site name is invalid to prompt user to edit this file
ifeq ($(SITE_NAME),)
#  SITE_NAME=CELERA
#  SITE_NAME=TIGR
  SITE_NAME=JCVI
  SITE_INC=/bioinfo/assembly/wgs_include
endif

ifeq ($(SITE_NAME), CELERA)
  $(error Bad SITE_NAME: CELERA. Please edit the file src/site_name.as!!)
endif
