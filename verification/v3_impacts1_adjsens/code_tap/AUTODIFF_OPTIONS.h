#ifndef AUTODIFF_OPTIONS_H
#define AUTODIFF_OPTIONS_H
#include "PACKAGES_CONFIG.h"
#include "CPP_OPTIONS.h"

CBOP
C !ROUTINE: AUTODIFF_OPTIONS.h
C !INTERFACE:
C #include "AUTODIFF_OPTIONS.h"

C !DESCRIPTION:
C *==================================================================*
C | CPP options file for AutoDiff (autodiff) package:
C | Control which optional features to compile in this package code.
C *==================================================================*
CEOP

#ifdef ALLOW_AUTODIFF
#ifdef ECCO_CPPOPTIONS_H

C-- When multi-package option-file ECCO_CPPOPTIONS.h is used (directly included
C    in CPP_OPTIONS.h), this option file is left empty since all options that
C   are specific to this package are assumed to be set in ECCO_CPPOPTIONS.h

#else /* ndef ECCO_CPPOPTIONS_H */
C   ==================================================================
C-- Package-specific Options & Macros go here

C o Include/exclude code in order to automatically differentiate MITgcm code
C   using TAF (Transformation of Algorithms in Fortran, http://www.FastOpt.de)
C   or using TAMC (Tangent Linear & Adjoint Model Compiler, needs both defined):
#undef ALLOW_AUTODIFF_TAMC
#undef AUTODIFF_TAMC_COMPATIBILITY

C       >>> Checkpointing as handled by TAMC
#undef ALLOW_TAMC_CHECKPOINTING

C       >>> Extract adjoint state
#define ALLOW_AUTODIFF_MONITOR
C       >>> and DYNVARS_DIAG adjoint state
#undef ALLOW_AUTODIFF_MONITOR_DIAG

C       >>> DO 2-level checkpointing instead of 3-level
#undef AUTODIFF_2_LEVEL_CHECKPOINT

C extend to 4-level checkpointing
#undef AUTODIFF_4_LEVEL_CHECKPOINT

C o use divided adjoint to split adjoint computations
#undef ALLOW_DIVIDED_ADJOINT

C o This flag is incredibly useful as it reduces the number of
C   tape-files on the disc. Maybe it should even be the default.
#define ALLOW_AUTODIFF_WHTAPEIO
C   and related to above:
#undef ALLOW_INIT_WHTAPEIO

C o use standard MDSFINDUINTS instead of local pkg/autodiff version for
C   WHTAPEIO code I/O.
C   Note: comment out the #define below (instead of having an #undef) to
C   enable to set this Option in CPP command line (from the optfile)
c#define AUTODIFF_USE_MDSFINDUNITS

C o use the deprecated autodiff_store/restore method where multiple fields
C   are collected in a single buffer field array before storing to tape.
C   This functionality has been replaced by WHTAPEIO method (see above).
C   Might still be used for OBCS since WHTAPEIO does not support OBCS fields.
#undef AUTODIFF_USE_STORE_RESTORE
#undef AUTODIFF_USE_STORE_RESTORE_OBCS

C o allow using viscFacInAd to recompute viscosities in AD
#undef AUTODIFF_ALLOW_VISCFACADJ

C o To remove part of MOM_CALC_VISC (better name would be: MOM_DISABLE_*)
#undef AUTODIFF_DISABLE_LEITH
#undef AUTODIFF_DISABLE_REYNOLDS_SCALE

C o for output of AD-variables (ALLOW_AUTODIFF_MONITOR), specific code (e.g.,
C   in addummy_in_stepping.F) relies on adexch_uv_xy_rs and adexch_xy_rs S/R
C   which might not always be generated by TAF (e.g., when controls do not
C   include any 2D forcing field). In those cases, defining this cpp-option
C   allows to circumvent this missing code issue.
#undef AUTODIFF_EXCLUDE_ADEXCH_RS

C   ==================================================================
#endif /* ndef ECCO_CPPOPTIONS_H */
#endif /* ALLOW_AUTODIFF */
#endif /* AUTODIFF_OPTIONS_H */
