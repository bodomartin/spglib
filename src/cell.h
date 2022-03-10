/* Copyright (C) 2008 Atsushi Togo */
/* All rights reserved. */

/* This file is part of spglib. */

/* Redistribution and use in source and binary forms, with or without */
/* modification, are permitted provided that the following conditions */
/* are met: */

/* * Redistributions of source code must retain the above copyright */
/*   notice, this list of conditions and the following disclaimer. */

/* * Redistributions in binary form must reproduce the above copyright */
/*   notice, this list of conditions and the following disclaimer in */
/*   the documentation and/or other materials provided with the */
/*   distribution. */

/* * Neither the name of the spglib project nor the names of its */
/*   contributors may be used to endorse or promote products derived */
/*   from this software without specific prior written permission. */

/* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS */
/* "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT */
/* LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS */
/* FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE */
/* COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, */
/* INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, */
/* BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; */
/* LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER */
/* CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT */
/* LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN */
/* ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE */
/* POSSIBILITY OF SUCH DAMAGE. */

#ifndef __cell_H__
#define __cell_H__

#include "mathfunc.h"

typedef struct {
  int size;
  int aperiodic_axis;
  double (*lattice)[3]; /* 3x3 matrix */
  int *types;
  double (*position)[3];
} Cell;

Cell *cel_alloc_cell(const int size);
void cel_free_cell(Cell * cell);
void cel_set_cell(Cell * cell,
                  SPGCONST double lattice[3][3],
                  SPGCONST double position[][3],
                  const int types[]);
void cel_set_layer_cell(Cell * cell,
                        SPGCONST double lattice[3][3],
                        SPGCONST double position[][3],
                        const int types[],
                        const int aperiodic_axis);
Cell * cel_copy_cell(const Cell * cell);
int cel_is_overlap(const double a[3],
                   const double b[3],
                   SPGCONST double lattice[3][3],
                   const double symprec);
int cel_is_overlap_with_same_type(const double a[3],
                                  const double b[3],
                                  const int type_a,
                                  const int type_b,
                                  SPGCONST double lattice[3][3],
                                  const double symprec);
int cel_any_overlap(const Cell * cell,
                    const double symprec);
int cel_any_overlap_with_same_type(const Cell * cell,
                                   const double symprec);
Cell * cel_trim_cell(int * mapping_table,
                     SPGCONST double trimmed_lattice[3][3],
                     const Cell * cell,
                     const double symprec);
int cel_layer_is_overlap(const double a[3],
                         const double b[3],
                         SPGCONST double lattice[3][3],
                         const int periodic_axes[2],
                         const double symprec);
int cel_layer_is_overlap_with_same_type(const double a[3],
                                        const double b[3],
                                        const int type_a,
                                        const int type_b,
                                        SPGCONST double lattice[3][3],
                                        const int periodic_axes[2],
                                        const double symprec);
int cel_layer_any_overlap_with_same_type(const Cell * cell,
                                         const int periodic_axes[2],
                                         const double symprec);

#endif
