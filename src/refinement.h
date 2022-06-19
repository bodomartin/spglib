/* Copyright (C) 2011 Atsushi Togo */
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

#ifndef __refinement_H__
#define __refinement_H__

#include "cell.h"
#include "spacegroup.h"
#include "symmetry.h"

typedef struct {
    Cell *bravais;
    Symmetry *symmetry;
    int *wyckoffs;
    char (*site_symmetry_symbols)[7];
    int *equivalent_atoms;
    int *crystallographic_orbits;
    int *std_mapping_to_primitive;
    double rotation[3][3];
} ExactStructure;

ExactStructure *ref_get_exact_structure_and_symmetry(Spacegroup *spacegroup,
                                                     const Cell *primitive,
                                                     const Cell *cell,
                                                     const int *mapping_table,
                                                     const double symprec);
void ref_free_exact_structure(ExactStructure *exstr);
int ref_find_similar_bravais_lattice(Spacegroup *spacegroup,
                                     const double symprec);
void ref_get_conventional_lattice(double lattice[3][3],
                                  SPGCONST Spacegroup *spacegroup);

#endif
