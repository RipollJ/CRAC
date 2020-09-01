/******************************************************************************
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
*  with additional contribution for the packaging of:	                      *
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
******************************************************************************/

#ifndef LOCATE_ON_GENOME
#define LOCATE_ON_GENOME

extern "C" {
#include <sys/time.h>
}

#include <ctime>
#include <fstream>
#include <iostream>
#include <utility>
#include "SSA.h"
#include "chrPosition.h"
#include "GenomeInfo.h"
using namespace std;


class LocateOnGenome {
 private:
  GenomeInfo *genomeInfo;
  ostream *out_uniq, *out_multiple, *out_none, *out_sam;
  TFMindex *index;
  ulong nb_tags_treated;
  ulong nb_displayed;
  ulong nbOccs_lastTag;
  time_t chrono_sec;
  suseconds_t chrono_usec;
  bool display_nb;

  void display_occ(ostream *stream, ulong *occ, ulong n, ulong tag_length, bool sens);


 public:
  /**
   * @param file_base: Basename of the index file (ie. the common prefix
   *                   to the .ssa index and to the .conf file).
   */ 
  LocateOnGenome(uchar *file_base);

  ~LocateOnGenome();

  void displaySingleLocation(uchar *tag, ulong tag_length, ulong n_sens, ulong n_antisens, ulong *occ_sens, ulong *occ_antisens);
  void displayMultipleLocations(uchar *tag, ulong tag_length, ulong n_sens, ulong n_antisens, ulong *occ_sens, ulong *occ_antisens);
  void displayNoLocation(uchar *tag, ulong tag_length);

  float getElapsedTime();


  void getLocations(uchar *tag, ulong tag_length, ulong &nb_sens, ulong &nb_antisens, ulong **poccs_sens, ulong **poccs_antisens);
  
  /**
   * @return all locations of a kmer at position_of_locations and its
   * number of occurrences.  In addition, position_of_locations is set
   * after pair<ChrPosition **,uint> which determines position of the
   * k-mer chosen for locations
   */
  pair<ChrPosition **,uint>getLocations(uchar *tag, uint klength, uint read_length, uint &position_of_locations);
  
 
  uint getPositionOfLocations();

  uint getNbLocations();

  /* 
   * @param chrPos : a ChrPosition object
   * @return the absolute position of chrPos
   */
  ulong getAbsolutePosition(ChrPosition *chrPos);

  /**
   * @param pos: absolute position
   * @return the chromosome which corresponds to that position
   */
  const uchar *getChromosome(ulong position);

  /**
   * @param chr_name: name of chr
   * @return an id of the chromosome which corresponds to that name.
   */
  ulong getIdChromosomeFromName(char *chr_name);

  /**
   * @param pos: absolute position
   * @return an id of the chromosome which corresponds to that position.
   * This id is mainly useful when one wants to compare two positions
   * and wants to know if they're on the same chromsome.
   */
  ulong getIdChromosome(ulong position);

  /**
   * @param pos: absolute position
   * @return the relative position on the genome.
   */
  ChrPosition *getChrPos(ulong pos, int strand=0);

  pair<uint, uint> getFMIndexRange(uchar *tag, ulong tag_length);

  pair<uint, uint> getFMIndexReverseRange(uchar *tag, ulong tag_length);

  void getOccurrencesFMIndexRange(ulong sp, ulong ep, ulong **occ, ulong m=0, ulong nb_disp=0);
  
  /**
   * @return A pointer to the GenomeInfo object for hvaing all the information
   *         on sequence length, name, IDs, ...
   */
  GenomeInfo * getGenomeInfo();

  uchar *getGenomeSubstring(ulong pos, ulong length, int strand=1,
                            char *chr_name = NULL);


  ulong getNbTagsTreated();
  void getNbOccurrences(uchar *tag, ulong tag_length, ulong &nb_sens, ulong &nb_antisens);
  ulong getNbOccsLastLocatedTag();

  void locateTags(fstream &file, ulong tag_length, ulong nb_threads = 1);

  void locateTag(uchar *tag, ulong tag_length);


  /**
   * @param disp: true iff read number precedes each entry in the output files
   */
  void setDisplayNumbers(bool disp);

  void setOutputStreams (ostream *uniq, ostream *multiple, ostream *none);
  void setOutputUniq(ostream *s);
  void setOutputMultiple(ostream *s);
  void setOutputNone(ostream *s);
  void setOutputSAM(ostream *s);

  void setNbLocations(ulong nb_displayed);

  bool startChrono();
  
};

#endif
