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

#ifndef SUPPORT_H
#define SUPPORT_H
#include <list>

#include "libSSA/locateOnGenome.h"
#include "Parameters.h"
#include "SupportBreak.h"
#include <gkArrays.h>

using namespace std;
using namespace gkarrays;

class TagInfo;

class Support {
  
  Parameters *parameters;	/* parameters to init constant values*/
  uint *support; 		/* the support by itself */
  uint tag_num; 			/* tag number */
  uint length;			/* support length */
  uint threshold;		/* size of the factors */
  uint tag_length;		/* length of the tag */
  char *current_tag;		/* tag corresponding to this support */
  LocateOnGenome *genome;	/* genome we are matching on */
  gkArrays *tags;		/* index of the tags */
  
/*   uint min_occ;			/\* minimal number of occurrences for the support *\/ */
/*   uint max_occ;			/\* maximal number of occurrences for the support *\/ */
  uint nb_single;		/* number of factors single located */
  uint nb_multiple;		/* number of factors located many times (>1) */
  uint nb_duplicate;		// number of duplicates 
  uint nb_pos_located;		/* number of distinct factors that are located 
				   (at most getLength()) */
  bool almostNormal;          

  // The range support[start_pos_repeat ... end_pos_repeat] is such that
  // each value is >= MIN_OCC_REPETITION
  int start_pos_repeat;	/* start position of the repeat in the support */
  int end_pos_repeat;		/* end position of the repeat in the support */
  

  ChrPosition *pos;		/* Position of an occurrence of one factor 
				   and its strand (if any) */
  uint position_of_location;    /* Position in the read of the k-mer location */
  uint position_of_locations;    /* Position in the read of k-mer locations for 
				    multiple treatment */

  uint *nb_locs;                /* number of occurrences of each factor are 
				   located on the genome */

  uint nb_locs_max;             /* max value among the number of occurrences 
				   computed in nb_locs */

  uint nb_breaks;		/* number of times consecutive factors are not
				   located on the genome */

     
  SupportBreak **breaks;


  pair <uint, uint> *ranges_forward; /* Range of occurrences, in the FM-index, 
					on the forward strand. */
  pair <uint, uint> *ranges_reverse; /* Range of occurrences, in the FM-index, 
					on the reverse strand. */

 public:

  Support(Parameters *p, uint *s, uint tag_num
	  ,LocateOnGenome *g, gkArrays *i);

  Support(const Support &);

  ~Support();

  // REQUESTS //

  /**
   * @pre i < j && i < nb_breaks && j < nb_breaks
   *      && i >= 0 && j >= 0
   * @return true when two breaks can be merged.
   * A break can be merged when we obvisouly have random
   * locations in the middle of a larger break.
   * This can be detected by checking if genome locations
   * of both breaks are in the same region.
   *
   * We also consider that we can merge two breaks
   * when one of them is on a border and (is a chimeric break
   * or does not have a long enough break).
   */
  bool canMergeBreak(uint i, uint j);

  /**
   * @param i: the SupportBreak to return
   * @param consider_strand: should we consider the strand when returning
   *                the break number i?
   * @pre 0 <= i < nb_breaks
   * @return the break number i in the read iff 
   *          - we don't consider the strand,
   *          - it is on forward strand,
   *          - it is unknown
   *         the break number getNbBreaks() - i - 1, otherwise.
   */
  SupportBreak *getBreak(uint i, bool consider_strand=false);

  /**
   * @return end_pos_repeat
   */
  int getEndPosRepeat();

  /**
   * @return tags 
   */
  gkArrays *getIndexTags();

  /**
   * @return genome
   */
  LocateOnGenome *getGenome();
  
  /**
   * @return the support length
   */
  uint getLength();

  /**
   * @return the 'best' location of a factor (of length getThreshold()) in the
   *         read. 
   *         This is the location in the middle of a "run" of 
   *         single-located factors. If no such location is found
   *         take a location in the middle of a run of located factors.
   *         If nothing is found return NULL
   */
  ChrPosition *getLocation();

  /**
   * @return all locations of a factor (of length getThreshold()) in
   *         the read and its numbers.  This is the first factor with
   *         the most of occurrences on the reference.  If too much
   *         locations return only getGenome()->getNbLocations()
   *         locations If nothing is found return a pair<ChrPosition
   *         **, uint> empty.  In addition, position_of_locations is
   *         set after pair<ChrPosition **,uint> which determines
   *         position of the k-mer chosen for locations
   */
  pair<ChrPosition **, uint> getLocations(uint &position_of_locations);

  /**
   * @return min_occ
   * @post min_occ > 0 iff getNbSingle()+getNbMultiple() == getLength()
   */
  /* uint getMinOcc();	        */

  /**
   * @return max_occ
   * @post max_occ == 0 iff isNone()
   */
/*   uint getMaxOcc(); */

  /**
   * @return nb_breaks
   */
  uint getNbBreaks();

  /**
   * @return nb_duplicate
   * @post getNbDuplicate() <= getNbMultiple()
   */
  uint getNbDuplicate();

  /**
   * @return nb_locs
   */
  uint *getNbLocs();

  /**
   * @return nb_locs at position i
   */
  uint getNbLocs(uint i);

  /**
   * @return the max value in nb_locs
   */
  uint getNbLocsMax();

  /**
   * @return nb_multiple
   */
  uint getNbMultiple();

  /**
   * @return the number of distinct positions where a factor of length
   *         getThreshold() is located.
   */
  uint getNbPositionsLocated();

  /**
   * @return nb_single
   */
  uint getNbSingle();

  /**
   * @return parameters
   */
  Parameters *getParameters();

  /**
   * getLocation() is the position of the k-mer at position
   * getPositionOfLocation() in the read.
   * Return getTagLength() if getLocation() == NULL
   */
  uint getPositionOfLocation();

  /** Return the range of the k-mer at position i of
   * the profile on a specific strand
   *
   * @param i the number of the kmer if the support
   * @param strand the strand on wich we get the range
   * @return the range of the kmer i in the genome index
   */
  pair<uint, uint> getKmerRange(uint i, uint strand);
  
  /**
   * @return the length of the repetitioon or 0 if not
   */
  int getRepeatLength();

  /**
   * @return start_pos_repeat, -1 if no repeat exists
   */
  int getStartPosRepeat();

  /**
   * @return support
   */
  uint *getSupport();

  /**
   * @return support at position i
   */
  uint getSupport(uint i);

  /**
   * @return the tag for this support
   */
  char *getTag();

  /**
   * @return the tag number
   */
  uint getTagNum();

  /**
   * @return the length of the current read
   */
  uint getTagLength();

  /**
   * @return threshold
   */
  uint getThreshold();

  /**
   * A repetition is built by a number of consecutive k-mers with 
   * a number of occurrence >= parameters->min_occ_repetition
   * @return (getRepeatLength() >= percent_min_unique_repetition*getNbPositionsLocated())
   * 
   */
  bool hasRepeat();

  /**
   * @return isAlmostNormal
   */
  bool isAlmostNormal();

  /**
   * return getLength() == getNbMultiple()+getNbSingle()
   */
  bool isContinuous();

  /**
   * @return getNbDuplicate() >= params->percent_min_duplicate*getNbPositionsLocated()
   */
  bool isDuplicate();

  /**
   * @return ! isDuplicate 
   * && (getNbSingle() < parameters->percent_max_unique_multiple*getNbPositionsLocated()
   */
  bool isMultiple();
  
  /**
   * @return (nb_single == 0 && nb_multiple == 0)
   * ||  (nb_breaks == 1 && breaks[0]->getPositionStartBreak() == 0
   *	&& breaks[0]->getPositionEndBreak() == getLength()-1)
   */
  bool isNone();
  
  /**
   * @return !isNone() && !isMultiple() && !isDuplicate()
   */ 
  bool isSingle();
 private:

  /**
   * Compute the best location either for k-mer only single located
   * or for any k-mer.
   * @param only_single if the parameter is true we only search for k-mers that 
   *        are located unambiguously (ie. only once)
   * @post getLocation() will return the location (if a suitable one was found)
   * NULL otherwise
   */
  void computeBestLocation(bool only_single);
  
  /**
   * Will merge breaks using the following algorithm:
   * If a break is at the beginning of the read and the next one
   * ! isNiceBreak()
   *   -> merge them (if there are less than MIN_BASES_BEFORE_BREAK nt between them) 
   * If a break is at the end of the read and the previous one
   * ! isNiceBreak()
   *   -> merge them (if there are less than MIN_BASES_BEFORE_BREAK nt between them) 
   * For each short chimeric break:
   *   if isNiceMerge() is true with the next or previous break and is not too far
   *     -> merge
   *   else if the next or previous break is chimeric
   *      merge with the closest one if it is not too far
   *      (< 2 * getThreshold)
   *   else if the next or previous break is at the start or at the 
   *        end of the read
   *     -> merge
   * For each chimeric break (not yet merged):
   *   if isNiceMerge() is true with the next or previous break and is not too far
   *     or if the next break is chimeric and not too far (< 2 * getThreshold())
   *        from the current break
   *     or if the next or previous break is at the start or at the 
   *        end of the read
   *     -> merge
   * For each break:
   *    Merge with the next one iff the resulting break isVeryNiceBreak(),
   *    the current break or the next one are not isVeryNiceBreak()
   *    and if the resulting break lengthis <= getThreshold(). 
   *    We also merge in the case two consecutive breaks are overlapping
   */
  void tryToMergeBreaks();

  /**
   * Merge break i vith break i+1.
   * The other arrays are updated accordingly.
   * @param no_chim: is it a chimeric break being merged?
   * @post merged_with[i+1] == i
   *       breaks[i+1] == breaks[i]
   *       chimeric_breaks[i+1] == chimeric_breaks[i]
   *       short_chimeric_breaks[i+1] == short_chimeric_breaks[i]
   */
  void merge(list<SupportBreak*> &break_list, list<bool> &chimeric_breaks,
             list<bool> &short_chimeric_breaks,
             list<SupportBreak*>::iterator &i, list<bool>::iterator &chim, 
             list<bool>::iterator &shortB, bool no_chim=false);
  
};

#endif

