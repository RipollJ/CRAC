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

#ifndef SUPPORT_BREAK_H
#define SUPPORT_BREAK_H
#include <utility>
#include "Parameters.h"
#include "CandidatBreak.h"
#include "libSSA/chrPosition.h"
#include "libSSA/locateOnGenome.h"
#include "types.h"

class Support;

class SupportBreak {

  CandidatBreak **candidats;
  ulong nb_candidats;
  CandidatBreak *candidat_chosen;

  uint k;                       /* Threshold, size of k-mer used */
  /* bool duplicated;             /\* true iff there are two possible causes *\/ */

  bool repeated;		/* true iff the break is in a repeated region.
  				   In that case, loc_{start,end}_break are not
  				   defined. */
  
  uint pos_start_break;		/* Position in the support of the start 
				   of the break. */
  uint pos_end_break;		/* Position in the support of the end 
				   of the break.*/

  uint override_pos_end_break; /* Position in read of the end 
				   of the break. (ie. pos_end_break + k - 1) */

  float score_outside_break;	/* the score (can be median, average
				   or whatever we can think about)
				   of the support, just before or
				   after the break. */
  
  float score_during_break;	/* Well, that's self-speaking, isn't it? */

  float average_low_inside;     /* Average of the lowest values inside the break */
  float average_high_inside;    /* Average of the highest values inside the break */
  float pc_ones_left_inside,    /* Percentage of 1s in the support on the */
    pc_ones_left_outside,       /* border of the breaks */
    pc_ones_right_inside,
    pc_ones_right_outside;
  float score_inside_averages;  /* Score between the two averages */

  Support *support;		/* The support we are working on. */
  
  Parameters* parameters;	/* Parameters chosen by the user or chosen 
				   by default by the great developers. */

  pair<uint, uint> *range_occs_fwd; /* Range of the occurrences on the 
				       forward strand. */

  pair<uint, uint> *range_occs_rev; /* Range of the occurrences on the 
				       reverse strand. */

  uint loc_end_break;   /* the first absolute location after the 
                            break */

  uint loc_start_break; /* the last absolute location before the 
                            break */
  
  uint nb_loc_end;       	/* the number of location at the current pos after the break */
  uint nb_loc_start; 		/* the number of location at the current pos before the break */

  int strand_end_break;         /* the strand  of the end_break position*/

  int strand_start_break;       /* the strand of the start_break position */

  uint inside_quartile1;        /* first quartile of the support
                                 * inside the break*/
  uint inside_quartile4;        /* last quartile of the support
                                 * inside the break*/
  uint nb_merges;       /* the number of merges that have been done to create this break */

 public:

  SupportBreak(uint start, uint end, 
	       Support *s, Parameters *p,
	       pair<uint, uint> *range_occs_forward,
	       pair<uint, uint> *range_occs_reverse,
         uint nb_merges = 0);

  /**
   * Merge the two breaks in parameters.
   */
  SupportBreak(SupportBreak &b1, SupportBreak &b2);

  ~SupportBreak();

  /**
   * Computes score inside and outside break. It depends of
   * pos_start_break and pos_end_break
   */
  void computeScore();

  /**
   * @param pos_start: a position to define a start break 
   * @param pos_end: a position to define an end break   
   * Computes score inside and outside break between a pos_start and a
   * pos_end defined in argument.
   */
  void computeScore(uint pos_start, uint pos_end);

  /**
   * @return the average low support inside the break
   */
  float getAverageLowInside();

  /**
   * @return the average high support inside the break
   */
  float getAverageHighInside();

  /**
   * @return the candidat chosen for the break
   */
  CandidatBreak *getCandidatChosen();

  /**
   * @return the candidat number i among the potential candidats
   */
  CandidatBreak *getPotentialCandidat(uint i);

  /**
   * @param: position of the location
   * @return the chromosome corresponding to the position
   * in the genome of the factor before
   * or after the break.
   * ~0 if p == START_BREAK && getPositionStartBreak() == 0
   * || p == END_BREAK && getPositionEndBreak() == support->getLength()
   */
  const uchar *getChr(pos_location_break);


  /**
   * @param: position of the location
   * @return the ID of the chromosome corresponding to the position
   * in the genome of the factor before
   * or after the break.
   * ~0 if p == START_BREAK && getPositionStartBreak() == 0
   * || p == END_BREAK && getPositionEndBreak() == support->getLength()
   */
  ulong getChrId(pos_location_break);

  /**
   * @param: position of the location
   * @return a ChrPosition object corresponding to the position
   * in the genome of the factor before
   * or after the break.
   * NULL if p == START_BREAK && getPositionStartBreak() == 0
   * || p == END_BREAK && getPositionEndBreak() == support->getLength()
   */
  ChrPosition *getChrPosition(pos_location_break p);


  /**
   * @return genome
   */
  LocateOnGenome *getGenome();

  /**
   * @return The gap length between the location of the factors on the
   * genome, before the break and after the break. If no factor is
   * located before or after the break return 0.
   */
  gap_size_t getGenomeGapLength();

  /**
   * @return the value delimiting the first and second quartile of the
   * support inside the break (ie. 25% of the values are lower or
   * equal than getInsideQuartile1())
   */
  uint getInsideQuartile1();

  /**
   * @return the value delimiting the third and fourth quartile of the 
   * support inside the break (ie. 25% of the values are greater or
   * equal than getInsideQuartile1())
   */
  uint getInsideQuartile4();

  /**
   * @return the score inside the break i.e. the average of support
   * of each pos inside the break.
   */
  float getInsideScore();
    
  /**
   * @return the location of the factor just after the break.
    */
  ChrPosition *getLocationEndBreak();

  /**
   * @return the location of the factor just before the break.
   */
  ChrPosition *getLocationStartBreak();


  /**
   * @return length of the support as seen using the current threshold getThreshold().
   * getLength() + support->getThreshold() - getThreshold() == support->getLength() 
   */
  uint getLength();
  
  /**
   * @return the number of candidats saved during the extension process
   */
  ulong getNbCandidats();

  /**
   * @return the number of merge processed to create this break
   */
  uint getNbMerges();

  /**
   * Return the number of indels on the genome.
   * @return getGenomeGapLength() - support->getTagBreakLength()
   *                                         , if getGenomeGapLength() > 0
   *    INT_MAX, otherwise.
   */
  int getNbGenomeIndels();

  /**
   * Return the number of indels in the tag.
   * @return getTagBreakLength() - support->getGenomeGapLength()
   */
  int getNbTagIndels();

  /**
   * @return the score outiside the break i.e.
   * the average of support for each position outside the break
   */
  float getOutsideScore();

  /**
   * @param strand: if -1 consider the break is on the other strand
   * @return pos_end_break (ie. position where the k-mer begins after the break)
   *         or getSupportLength() - 1 - getPositionStartBreak() (if strand == -1)
   */
  uint getPositionEndBreak(int strand=1);

  /**
   * @return pos_end_break + k -1 or override_pos_end_break (ie. the
   * real pos_end_break in more detailed search)
   */
  uint getPositionEndBreakOnTheRead();

  /**
   * @param: position of the location
   * @return absolute posittion in the genome
   * corresponding the factor before
   * or after the break.
   * ~0 if p == START_BREAK && getPositionStartBreak() == 0
   * || p == END_BREAK && getPositionEndBreak() == support->getLength()
   */
  ulong getPositionLocation(pos_location_break);

  /**
   * @param strand: the strand we must consider to return the position
   * @return pos_start_break or getSupportLength() - 1 - getPositionEndBreak()
   *         if strand == -1
   */
  uint getPositionStartBreak(int strand=1);

  /**
   * @param i: i == FORWARD_STRAND || i == REVERSE_STRAND
   * @return range_occurrences[END_BREAK][i]
   */
  pair<uint, uint> getRangeEndBreak(uint i);

  /**
   * @param i: i == FORWARD_STRAND || i == REVERSE_STRAND
   * @return range_occurrences[START_BREAK][i]
   */
  pair<uint, uint> getRangeStartBreak(uint i);

  /**
   * Compute the score to evaluate the coverage between two exons from
   * getScoreInsideBreak() and getScoreOutsideBreak().  @return a
   * float between -1 and 1 (-1: biological reason, 0: not different, 1: error)
   */
  float getScoreComputedInterExon();

  /**
   * Compute the score to evaluate the coverage inside an exon from
   * getScoreInsideBreak() and getScoreOutsideBreak().  @return a
   * float between -1 and 1 (-1: biological reason, 0: not different, 1: error)
   */
  float getScoreComputedIntraExon();

  /**
   * @return the score computed between the average of the 50% highest values and
   * the average of the 50% lowest values, inside the break
   */
  float getScoreInsideAverages();

  /**
   * @return score_inside_break
   */ 
  float getScoreInsideBreak();

  /**
   * @return score_outside_break
   */
  float getScoreOutsideBreak();

  /**
   * @param: position of the location
   * @return strand in the genome
   * corresponding the factor before
   * or after the break.
   * 0 if p == START_BREAK && getPositionStartBreak() == 0
   * || p == END_BREAK && getPositionEndBreak() == support->getLength()
   */
  int getStrandLocation(pos_location_break);

  /**
   * The break length between the the break and after the break, in the tag.
   * @return getPositionEndBreak()-getPositionStartBreak() + 1
   */
  uint getTagBreakLength();

  /**
   * @return the length of the current read
   */
  uint getTagLength();

  /**
   * @param : the score outside the break 
   * @return a score computed by a
   * function from the score to define the threshold between errors
   * and biological reasons like chimera or splice.
   */
  float getThresholdForScoreInterExon(float score);

  /**
   * @return the threshold (ie. the k-mer length)
   */
  uint getThreshold();
    
    /**
   * @param : the score outside the break
   * @return a score computed by a function from the score to
   * define the threshold between errors and biological reasons
   * like snp or short indel. 
   */
  float getThresholdForScoreIntraExon(float score);

  /**
   * Return true if the support has no cover
   * @return score_outside_break <= 1 + EPSILON
   *  && score_inside_break <= 1 + EPSILON
   */
  bool hasNoCover();

  /**
   *  @return true if the getTagBreakLength() is long enough
   *  or if there is no start/end break (in that case we don't know
   *  where the break starts/finishes).
   *  There is no good reason for it to be less than 
   *  getThreshold() - parameters->max_bases_randomly_matched.
   */
  bool hasLongEnoughTagBreak();
  
  /**
   *  @return true if getPositionEndBreak() + 1 == support->getLength()
   *             or if hasNoStartBreak() (when strand == 1)
   */
  bool hasNoEndBreak(int strand=1);
  
  /**
   * @return true if getPositionStartBreak() == 0 (when strand == 1)
   *             or if hasNoEndBreak() (when strand == -1)
   */
  bool hasNoStartBreak(int strand=1);

  /**
   * @return true iff we have set the end position of the break on the read 
   *         using setPositionEndBreakOnTheRead()
   */
  bool hasPositionEndBreakOnTheRead(int strand=1);
  
  /**
   * @return return ! hasNoCover() 
   *                && (score_during_break >= score_outside_break 
   *                || isBiologicalIntraEvent()
   *                || isBiologicalInterEvent())
   */
  bool isBiologicalEvent();
  
  /**
   * @return getScoreComputedIntraExon() 
   *         <= parameters->p_value_variation_biological
   *        && !isTagIndel() && !isTagSubstitution()	 
   */
  bool isBiologicalIntraEvent();

  /**
   * @return getScoreComputedInterExon() <= 
   *         parameters->p_value_variation_biological	 
   *        && !isChimera() 
   */
  bool isBiologicalInterEvent();

  /**
   * For a chimera the break length in the tag
   * cannot be greater than threshold - 1.
   * There is a chimera is two parts of the read
   * are not located on the same strand or not on 
   * the same chromosome or when the break on the 
   * genome is larger than parameters->max_splice_length
   * or when the break length is long enough
   * @return true if the break is a chimera
   */
  bool isChimera();

  /**
   * @return true if a part of the support is low because of an error
   *         occuring at the same time as a mutation
   */
  bool isDeviated();

  /**
   * @return true if the cause is duplicated, i.e. there are two
   *         possible splice or SNV
   */
  bool isDuplicated();

  /**
   * A genome deletion must be larger than
   * parameters->max_bio_ins_del bases
   * (otherwise we consider it to be a tag deletion).
   * A genome deletion is characterised by a larger break
   * in the tag than in the genome.
   * @return true if the break is a genome deletion
   */
  bool isGenomeDeletion();

  /**
   * For a genome insertion the break length in the tag
   * cannot be greater than threshold - 1.
   * Moreover, the insertion length must be greater than
   * parameters->max_bio_ins_del (otherwise it is
   * considered as an insertion in the tag).
   * @return true if the break is a genome insertion
   */
  bool isGenomeInsertion();

  /**
   * @return true iff 
   *  - locations are defined both before and after the break
   *  - both locations are on the same strand and on the same chromosome
   *  - genome gap length is not larger than parameters->max_splice_length
   *    and is coherent with the strand.
   */
  bool isNiceBreak();

  /**
   * @return true iff merging *this and *sb will produce a break
   *         where isNiceBreak() == true.
   */
  bool isNiceMerge(SupportBreak *sb);

  /**
   * @return true iff the break is in a repeated region and we didn't find 
   * a ``nice'' location (ie. locations before and after the break must be close).
   */
  bool isRepeated();

  /**
   * @return true if the support falls on the left or on the right border of
   * the break between the inside and outside.
   */
  bool isSupportFalling();

  /**
   * @return true if the support falls on the left border of the break
   * between the inside and outside (when strand == 1).
   * Return isSupportFallingRight() when strand == -1
   */
  bool isSupportFallingLeft(int strand=1);

  /**
   * @return true if the support falls on the right border of the break
   * between the inside and outside (when strand == 1)
   * Return isSupportFallingLeft() when strand == -1
   */
  bool isSupportFallingRight(int strand = 1);

  /**
   * The TagBreakLength must be equal to the GenomeGapLength. 
   * We don't compare the tag break length to the threshold since
   * we may have several substitutions in a raw leading to a larger gap.
   * @return true if the break is a substitution
   */
  bool isTagSubstitution();

  /**
   * An indel in a tag corresponds to an indel (length != 0)
   * whose maximal length is parameters->max_bio_ins_del.
   * In case of tag deletion the break length of tag must be less
   * than threshold.
   * With tag indel, the gap on the genome can't be negative. It would be
   * a chimera.
   * @return true if the break is a short indel
   */
  bool isTagIndel();

  /**
   * @return true iff 
   *  - locations are defined both before and after the break
   *  - both locations are on the same strand and on the same chromosome
   *  - genome gap length is not larger than parameters->max_intron_length
   *    and is coherent with the strand.
   */
  bool isVeryNiceBreak();

  /**
   * @return true iff merging *this and *sb will produce a break
   *         where isVeryNiceBreak() == true.
   */
  bool isVeryNiceMerge(SupportBreak *sb);


  /**
   * Compute the inside score between pos_start_break and
   * pos_end_break. An average support values inside the break is
   * computed (only pos with a getNbLocs() == 0 are considered).
   */
  void setInsideScore();


  /**
   * @param pos_start: a position to define a start break 
   * @param pos_end: a position to define an end break 
   * Compute the inside score between pos_start and pos_end. An average
   * support values inside the break is computed (only pos with a
   * getNbLocs() == 0 are considered because it may bias the
   * computation).
   */
  void setInsideScore(uint pos_start, uint pos_end);

  /**
   * loc_end_break = loc_end.
   */
  void setLocationEndBreak(uint loc_end);

  /**
   * loc_start_break = loc_start.
   */
  void setLocationStartBreak(uint loc_start);

  /**
   * pos_end_break = pos_end
   */
  void setPositionEndBreak(uint pos_end);

  /**
   * override_pos_end_break = pos_end
   */
  void setPositionEndBreakOnTheRead(uint pos_end);

  /**
   * pos_start_break = pos_start
   */   
  void setPositionStartBreak(uint pos_start);

  /**
   * @post getThreshold() == k
   */
  void setThreshold(uint k);

  /**
   * Compute the outside score between the pos_start_break and the
   * pos_end_break. An average of (at most)
   * parameters->support_score_window_length support values outside
   * the break is computed.
   */
  void setOutsideScore();

  /**
   * @param pos_start: a position to define a start break 
   * @param pos_end: a position to define an end break 
   * Compute the outside score between pos_start and pos_end. An
   * average of (at most) parameters->support_score_window_length
   * support values outside the break is computed.
   */
  void setOutsideScore(uint pos_start, uint pos_end);

  /**
   * Sets the forward range.
   */
  void setRangesForward(pair<uint, uint> *);
  
  /**
   * strand_end_break = strand_end.
   */
  void setStrandEndBreak(int strand_end);
  
  /**
   * strand_start_break = strand_start.
   */
  void setStrandStartBreak(int strand_start);

  /**
   * Sets the reverse range.
   */
  void setRangesReverse(pair<uint, uint> *);

  /**
   * Sets the support.
   */
  void setSupport(Support *);

 private:

  /**
   * Adjust the boundaries of getCandidatChosen because 
   * it has not necessary the good ones after the extension process
   */
  void adjustBoundaries();


  /**
   * Return true if we can check locations at positions start and end
   */
  bool canCheckPositions(ulong start, ulong end);


  /**
   * Merge the two breaks and check if the result is nice (or very nice if
   * the boolean is set at true)
   */
  bool check_nice_merge(SupportBreak *sb1, SupportBreak *sb2, 
                        bool very_nice);


  /**
   * @param pos_start: a position to define a start break 
   * @param pos_end: a position to define an end break
   * Sort the inside support between pos_start and pos_end. Compute
   * the average on the low and on the high values and compute the
   * score between those two averages.
   */
  void computeHighAndLowAverages(uint pos_start, uint pos_end);


  /**
   * Method to choose the best location.
   */
  void fillWithClosestMatch();


  /**
   * @param a table of candidats
   * @return the position of best candidat among the table of potential candidats during a step of the extension
   * the best candidat is :
   * - nice_candidat (at least a splice)
   * - single consistent
   * - the minimal gap on the genome
   * And notify a duplicated candidat if the table contains several good candidats
   */
  ulong getBestCandidat(CandidatBreak** potential_candidats, ulong nb_potential_candidats, bool after_extension = false);
  
  /**
   * @return getBestCandidat(this->candidats,this->nb_candidats,true) after the extension process
   */
  ulong getBestCandidat();
  
  /**
   * Initialises the members (called by the different constructors);
   */
  void init();

  /**
   * @return (potential_candidats[pos_i] eq potential_candidats[pos_j]) ? true : false
   */
  bool isSameCandidat(CandidatBreak **potential_candidats, ulong pos_i, ulong pos_j);

  /**
   * @return (potential_candidats[pos_i].end eq potential_candidats[pos_j].end) ? true : false
   */
  bool isSameCandidatEnd(CandidatBreak **potential_candidats, ulong pos_i, ulong pos_j);

  /**
   * @return (potential_candidats[pos_i].start eq potential_candidats[pos_j].start) ? true : false
   */
  bool isSameCandidatStart(CandidatBreak **potential_candidats, ulong pos_i, ulong pos_j);

  /**
   * Save the best candadidat for each extension step of the extension
   */
  void saveCandidates(uint pos_before_break, uint pos_after_break,
		      ulong *pos_occ_start, uint nb_loc_start, 
		      int strand_start, ulong *pos_occ_end, 
		      uint nb_loc_end, int strand_end, 
		      CandidatBreak ***potential_candidats, 
		      ulong &nb_potential, ulong &nb_potential_max);

  /**
   * Method to select potential candidats during the extension process.
   */
  void updateLocation(uint pos_before_break,
		      uint pos_after_break);

  
};

ostream &operator<<(ostream &, SupportBreak &);
#endif
