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

#include <stdlib.h>
#include <utility>
#include "Support.h"
#include "utils.h"

Support::Support(Parameters *p, uint *s, uint num
		,LocateOnGenome *g, gkArrays *ind):
  parameters(p), support(s), tag_num(num), length(ind->getSupportLength(num)),
  threshold(ind->getThreshold()), tag_length(ind->getTagLength(num)),
  current_tag(ind->getTag(num)), genome(g), tags(ind), nb_pos_located(0) {


  //   cerr << "tag length: "<< tag_length << endl;
  // cerr << "current_tag: " << current_tag << endl << endl;
  ulong nb_occs_sens, nb_occs_antisens;
//   ulong *occs_sens, *occs_antisens;
  ulong total_occ;
  bool located = true;

  ulong max_nb_breaks = 10;
  breaks = (SupportBreak**)malloc(sizeof(SupportBreak*)*max_nb_breaks);
  ulong pos_break=0;
  

//   pos_occ = new uint[tag_length-threshold+1];
//   strand_occ = new int[tag_length-threshold+1];


  nb_locs = new uint[length];

  start_pos_repeat = -1;
  end_pos_repeat = -1;
  nb_breaks = 0;
  nb_single = 0;
  nb_multiple = 0;
  nb_duplicate = 0;
  almostNormal = false;
  nb_locs_max = 0;

  ranges_forward = new pair<uint, uint>[length];
  ranges_reverse = new pair<uint, uint>[length];
  
  for (ulong j=0; j < length; j++) {
    ranges_forward[j] = genome->getFMIndexRange((uchar *)&current_tag[j],threshold);
    nb_occs_sens = (ranges_forward[j].first <= ranges_forward[j].second) ? 
      (ranges_forward[j].second - ranges_forward[j].first + 1):0;

    ranges_reverse[j] = genome->getFMIndexReverseRange((uchar *)&current_tag[j],threshold);
    nb_occs_antisens = (ranges_reverse[j].first <= ranges_reverse[j].second) ? 
      (ranges_reverse[j].second - ranges_reverse[j].first + 1):0;
    nb_locs[j] = nb_occs_antisens + nb_occs_sens;
    if (nb_locs_max < nb_locs[j])
      nb_locs_max = nb_locs[j];
  }
  
  computeBestLocation(true);

  // temporary repeat start, end and length
  int minRepeat = 0;
  int start_pos_r = -1;
  int end_pos_r = -1;
  for (ulong j=0; j < length; j++) {    
    total_occ = nb_locs[j];
    
    // almost normal
    if (!almostNormal && (j+1) < length) {
      almostNormal =  ((support[j] > support[j+1] 
			  && support[j]*(1.0 - parameters->percent_support_variation_almost_normal) > support[j+1]) 
			 || (support[j] < support[j+1] && support[j]*1/(1.0 - parameters->percent_support_variation_almost_normal) < support[j+1]));    
    }    
    if (total_occ > 0) {
      nb_pos_located++;
      // end of repetition in case of repetition (we suppose one
      // repetition per tag)
      if (total_occ < parameters->min_occ_repetition 
	  && start_pos_repeat != -1 && end_pos_repeat == -1 ){
	end_pos_repeat = j-1;
      }
      // only the longest repetition is saved
      if (end_pos_repeat == (int)(j-1)) {
	if (getRepeatLength() > minRepeat){
	  minRepeat = getRepeatLength();
	  start_pos_r = start_pos_repeat;
	  end_pos_r = end_pos_repeat;
	}
	start_pos_repeat = -1;
	end_pos_repeat = -1;
      }
      
      // duplication
      if (total_occ <= parameters->max_localisation_duplication) {
	
	if (total_occ >= parameters->min_localisation_duplication
	    && total_occ <= parameters->max_localisation_duplication) {
	  nb_duplicate++;
	}
      } else {
     	nb_multiple++;
	if(start_pos_repeat == -1 
	   && total_occ >= parameters->min_occ_repetition){
	  start_pos_repeat = j;
	}
      }
          
      if (total_occ == 1){
	nb_single++;
      }
      
      // We didn't attribute the score of the previous break 
      // and we're yet quite far from the break, so we attribute it right now!
      if (! located) {
	// We have more breaks than expected.
	// We expand the array.
	if (nb_breaks > max_nb_breaks) {
	  max_nb_breaks *= 2;
	  breaks = (SupportBreak**)realloc(breaks, sizeof(SupportBreak*)*max_nb_breaks);
	}

	  breaks[nb_breaks-1] = new SupportBreak(pos_break, j-1, 
						 this, parameters,ranges_forward,
						 ranges_reverse);
      }

      located = true;

    }else{
      if (located) {
	nb_breaks++;
	pos_break = j;
      }
      located = false;
    } 
  }
  
  // save the biggest repetition
  if (getRepeatLength() <= minRepeat){
    start_pos_repeat = start_pos_r;
    end_pos_repeat = end_pos_r;
  }
    
  if (! located) {
    if (nb_pos_located) {
      // We have more breaks than expected.
      // We expand the array.
      if (nb_breaks > max_nb_breaks) {
	max_nb_breaks *= 2;
	breaks = (SupportBreak**)realloc(breaks, sizeof(SupportBreak*)*max_nb_breaks);
      }
      breaks[nb_breaks-1] = new SupportBreak(pos_break, length-1,
					     this, parameters,ranges_forward,
					     ranges_reverse);
    } else {
      nb_breaks--;
    }
  } 

  // Now we try to merge breaks
  if (nb_breaks > 1) {
    tryToMergeBreaks();
  }

  if (getLocation() == NULL) {
    computeBestLocation(false);
  }

  almostNormal = almostNormal && ! hasRepeat() && isContinuous();
}

Support::Support(const Support &s): parameters(s.parameters), tag_num(s.tag_num),
                                    length(s.length), threshold(s.threshold), 
                                    tag_length(s.tag_length), genome(s.genome),
                                    tags(s.tags), nb_single(s.nb_single),
                                    nb_multiple(s.nb_multiple), nb_duplicate(s.nb_duplicate),
                                    nb_pos_located(s.nb_pos_located),
                                    almostNormal(s.almostNormal),
                                    start_pos_repeat(s.start_pos_repeat),
                                    end_pos_repeat(s.end_pos_repeat), nb_breaks(s.nb_breaks)
                           
{
  support = new uint[length];
  nb_locs = new uint[length];
  ranges_reverse = new pair<uint,uint>[length];
  ranges_forward = new pair<uint,uint>[length];
  for (uint i = 0; i < length; i++) {
    support[i] = s.support[i];
    nb_locs[i] = s.nb_locs[i];
    ranges_forward[i] = s.ranges_forward[i];
    ranges_reverse[i] = s.ranges_reverse[i];
  }

  current_tag = new char[tag_length+1];
  for (uint i = 0; i <= tag_length; i++) {
    current_tag[i] = s.current_tag[i];
  }

  pos = new ChrPosition(*(s.pos));
  breaks = new SupportBreak*[nb_breaks];
  for (uint i=0; i < nb_breaks; i++) {
    breaks[i] = new SupportBreak(*(s.breaks[i])); 
    breaks[i]->setSupport(this);
    breaks[i]->setRangesForward(ranges_forward);
    breaks[i]->setRangesReverse(ranges_reverse);
  }
  
}

Support::~Support() {
  for (uint i=0; i < nb_breaks; i++) {
    delete breaks[i];
  }
  delete [] ranges_forward;
  delete [] ranges_reverse;
//   delete [] pos_occ;
//   delete [] strand_occ;
  delete [] current_tag;
  if (pos)
    delete pos;
  delete [] nb_locs;
  delete [] support;
  free(breaks);
}

SupportBreak *Support::getBreak(uint i, bool consider_strand) {
  if (! consider_strand || ! getLocation() 
      || getLocation()->getStrand() == 1)
    return breaks[i];
  else
    return breaks[getNbBreaks() - i - 1];
}

int Support::getEndPosRepeat() {
  return end_pos_repeat;
}

gkArrays *Support::getIndexTags(){
  return tags;
}

LocateOnGenome *Support::getGenome() {
  return genome;
}

uint Support::getLength() {
  return length;
}

ChrPosition *Support::getLocation(){
  return pos;
}

pair<ChrPosition **, uint> Support::getLocations(uint &position_of_locations){
  ChrPosition **locations = NULL;
  position_of_locations = 0;
  ulong *occs_fwd;
  uint nb_occs_fwd=0;
  ulong *occs_rev;
  uint nb_occs_rev=0;
  // multiple is initialized to nb_max_diplayed (--max-locs)
  uint nb_locs_multiple = getGenome()->getNbLocations(); 

  // we search the k-mer with the most occurrences
  for (uint i=1 ; i<getLength() ; i++){ 
    if (nb_locs[i] > nb_locs[position_of_locations])
      position_of_locations = i;
  }

  // First, on forward strand
  if (ranges_forward[position_of_locations].first <= ranges_forward[position_of_locations].second) {
    genome->getOccurrencesFMIndexRange(ranges_forward[position_of_locations].first,
				       ranges_forward[position_of_locations].second,
				       &occs_fwd);
    nb_occs_fwd = ranges_forward[position_of_locations].second - ranges_forward[position_of_locations].first + 1;
  }
  // We check that the number of occurrences is smaller than nb_locs_multiple before the reverse strand treatment
  if (nb_occs_fwd < nb_locs_multiple) {
    if (ranges_reverse[position_of_locations].first <= ranges_reverse[position_of_locations].second) {
      genome->getOccurrencesFMIndexRange(ranges_reverse[position_of_locations].first,
					 ranges_reverse[position_of_locations].second,
					 &occs_rev);
      nb_occs_rev = ranges_reverse[position_of_locations].second - ranges_reverse[position_of_locations].first + 1;
    }
    
    // Again, we check that the total number of occurrences is smaller than nb_locs_multiple
    // and we adjust the nb_locs_multiple if the limit is not reached
    if ((nb_occs_fwd+nb_occs_rev) < nb_locs_multiple)
      nb_locs_multiple=nb_occs_rev+nb_occs_fwd;
  }

  // We save locations
  if (nb_locs_multiple>0){
    locations = (ChrPosition **)malloc(sizeof (ChrPosition *)*nb_locs_multiple);
    
    uint nb_displayed = min(nb_occs_fwd,nb_locs_multiple);
    for (uint i=0 ; i<nb_displayed ; i++){
      locations[i] = getGenome()->getChrPos(occs_fwd[i],1);
    }
    for (uint i=nb_displayed ; i<nb_locs_multiple ; i++){
      locations[i] = getGenome()->getChrPos(occs_rev[i-nb_displayed],-1);
    }
  }else{
    position_of_locations = getTagLength();
  }
  
  // We delete tempory tables
  if (nb_occs_fwd)
    free(occs_fwd);
  if (nb_occs_rev)
    free(occs_rev);
  
  return pair<ChrPosition **, uint>(locations,nb_locs_multiple);
}

uint Support::getNbBreaks() {
  return nb_breaks;
}

uint Support::getNbDuplicate() {
  return nb_duplicate;
}

uint *Support::getNbLocs() {
  return nb_locs;
}

uint Support::getNbLocs(uint i) {
  return nb_locs[i];
}

uint Support::getNbLocsMax() {
  return nb_locs_max;
}

uint Support::getNbMultiple() {
  return nb_multiple;
}

uint Support::getNbPositionsLocated() {
  return nb_pos_located;
}

uint Support::getNbSingle() {
  return nb_single;
}

Parameters *Support::getParameters() {
  return parameters;
}

uint Support::getPositionOfLocation() {
  return (getLocation() == NULL) ? getTagLength() : position_of_location;
}

pair<uint, uint> Support::getKmerRange(uint i, uint strand) {
  if (i<length && i>=0) {
    if (strand == FORWARD_STRAND)
      return ranges_forward[i];
    else
      return ranges_reverse[i];
  } else {
    throw ILLEGAL_STATE_EXCEPTION;
  }
}

int Support::getRepeatLength() {
  if (getEndPosRepeat() == -1 && getStartPosRepeat() == -1)
    return 0;
  else if (getStartPosRepeat() == -1)
    return getEndPosRepeat() + 1;
  else if (getEndPosRepeat() == -1)
    return (int)(getLength() - getStartPosRepeat() + 1);
  else
    return (int)(getEndPosRepeat() -getStartPosRepeat() + 1);
}


int Support::getStartPosRepeat() {
  return start_pos_repeat;
}

uint *Support::getSupport() {
  return support;
}

uint Support::getSupport(uint i) {
  return support[i];
}

char *Support::getTag() {
  return current_tag;
}

uint Support::getTagNum(){
  return tag_num;
}

uint Support::getTagLength() {
  return tag_length;
}

uint Support::getThreshold(){
  return threshold;
}

bool Support::hasRepeat() {
  return (getStartPosRepeat() != -1 || getEndPosRepeat() != -1) 
    && getRepeatLength() >= parameters->percent_min_unique_repetition*getNbPositionsLocated();
}

bool Support::isAlmostNormal(){
  return almostNormal;
}

bool Support::isContinuous() {
  return getLength() == (getNbMultiple() + getNbSingle() + getNbDuplicate());
}

bool Support::isDuplicate() {
  return !isNone() 
    && getNbDuplicate() > 0
    && getNbDuplicate() >= (uint) (parameters->percent_min_duplicate 
                                   * getNbPositionsLocated());
} 

bool Support::isMultiple() {
  return ! isNone() 
    && getNbMultiple() > 0
    && getNbMultiple() >= (uint) (parameters->percent_min_multiple 
				  * getNbPositionsLocated());
}
  
bool Support::isNone() {
  bool is_none = (getNbSingle() == 0 && getNbMultiple() == 0 && getNbDuplicate() == 0);
  if (nb_breaks == 1 && !is_none)
    is_none = (breaks[0]->hasNoEndBreak()
	       && breaks[0]->hasNoStartBreak());
  return is_none;
}

bool Support::isSingle() {
  return !isNone() 
    && getNbSingle() > 0
    && getNbSingle() >= (uint) (parameters->percent_min_single
                                * getNbPositionsLocated());
}

// PRIVATE
void Support::computeBestLocation(bool only_single) {
    // Look for a not isolated 1 in nb_locs
  uint i;
  uint longest_run = 0, current_run = 0;
  uint start_longest_run = 0, start_current_run = 0;
  
  bool in_run = false;
  for (i=0; i < getLength(); i++) {
    if ((only_single && nb_locs[i] == 1)
        || (! only_single && nb_locs[i] > 0)){
      if (in_run) {
	current_run++;
      } else {
	in_run = true;
	start_current_run = i;
	current_run = 1;
      }
    } else if (in_run) {
      in_run = false;
      if (current_run > longest_run) {
	longest_run = current_run;
	start_longest_run = start_current_run;
      }
    }
  }

  if (in_run) {
      if (current_run > longest_run) {
	longest_run = current_run;
	start_longest_run = start_current_run;
      }
  }
  
  // Defines a chrPosition for the support so that we can display at least one
  // genome location (if any).
  if (longest_run > 0) {
    // Best position of the 1 (in the middle of the run)
    uint i = start_longest_run + (longest_run - 1)/2;
    position_of_location = i;
    ulong *occs;
    if (ranges_forward[i].first <= ranges_forward[i].second) {
      genome->getOccurrencesFMIndexRange(ranges_forward[i].first,
					 ranges_forward[i].second,
					 &occs);
      pos = getGenome()->getChrPos(occs[0], 1);
    } else {
      genome->getOccurrencesFMIndexRange(ranges_reverse[i].first,
					 ranges_reverse[i].second,
					 &occs);
      pos = getGenome()->getChrPos(occs[0], -1);
    }
    free(occs);
  } else
    pos = NULL;
}

void Support::tryToMergeBreaks() {
  list<SupportBreak *> break_list;
  list<bool> chimeric_breaks, short_chimeric_breaks;
  list<SupportBreak *>::iterator it_b;
  list<bool>::iterator it_chim;
  list<bool>::iterator it_short;

  for (uint i=0; i < nb_breaks; i++) {
    break_list.push_back(breaks[i]);
  }

  // Filling arrays
  for (uint i = 0; i < nb_breaks; i++) {
    bool value = ! breaks[i]->isNiceBreak()
      && ! breaks[i]->hasNoStartBreak() && ! breaks[i]->hasNoEndBreak();
    chimeric_breaks.push_back(value);
    short_chimeric_breaks.push_back(value && 
                                    (! breaks[i]->hasLongEnoughTagBreak()));
  }

  it_b = break_list.begin();
  it_chim = chimeric_breaks.begin();
  it_short = short_chimeric_breaks.begin();

  
  // Is the first break located?
  // Should we do a loop?
  if ((*it_b)->hasNoStartBreak() && ! (*((++it_b)--))->isNiceBreak()
      && (*((++it_b)--))->getPositionStartBreak() 
      - (*it_b)->getPositionEndBreak() - 1 
      < parameters->min_bases_before_break) {
    merge(break_list, chimeric_breaks, short_chimeric_breaks,    
          it_b, it_chim, it_short);
  }

  // Is the last break located?
  // Should we do a loop?
  it_b = break_list.begin();
  it_chim = chimeric_breaks.begin();
  it_short = short_chimeric_breaks.begin();
  advance(it_b, break_list.size()-1);
  if (break_list.size() > 2 ) {
    advance(it_chim, break_list.size()-2);
    advance(it_short, break_list.size()-2);
  }

  if (break_list.size() > 1 && (*it_b)->hasNoEndBreak() 
      && (*(--it_b)++)->isNiceBreak()
      && (*it_b)->getPositionStartBreak() 
      - (*(--it_b)++)->getPositionEndBreak() - 1 
      < parameters->min_bases_before_break) {
    merge(break_list, chimeric_breaks, short_chimeric_breaks,    
          --it_b, it_chim, it_short);
  }

  // Merging short chimeric breaks
  it_b = break_list.begin();
  it_chim = chimeric_breaks.begin();
  it_short = short_chimeric_breaks.begin();

  while (it_b != break_list.end()) {
    if (*it_short) {
      if (*it_b != break_list.front() 
          && (*it_b)->getPositionEndBreak()
            - (*((--it_b)++))->getPositionStartBreak() + 1
          <= parameters->max_nb_overlapping_breaks * getThreshold()
          && (*((--it_b)++))->isNiceMerge(*it_b)) {
        merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              (--it_b), --it_chim, --it_short);
      } else if (*it_b != break_list.back() 
                 && (*((++it_b)--))->getPositionEndBreak()
                 - (*it_b)->getPositionStartBreak() + 1
                 <= parameters->max_nb_overlapping_breaks * getThreshold()
                 && (*it_b)->isNiceMerge(*(++it_b)--))
        merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              it_b, it_chim, it_short);

      // else if (*it_b != break_list.front() && *(--it_short)++) {
      //   merge(break_list, chimeric_breaks, short_chimeric_breaks,    
      //         --it_b, --it_chim, --it_short);
      // }
      else {
        uint distance_left = 0, distance_right = 0;
        uint short_distance_left = 0, short_distance_right = 0;
        if (*it_b != break_list.front() && *(--it_chim)++) {
          distance_left = (*it_b)->getPositionEndBreak()
            - (*((--it_b)++))->getPositionStartBreak() + 1;
          short_distance_left = (*it_b)->getPositionStartBreak()
            - (*((--it_b)++))->getPositionEndBreak() + 1;
        }
        if (*it_b != break_list.back() && *(++it_chim)--) {
          distance_right = (*((++it_b)--))->getPositionEndBreak()
            - (*it_b)->getPositionStartBreak() + 1;
          short_distance_right = (*((++it_b)--))->getPositionStartBreak()
            - (*it_b)->getPositionEndBreak() + 1;
        }
        if (distance_left 
            && (distance_left 
                <= parameters->max_nb_overlapping_breaks * getThreshold()
                || short_distance_left <= parameters->max_bases_randomly_matched)
            && (distance_left < distance_right || ! distance_right))
          merge(break_list, chimeric_breaks, short_chimeric_breaks,    
                (--it_b), --it_chim, --it_short);
        else if (distance_right 
                 && (distance_right 
                     <= parameters->max_nb_overlapping_breaks * getThreshold()
                     || short_distance_right 
                     <= parameters->max_bases_randomly_matched)
                 && (distance_right < distance_left || ! distance_left))
          merge(break_list, chimeric_breaks, short_chimeric_breaks,    
                it_b, it_chim, it_short);

        else if (*it_b != break_list.front() 
                 && (*((--it_b)++))->hasNoStartBreak()
                 && (*it_b)->getPositionStartBreak()
                 - (*((--it_b)++))->getPositionEndBreak() + 1 
                 <= parameters->max_bases_randomly_matched) {
          merge(break_list, chimeric_breaks, short_chimeric_breaks,    
                (--it_b), --it_chim, --it_short);
        } else if (*it_b != break_list.back() && (*(++it_b)--)->hasNoEndBreak()
                   && (*((++it_b)--))->getPositionStartBreak()
                   - (*it_b)->getPositionEndBreak() + 1 
                   <= parameters->max_bases_randomly_matched)
          merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              it_b, it_chim, it_short);
        else {
          it_b++;
          it_chim++;
          it_short++;
        }       
      } 
    } else {
      it_b++;
      it_chim++;
      it_short++;
    }        
  }

  // for (uint i = 0; i < nb_breaks; i++) {
  //   printf("%d %d %d %d\n",i, chimeric_breaks[i], short_chimeric_breaks[i],
  //          merged_with[i]);
  // }

  // Merging chimeric_breaks
  it_b = break_list.begin();
  it_chim = chimeric_breaks.begin();
  it_short = short_chimeric_breaks.begin();
  while (it_b != break_list.end()) {
    if (*it_chim) {
      if (*it_b != break_list.front()&& (*((--it_b)++))->isNiceMerge(*it_b)
          && ((*it_b)->getPositionEndBreak()
              - (*((--it_b)++))->getPositionStartBreak() + 1 
              <= parameters->max_nb_overlapping_breaks * getThreshold()
              || (*it_b)->getPositionStartBreak()
              - (*((--it_b)++))->getPositionEndBreak() + 1
              <= parameters->max_bases_randomly_matched)) {
        merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              (--it_b), --it_chim, --it_short);
      } else if (*it_b != break_list.back()&& (*it_b)->isNiceMerge(*(++it_b)--)
                 && ((*((++it_b)--))->getPositionEndBreak()
                     - (*it_b)->getPositionStartBreak() + 1 
                     <= parameters->max_nb_overlapping_breaks * getThreshold()
                     || (*((++it_b)--))->getPositionStartBreak()
                     - (*it_b)->getPositionEndBreak() + 1 
                     <= parameters->max_bases_randomly_matched))
        merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              it_b, it_chim, it_short);

      // If we have a chimeric break next to us, we merge it if it is not too
      // far.
      else if (*it_b != break_list.back() && *(++it_chim)--
               && ((*((++it_b)--))->getPositionEndBreak()
                   - (*it_b)->getPositionStartBreak() + 1 
                   <= parameters->max_nb_overlapping_breaks * getThreshold()
                   || (*((++it_b)--))->getPositionStartBreak()
                   - (*it_b)->getPositionEndBreak() + 1 
                   <= parameters->max_bases_randomly_matched))
        merge(break_list, chimeric_breaks, short_chimeric_breaks, 
              it_b, it_chim, it_short);
               
      else if (*it_b != break_list.front()&& (*((--it_b)++))->hasNoStartBreak()
               && (*it_b)->getPositionStartBreak()
               - (*((--it_b)++))->getPositionEndBreak() + 1 
               <= parameters->max_bases_randomly_matched) {
        merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              (--it_b), --it_chim, --it_short);
      } else if (*it_b != break_list.back()&& (*(++it_b)--)->hasNoEndBreak()
                 && (*((++it_b)--))->getPositionStartBreak()
                 - (*it_b)->getPositionEndBreak() + 1 
                 <= parameters->max_bases_randomly_matched)
        merge(break_list, chimeric_breaks, short_chimeric_breaks,    
              it_b, it_chim, it_short);
      else {
        it_b++;
        it_chim++;
        it_short++;
      }
    } else {
      it_b++;
      it_chim++;
      it_short++;
    }
  }

  // Merging other breaks
  it_b = break_list.begin();
  while (*it_b != break_list.back()) {
    if (! (*it_b)->hasNoStartBreak() && ! (*((++it_b)--))->hasNoEndBreak()
        && (*((++it_b)--))->getPositionEndBreak() 
        - (*it_b)->getPositionStartBreak() +1 <= getThreshold()
        && (! (*it_b)->isVeryNiceBreak() || ! (*((++it_b)--))->isVeryNiceBreak())
        && (*it_b)->isVeryNiceMerge(*((++it_b)--))) {
      merge(break_list, chimeric_breaks, short_chimeric_breaks,    
            it_b, it_chim, it_short, true);
    } else if ((*it_b)->getPositionEndBreak() 
               >= (*((++it_b)--))->getPositionStartBreak()) {
      // If the breaks overlap, we merge them.
        merge(break_list, chimeric_breaks, short_chimeric_breaks, 
              it_b, it_chim, it_short, true);
    } else
      it_b++;
  }

  // Just keep the necessary breaks
  it_b = break_list.begin();
  uint current_pos = 0;
  while (it_b != break_list.end()) {
    breaks[current_pos++] = *it_b;
    it_b++;
  }
    
  if (break_list.size() < nb_breaks) {
    // Realloc memory so that we don't take unnecessary memory
    breaks = (SupportBreak **)realloc(breaks, 
                                      break_list.size() * sizeof(SupportBreak *));
    nb_breaks = break_list.size();
  }
}

void Support::merge(list<SupportBreak*> &break_list, list<bool> &chimeric_breaks,
                    list<bool> &short_chimeric_breaks,
                    list<SupportBreak*>::iterator &i, 
                    list<bool>::iterator &chim, 
                    list<bool>::iterator &shortB,
                    bool no_chim) {
  uint pos_break = (*i)->getPositionStartBreak();
  uint pos_end = (*((++i)--))->getPositionEndBreak();
  uint nb_merges = (*i)->getNbMerges() + (*((++i)--))->getNbMerges() + 1;

  // In the case where break i is at the start of the read
  // maybe pos_break has been modified, so we need to find the original one.
  if ((*i)->hasNoStartBreak() && pos_break == 0) {
    while (pos_break < pos_end && getNbLocs()[pos_break] > 0)
      pos_break++;    
  }
  if ((*((++i)--))->hasNoEndBreak() && pos_end == getLength() - 1) {
    while (pos_end > pos_break && getNbLocs()[pos_end] > 0)
      pos_end--;
  }

  // cout << "Merging " << *i << endl;
  // PRINT_VAR(pos_break);
  // PRINT_VAR(pos_end);
  delete (*i);
  i = break_list.erase(i);
  delete (*i);
  i = break_list.erase(i);
  i = break_list.insert(i, new SupportBreak(pos_break, pos_end, this, parameters,
                                            ranges_forward, ranges_reverse,
                                            nb_merges));

  if (! no_chim) {
    bool value = ! (*i)->isNiceBreak()
      && ! (*i)->hasNoStartBreak() 
      && ! (*i)->hasNoEndBreak();

    chim = chimeric_breaks.erase(chim);
    chim = chimeric_breaks.erase(chim);
    chim = chimeric_breaks.insert(chim, value);

    shortB = short_chimeric_breaks.erase(shortB);
    shortB = short_chimeric_breaks.erase(shortB);
    shortB = short_chimeric_breaks.insert(shortB, value 
                                          && (! (*i)->hasLongEnoughTagBreak()));
  }
}
