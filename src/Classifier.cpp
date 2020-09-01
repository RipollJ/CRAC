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
*                   Jérôme AUDOUX    <jerome.audoux@univ-montp2.fr            *
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

#include <Classifier.h>
#include <iostream>
#include <sstream>
#include "utils.h"

using namespace std;

/****************************************
*         SingleReadClassifier          *
****************************************/

SingleReadClassifier::SingleReadClassifier(uint tagId, 
                      LocateOnGenome *genome, 
                      gkArrays *tags, 
                      Parameters *parameters):
  tags(tags),genome(genome),parameters(parameters),tagId(tagId){
  uint *S = tags->getSupport(tagId);
  this->suppo = new Support(parameters, S, tagId, genome, tags);
  this->taginfo = new TagInfo(genome, suppo);
}

SingleReadClassifier::~SingleReadClassifier() {
  if (taginfo != NULL)
    delete taginfo;
  if (suppo != NULL)
    delete suppo;
}

void SingleReadClassifier::classify() {
  char *tag = tags->getTag(tagId);
  int nb_tag_indel;
  int nb_genome_indel;

  if (tagId >= tags->getNbTags()) {
    // TODO return NULL;
  }

  if (! suppo->isContinuous() 
      && ( suppo->isSingle() || suppo->isDuplicate() || suppo->isMultiple() )
      ) {
    // all cases with at least an interruption and where the read is not classify "single" or "duplicate"
    // we have to do an other pass to find the interruption
    for (uint i= 0 ; i < suppo->getNbBreaks() ; i++){
      taginfo->setCurrentBreak(i);
      SupportBreak *suppoBreak = suppo->getBreak(i);
      nb_tag_indel = suppoBreak->getNbTagIndels();
      nb_genome_indel = suppoBreak->getNbGenomeIndels();

      // if (suppoBreak->isDuplicated()){
      // 	cout << "DUPLICATION" << endl;
      // }

      if (! suppoBreak->hasLongEnoughTagBreak()) {
	  taginfo->addUndeterminedError("The break in the tag is too small (break length: %d)",
					suppoBreak->getTagBreakLength());
      } else if (suppoBreak->hasNoCover()) {
	if (suppoBreak->isGenomeInsertion()
	    && parameters->max_bio_ins_del <= labs(nb_genome_indel)
	    && parameters->max_splice_length >= labs(nb_genome_indel)) {
	  // Splice with little support	
	  // Todo: change the category
	  taginfo->addSpliceNoCover(nb_genome_indel);
	} else {
	  taginfo->addUndeterminedError("The read (or the part of read) has no cover, we cannot deduce anything! (score in break:%.2f, score outside break:%.2f)",
					suppoBreak->getScoreInsideBreak(), suppoBreak->getScoreOutsideBreak());
	}
      }
      // Biological event
      else if (suppoBreak->isBiologicalEvent()) {
        // In this case we have a combination of biological reason and error
        // But we can't really say anything on that.
        if (suppoBreak->isDeviated()){
          taginfo->addBioUndetermined(NULL,
                                      suppoBreak->getPositionEndBreak(),
                                      "Combination of error and biological reason (deviation score: %.2f)", 
                                      suppoBreak->getScoreInsideAverages());
        }  else if (suppoBreak->isRepeated()) {
	// ambiguous cases (no limit start or end)
	  taginfo->addBioUndetermined(NULL, suppoBreak->getPositionStartBreak(),
				      "Repetition around the break");
	} 
	else if (suppoBreak->hasNoStartBreak()){ 
	  taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak()
				      ,suppoBreak->getPositionEndBreak()
				      ,"ambiguous case : no start break");
          if (parameters->deep_snp_search) {
            perform_deep_snp_search(suppoBreak, i);
          }
	}	
	else if (suppoBreak->hasNoEndBreak()){
	  taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak()
				      ,suppoBreak->getPositionStartBreak()
				      ,"ambiguous case : no end break");
          if (parameters->deep_snp_search) {
            perform_deep_snp_search(suppoBreak, i);
          }
	}      
	// unambiguous cases
	else {
	  // Check if the candidat is ambiguous or not
	  if (parameters->no_ambiguity && suppoBreak->isDuplicated()){
      string message = "Biological event reclassication because of #no-ambiguity option. Probably several matches of the candidat on the genome.";
      if(suppoBreak->isTagSubstitution()) {
        message += " #snp";
      } else if(suppoBreak->isTagIndel() || suppoBreak->isGenomeDeletion() ||
                (suppoBreak->isGenomeInsertion() && parameters->max_bio_ins_del > labs(nb_genome_indel))) {
        message += " #indel";
      } else if(suppoBreak->isGenomeInsertion() && parameters->max_splice_length > labs(nb_genome_indel)) {
        message += " #splice";
      } else if(suppoBreak->isChimera()) {
        message += " #chimera";
      }
	    taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
					suppoBreak->getPositionStartBreak(),
					message.c_str());
	  }else{
	    // SNP 
	    if (suppoBreak->isTagSubstitution() 
		&& suppoBreak->isBiologicalIntraEvent()){
	      // two substituions
	      if (suppoBreak->getTagBreakLength() > suppo->getThreshold()) {
		taginfo->addSNP(tag[suppoBreak->getPositionEndBreak()
				    -(suppoBreak->getTagBreakLength() - suppo->getThreshold())],
				SECOND_SUBSTITUTION);
	      }
	      
	      // one substitution
	      taginfo->addSNP(tag[suppoBreak->getPositionEndBreak()],
			      FIRST_SUBSTITUTION);  
	    } // end if snp
	    // case biological tag ins/del
	    else if (suppoBreak->isTagIndel()
		     && suppoBreak->isBiologicalIntraEvent()){
	      // snp (insertion) at potsition j-1 into the tag (a snp is one sub, one ins or one del) 
	      if (nb_tag_indel == 1){ 
		taginfo->addSNP(tag[suppoBreak->getPositionEndBreak()], 
				INSERTION);
	      }
	      // insertion at position j-1 into the tag 
	      else if (nb_tag_indel > 0){
		taginfo->addBioTagIndel(nb_tag_indel,0);
	      }
	      // snp (deletion) at position j-1 into the tag 
	      else if (nb_tag_indel == -1){ 
		taginfo->addSNP(0, DELETION);
	      }
	      //deletion of lg ins_del_tag at position j-1 into the tag
	      else {
		taginfo->addBioTagIndel(0,labs(nb_tag_indel));
	      } 
	    } // end if tag indel
	    // case biological deletion in the genome
	    else if (suppoBreak->isGenomeDeletion()
		     && suppoBreak->isBiologicalIntraEvent()){
	      taginfo->addBioTagIndel(labs(nb_genome_indel),0);
	    } // end if genomeDeletion
	    // case biological insertion in the genome
	    else if (suppoBreak->isGenomeInsertion()){
	      // case insertion genome < max_bio_ins_del
	      if (parameters->max_bio_ins_del > labs(nb_genome_indel)){
		taginfo->addBioTagIndel(0,labs(nb_genome_indel));
	      } 
	      // case splice
	      else if (parameters->max_splice_length >= labs(nb_genome_indel)){
		taginfo->addSplice(nb_genome_indel);
	      } else {
		taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
					    suppoBreak->getPositionStartBreak(),
					    "Unknown type of genome insertion");
	      }
	    } // end if genomeInsertion
	    // case inter/intra transplicing (chimera)
	    else if (suppoBreak->isChimera()){
    uint stringentChim = 0;
    string message = "";
	      // If we are using stringent chimera option
	      if(parameters->stringent_chimera) {
		stringentChim = isStringentChimera(suppoBreak,i);
    message = "Chimera reclassication because of stringent-chimera option. Cause is: ";
		switch(stringentChim) {
		case 1:
      message += "break_length";
		  break;
		case 2:
      message += "nb_merges";
		  break;
		case 3:
      message += "repeated_read";
		  break;
		case 4:
      message += "support_variation";
		  break;
		case 5:
      message += "kmer_verification";
		  break;
		case 6:
      message += "no_single_loc";
		  break;
		}
        }
      if(stringentChim > 0) {
		  taginfo->addBioUndetermined(suppoBreak->getLocationStartBreak(),
					      suppoBreak->getPositionStartBreak(),
					      message.c_str());
	      } else {
		taginfo->addSpliceInter();
	      }
	    } // end if chimera
	    // is continuous but no known biological reasons
	    // different cases of biological undetermined reasons
	    else{
	      if (suppoBreak->getTagBreakLength() >= suppo->getThreshold())
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
					    "The break in the tag is too large (%d). There may have several biological causes",
					    suppoBreak->getTagBreakLength());
	      // else if (suppoBreak->getBreakLength() > suppo->getThreshold())
	      //   taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
	      // 				  suppoBreak->getPositionEndBreak(),
	      // 				  "The break is too large (%d). There may have several biological causes",
	      // 				  suppoBreak->getBreakLength());
	      else if (suppoBreak->getTagBreakLength() < parameters->min_break_length)
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
					    "The break is too short (%d).",
					    suppoBreak->getTagBreakLength());
	      else if (! suppoBreak->isBiologicalIntraEvent() 
		       && suppoBreak->isBiologicalInterEvent()) {
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
					    "Support is ok for an inter event but not for an intra.");
	      } else {
		taginfo->addBioUndetermined(suppoBreak->getLocationEndBreak(),
					    suppoBreak->getPositionEndBreak(),
					    "This case was not expected.");
	      }
	    }
	  } 
	}
      } // end if (isBiologicalEvent)
      // --------> start of erronneous reason
      else{
	if (suppoBreak->isRepeated()) {
	  // In this case we have generally a few positions which are not located
	  // But we can't really say anything on that.
	  taginfo->addUndeterminedError("pos=%d-%d Repetitions around the break",
					suppoBreak->getPositionStartBreak(),
					suppoBreak->getPositionEndBreak());
	}
	// ambiguous cases (no limit start or end)
	else if (suppoBreak->hasNoStartBreak()){
	  taginfo->addSeqErr(UNKNOWN_START_MISSING,
	  		     suppoBreak->getPositionEndBreak(),~0,~0,
	  		     suppoBreak->getScoreComputedIntraExon());	    
	}else if(suppoBreak->hasNoEndBreak()){
	  uint pos = suppoBreak->getPositionStartBreak() + suppo->getThreshold() - 1;
	  taginfo->addSeqErr(UNKNOWN_END_MISSING,
			     pos >= tags->getTagLength(i) ? tags->getTagLength(i)-1 : pos,~0,~0,
			     suppoBreak->getScoreComputedIntraExon());	    
	}
	// seq error ---------------------------> (unambiguous choice)
	else if (suppoBreak->isTagSubstitution()){
	  //two errors
	  if (suppoBreak->getTagBreakLength() > suppo->getThreshold()) {
	    uint shift = suppoBreak->getTagBreakLength() - suppo->getThreshold();
	    taginfo->addSeqErr(SECOND_SUBSTITUTION
			       ,suppoBreak->getPositionEndBreak()-shift,0,0,
			       suppoBreak->getScoreComputedIntraExon(),
			       1,
			       getTagAtEndBreak(suppoBreak, tag, 1, -shift), 1);
	  }
	  // one error
	  taginfo->addSeqErr(FIRST_SUBSTITUTION,
			     suppoBreak->getPositionEndBreak(),0,0,
			     suppoBreak->getScoreComputedIntraExon(),
			     1, 
			     getTagAtEndBreak(suppoBreak, tag, 1), 1);
	  
	}
	// case ins/del error into tag
	else if (suppoBreak->isTagIndel()) {
	  // insertion of lg ins_del_tag 
	  if (nb_tag_indel >= 0){
	    taginfo->addSeqErr(INSERTION,
			       suppoBreak->getPositionEndBreak(),
			       nb_tag_indel,0,
			       suppoBreak->getScoreComputedIntraExon(),
                               0,
			       getTagAtEndBreak(suppoBreak, tag, nb_tag_indel),
			       nb_tag_indel);
	  }
	  //suppression of lg ins_del_tag 
	  else{
	    taginfo->addSeqErr(DELETION
			       ,suppoBreak->getPositionEndBreak(),
			       0,labs(nb_tag_indel),
			       suppoBreak->getScoreComputedIntraExon(),
			       labs(nb_tag_indel), NULL, ~0);
	  } 
	}	
	//insertion of lg ins_del_genome 
	else if (suppoBreak->isGenomeInsertion()) {
	  if (parameters->max_bio_ins_del < labs(nb_genome_indel)
	      && parameters->max_splice_length >= labs(nb_genome_indel)
              // The following condition is there just to select
              // splices with no cover.
              // Otherwise it may be a false positive due to duplications
              // and an error
	      && suppoBreak->getScoreOutsideBreak() <= parameters->max_support_out_no_cover) {
	    // Splice with little support	
	    taginfo->addSpliceNoCover(nb_genome_indel);
	  }else {
            ChrPosition *chrPos = suppoBreak->getLocationEndBreak();
            taginfo->addUndeterminedError("position %d. Too many indels to be an error (probably a splice gap=%d, location=%s|%d,%u).", suppoBreak->getPositionEndBreak(), nb_genome_indel,suppoBreak->getChr(END_BREAK),
                                          suppoBreak->getStrandLocation(END_BREAK), chrPos->getRelativePosition());
            delete chrPos;
          }
	}
	//suppresion of lg labs(ins_del_genome) 
	else if (suppoBreak->isGenomeDeletion()) { 
	  taginfo->addSeqErr(INSERTION,
                             suppoBreak->getPositionEndBreak()
			     ,labs(nb_genome_indel),0,
			     suppoBreak->getScoreComputedInterExon(),
			     0,
			     getTagAtEndBreak(suppoBreak, tag, labs(nb_genome_indel)),
			     labs(nb_genome_indel));
	}
	// is not continuous but no erroneous reason
	else{
	  taginfo->addUndeterminedError("position = %d -> No location neither biological nor erroneous cause found (ins_del_tag: %d and ins_del_genome: %d, chr_before_break: %d and  chr_after_break: %d, strand_before_break: %d and strand_after_break: %d)"
					,suppoBreak->getPositionEndBreak()
					,nb_tag_indel
					,nb_genome_indel
					,suppoBreak->getChrId(START_BREAK)
					,suppoBreak->getChrId(END_BREAK)
					,suppoBreak->getStrandLocation(START_BREAK)
					,suppoBreak->getStrandLocation(END_BREAK));
	}
      }
    }
    taginfo->setCurrentBreak(suppo->getNbBreaks());
  }
  delete [] tag;
}

uint SingleReadClassifier::getTagId() {
  return tagId;
}

TagInfo *SingleReadClassifier::getTagInfo() {
  return taginfo;
}

readIterator *SingleReadClassifier::setInfos(readIterator *readIt) {
  taginfo->setReadName(readIt->getName());
  taginfo->setReadQuality(readIt->getQuality());
  taginfo->setReadNumber(readIt->getReadNumber());  
  ++(*readIt);
  return readIt;
}

ostream &SingleReadClassifier::samOutput(ostream &os) {
  taginfo->samOutput(os);
  return os;
}

void SingleReadClassifier::writeOutputs(ostream *snp
      , ostream *bioTagIndel
      , ostream *seqErr 
      , ostream *splice
      , ostream *spliceNoCover
      , ostream *spliceInter
      , ostream *undetermined
      , ostream *repetition
      , ostream *duplication
      , ostream *nothing
      , ostream *normal
      , ostream *almostNormal
      , ostream *multiple
      , ostream *none
      , ostream *bioUndermined
      , ostream *single) {

  uint i = taginfo->getReadNumber();

  // First step : the mapping output
  if (taginfo->hasNothing()){
    if (nothing != NULL)
      *nothing << i << " " << *taginfo <<  endl;
  }else{
    if (single != NULL)
      if (taginfo->isSingle())
        *single << i << " " << *taginfo << endl;
    if (duplication != NULL)
      if (taginfo->isDuplication())
        *duplication << i << " " << *taginfo <<  endl ;
    if (multiple != NULL)
      if (taginfo->isMultiple())
        *multiple << i << " " << *taginfo << endl;
    if (none != NULL)
      if (taginfo->isNone())
        *none << i << " " << *taginfo << endl;
  }
    
  //Second step : the classification output
  // if it is an erroneous tag then we uniquely write in the output errors file
  if (taginfo->getNbSeqErr() > 0) {
    if (seqErr != NULL)
      for (uint j=0; j < taginfo->getNbSeqErr(); j++) 
        *seqErr << i << " " << *taginfo->getInfosSeqErr()[j] << " " << *taginfo <<  endl ;
  }
  //     else{
  if (normal != NULL)
    if (taginfo->isNormal())
      *normal << i << " " << *taginfo <<  endl;
  if (almostNormal != NULL)
    if (taginfo->isAlmostNormal())
      *almostNormal << i << " " << *taginfo <<  endl;
  if (snp  != NULL)
    for (uint j=0; j < taginfo->getNbSNP(); j++) 
      *snp << i << " " << *taginfo->getInfosSNP()[j] << " " << *taginfo <<  endl;
  if (bioTagIndel  != NULL)
    for (uint j=0; j < taginfo->getNbBioTagIndel(); j++) 
      *bioTagIndel << i << " " << *taginfo->getInfosBioTagIndel()[j] << " " << *taginfo <<  endl;
  if (repetition != NULL)
    for (uint j=0; j < taginfo->getNbRepetition() ; j++) 
      *repetition << i << " " << *taginfo->getInfosRepetition()[j] << " " << *taginfo <<  endl;
  if (undetermined != NULL)
    for (uint j=0; j < taginfo->getNbUndeterminedError(); j++) 
      *undetermined << i << " " << *taginfo->getInfosUndeterminedError()[j] << " " << *taginfo <<  endl;
  if (splice != NULL)
    for (uint j=0; j < taginfo->getNbSplice(); j++)
      *splice << i << " " << *taginfo->getInfosSplice()[j] << " " << *taginfo <<  endl ;
  if (spliceNoCover != NULL)
    for (uint j=0; j < taginfo->getNbSpliceNoCover(); j++)
      *spliceNoCover << i << " " << *taginfo->getInfosSpliceNoCover()[j] << " " << *taginfo <<  endl;
  if (spliceInter != NULL)
    for (uint j=0; j < taginfo->getNbSpliceInter(); j++)
      *spliceInter << i << " " << *taginfo->getInfosSpliceInter()[j] << " " << *taginfo <<  endl;
  if (bioUndermined != NULL)
    for (uint j=0; j < taginfo->getNbBioUndetermined(); j++)
      *bioUndermined << i << " " << *taginfo->getInfosBioUndetermined()[j] << " " << *taginfo << endl ;
}

void SingleReadClassifier::updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable) {
  int code = taginfo->getCode();
  for (uint j=0; j < NB_MASKS; j++) {
    if (code & (1<<j))
      (*nb_classes)[j]++;
  }
  if (taginfo->isExplainable())
    (*nb_explainable)++;
  if (code == 0)
    (*nb_classes)[NB_MASKS]++;
}

/// PRIVATE ///
char *SingleReadClassifier::getTagAtEndBreak(SupportBreak *suppoBreak, 
                            char *tag,
				                    uint length, 
                            uint shift) {
  if (length > parameters->max_bases_retrieved) 
    return NULL;

  uint pos = suppoBreak->getPositionEndBreak();
  char *dna = new char[length+1];
  dna[length] = 0;
  strncpy(dna, &tag[pos+shift-length+1], length);
  return dna;
}

void SingleReadClassifier::perform_deep_snp_search(SupportBreak *sbreak, uint break_num) {
  
  if (! sbreak->isBiologicalEvent() 
      || sbreak->getScoreComputedIntraExon() > parameters->p_value_variation_biological)
    return;
  // For now, we currently search substitutions.

  uint loc_break;               // location of the k-mer before or after the break
  uint loc_snp;                 // location of the putative SNP
  // number of nucleotides for the comparison (the minimum is specified here, 
  // but that can be larger)
  uint nb_nuc = parameters->number_nucleotides_snp_comparison; 
  
  uchar *test_qmer;      // test q-mer from the read that (may) contain
  // the SNP (the SNP is at the end of the q-mer in the case of a start break
  // and at the start of the q-mer in the case of an end break)
  uint loc_test_qmer;           /* putative location of the test q-mer on the 
                                 * genome */
  uchar *real_content;      // real DNA from the genome, including the SNP that should correspond to test_qmer
  uint pos_snp_in_qmer;         // position of the putative SNP in the test q-mer
  uint pos_of_test_qmer;        // position of the modified q-mer in the read
  int strand = 0;

  // Variables used to update the break's informations
  uint pos_start_break;
  uint pos_end_break;
  uint loc_start_break;
  uint loc_end_break;

  const char *read = taginfo->getTag();
  uint read_length = taginfo->getSupportObject()->getTagLength();
  if (sbreak->hasNoStartBreak() 
      && sbreak->getPositionEndBreak() >= nb_nuc ) {
    
    nb_nuc = max(nb_nuc, sbreak->getPositionEndBreak());
    test_qmer = new uchar[nb_nuc+2];

    // End break's most likely location
    // SNP should be one position before (on fwd strand)
    pos_end_break = sbreak->getPositionEndBreak();
    loc_end_break = loc_break = taginfo->getLocationEndBreak(break_num);
    strand = taginfo->getChrPosEndBreak(break_num)->getStrand();
    if (strand == -1) {
      loc_end_break += taginfo->getThreshold() - nb_nuc;
    }
    pos_of_test_qmer = sbreak->getPositionEndBreak() - nb_nuc;
    pos_start_break = pos_of_test_qmer + 1;
    strncpy((char *)test_qmer,
            &read[pos_of_test_qmer],
            nb_nuc+1);
    pos_snp_in_qmer = nb_nuc;
    if (strand == 1) {
      loc_test_qmer = loc_break - (1 + nb_nuc);
      loc_snp = loc_break - 1;
      loc_start_break = loc_test_qmer; // starts at the same position as the test-qmer
    } else {
      loc_test_qmer = loc_break + taginfo->getThreshold();
      loc_snp = loc_test_qmer;
      loc_start_break = loc_test_qmer + 1; // ends one position before the test q-mer

    }
  } else if (sbreak->hasNoEndBreak()
             && sbreak->getPositionStartBreak() + taginfo->getThreshold() + nb_nuc <= read_length) {
    nb_nuc = read_length - (sbreak->getPositionStartBreak() + taginfo->getThreshold());
    test_qmer = new uchar[nb_nuc+2];

    loc_break = taginfo->getLocationStartBreak(break_num);
    pos_of_test_qmer = sbreak->getPositionStartBreak() + taginfo->getThreshold() - 1;
    pos_end_break = pos_of_test_qmer;
    pos_start_break = pos_end_break - nb_nuc + 1;
    strand = taginfo->getChrPosStartBreak(break_num)->getStrand();
    loc_start_break = loc_break + strand*(pos_start_break - sbreak->getPositionStartBreak());
    if (strand == -1)
      // If we are on the reverse strand, we must correct the position since the threshold changed
      // (on the reverse strand, the position depends on the threshold)
      loc_start_break += taginfo->getThreshold() - nb_nuc; 

    strncpy((char *)test_qmer,
            &read[pos_of_test_qmer],
            nb_nuc+1);
    pos_snp_in_qmer = 0;

    if (strand == 1) {
      loc_test_qmer = loc_break + taginfo->getThreshold();
      loc_snp = loc_test_qmer;
      loc_end_break = loc_test_qmer + 1; // starts one position after the snp
    } else {
      loc_test_qmer = loc_break - 1 - nb_nuc;
      loc_snp = loc_break - 1;
      loc_end_break = loc_test_qmer; // ends at the same pos as the q-mer with the SNP
    }
  } else {
    return;
  }

  test_qmer[nb_nuc + 1]=0;
  real_content = genome->getGenomeSubstring(loc_test_qmer, nb_nuc+1, strand);
  
  snp_type type = SNP_UNKNOWN;
  bool match = true;
  for (uint i=0, j=0; match && i <= nb_nuc && j <= nb_nuc; i++, j++) {
    if (real_content[i] != test_qmer[j]) {
      if (i == j && i == pos_snp_in_qmer)
        type = SNP_SUBSTITUTION;
      else
        match = false;
    }
  }

  ChrPosition *pchrPos = (genome->getChrPos(loc_snp, strand));
  ChrPosition chrPos(*pchrPos);
  delete pchrPos;
  if (type == SNP_SUBSTITUTION && match) {
    taginfo->changeGenericElement(INFO_BIOLOGICAL_UNDETERMINED, INFO_SNP,
                                  MASK_BIOLOGICAL_UNDETERMINED, MASK_SNP,
                                  break_num, 0,
                                  new SNPInfo(chrPos, pos_of_test_qmer + pos_snp_in_qmer, 
                                              sbreak->getScoreComputedIntraExon(),
                                              taginfo->isDuplicated(break_num),
                                              test_qmer[pos_snp_in_qmer], 
                                              real_content[pos_snp_in_qmer]));
    // Update break's informations
    sbreak->setPositionStartBreak(pos_start_break);
    sbreak->setPositionEndBreak(pos_end_break);
    sbreak->setThreshold(nb_nuc);
    sbreak->setLocationStartBreak(loc_start_break);
    sbreak->setLocationEndBreak(loc_end_break);
    sbreak->setStrandStartBreak(strand);
    sbreak->setStrandEndBreak(strand);
    sbreak->computeScore(pos_start_break, pos_end_break);
  }
  delete [] test_qmer;
  free(real_content);
}

uint SingleReadClassifier::isStringentChimera(SupportBreak *sbreak, uint break_num) {

  // TODO add tests for strigent-chimera option
  uint is_strigent_chimera = 0;

  // Check the size of support with higer expectations than usual
  if (sbreak->getTagBreakLength() < parameters->stringent_chimera_min_break_length) {
    is_strigent_chimera = 1;
  }

  // Check if the break where the chimera is from has a number of merges < to the related parameter
  if(is_strigent_chimera == 0 && sbreak->getNbMerges() > parameters->stringent_chimera_max_number_of_merges) {
    is_strigent_chimera = 2;
  }

  // Check if the break is in a repeated region
  // TODO naive check, we should try to know how far is the repeat region
  // from the chimera and add a min length parameter in const.h
  if (is_strigent_chimera == 0 && taginfo->getSupportObject()->hasRepeat()) {
    is_strigent_chimera = 3;
  }

  // Check interquartile range of support values in the break
  if(is_strigent_chimera == 0 && (float)(sbreak->getInsideQuartile4() - sbreak->getInsideQuartile1()) > sbreak->getAverageHighInside()*PERCENT_SUPPORT_VARIATION_ALMOST_NORMAL) {
    is_strigent_chimera = 4;
  }

  // Check if we found a single loc before and after the break within a window of k-1
  if(is_strigent_chimera == 0){ 
    int first_single_loc_before_break, first_single_loc_after_break;
    uint startBreak = sbreak->getPositionStartBreak();
    uint endBreak = sbreak->getPositionEndBreak();
    int max_ext_length = parameters->threshold - (endBreak - startBreak + 1) - 1;

    if(max_ext_length > 0) {
      // First we try to find a single loc before the break
      first_single_loc_before_break = startBreak - 1;
      while (first_single_loc_before_break >= 0 &&
             suppo->getNbLocs()[first_single_loc_before_break] != 1 &&
             (int) (startBreak - first_single_loc_before_break) <= max_ext_length) {
        first_single_loc_before_break--;
      }

      // First we try to find a single loc before the break
      first_single_loc_after_break = endBreak + 1;
      while (first_single_loc_after_break < (int) suppo->getLength() &&
             suppo->getNbLocs()[first_single_loc_after_break] != 1 &&
             (int) (first_single_loc_after_break - endBreak) <= max_ext_length) {
        first_single_loc_after_break++;
      }

      if(first_single_loc_before_break < 0 || first_single_loc_after_break >= (int) suppo->getLength() ||
         (first_single_loc_after_break - first_single_loc_before_break) >= (int) parameters->threshold) {
        is_strigent_chimera = 6;
      }
    }
  }

  // Check a k_mer location profile after and before the chimera
  if(is_strigent_chimera == 0) {
    uint start_break = sbreak->getPositionStartBreak();
    uint end_break = sbreak->getPositionEndBreak();
    Support *suppo = taginfo->getSupportObject();
    int before_pos = start_break - parameters->max_extension_length * 2;
    bool loc_match_found = false;
    uint nb_locs;

    // if the k-mer before the chimera is in the support
    if ( before_pos >= 0 && before_pos < (int) suppo->getLength()) {
      if(suppo->getNbLocs()[before_pos] > 0) {
        ChrPosition *start_break_pos = sbreak->getChrPosition(START_BREAK);
        if(start_break_pos) {
          ChrPosition *before_chrPos;
          pair <uint, uint> before_range;
          ulong *locs_before_pos;

          if(start_break_pos->getStrand() == 1) {
            before_range = suppo->getKmerRange(before_pos,FORWARD_STRAND);
          }
          else {
            before_range = suppo->getKmerRange(before_pos,REVERSE_STRAND);
          }

          genome->getOccurrencesFMIndexRange(before_range.first,
                                               before_range.second,
                                               &locs_before_pos);

          nb_locs = min(before_range.second - before_range.first + 1,genome->getNbLocations());
          uint i = 0;
          while(!loc_match_found && i < nb_locs) {
            before_chrPos = genome->getChrPos(locs_before_pos[i],start_break_pos->getStrand());
            long long int relativeDist = labs((long long int)before_chrPos->getRelativePosition() - (long long int)start_break_pos->getRelativePosition());
            if (((before_chrPos->getRelativePosition() < start_break_pos->getRelativePosition()) == (start_break_pos->getStrand() == 1))
                && strcmp(before_chrPos->getChrPosition(),start_break_pos->getChrPosition()) == 0
                && relativeDist < parameters->max_splice_length) {
              loc_match_found = true;
            } 
            delete before_chrPos;
            i++;
          }
          free(locs_before_pos);
          if (!loc_match_found)
            is_strigent_chimera = 5;
        }
        if (start_break_pos) 
          delete start_break_pos;
      }
    }
    // if the k-mer after the chimera is still in the support
    if (is_strigent_chimera == 0) {
      int after_pos = end_break + parameters->max_extension_length * 2;
      if (after_pos < (int) suppo->getLength() && after_pos >= 0) {
        if(suppo->getNbLocs()[after_pos] > 0) {
          ChrPosition *end_break_pos = sbreak->getChrPosition(END_BREAK);
          if(end_break_pos) {
            ChrPosition *after_chrPos;
            pair <uint, uint> after_range;
            ulong *locs_after_pos;

            if(end_break_pos->getStrand() == 1) {
              after_range = suppo->getKmerRange(after_pos,FORWARD_STRAND);
            }
            else {
              after_range = suppo->getKmerRange(after_pos,REVERSE_STRAND);
            }

            genome->getOccurrencesFMIndexRange(after_range.first,
                                                 after_range.second,
                                                 &locs_after_pos);

            nb_locs = min(after_range.second - after_range.first + 1,genome->getNbLocations());
            uint i = 0;
            loc_match_found = false;
            while(!loc_match_found && i < nb_locs) {
              after_chrPos = genome->getChrPos(locs_after_pos[i],end_break_pos->getStrand());
              long long int relativeDist = labs((long long int)after_chrPos->getRelativePosition() - (long long int)end_break_pos->getRelativePosition());
              if (((after_chrPos->getRelativePosition() > end_break_pos->getRelativePosition()) == (end_break_pos->getStrand() == 1))
                  && strcmp(after_chrPos->getChrPosition(),end_break_pos->getChrPosition()) == 0
                  && relativeDist < parameters->max_splice_length) {
                loc_match_found = true;
              } 
              delete after_chrPos;
              i++;
            }
            free(locs_after_pos);
            if (!loc_match_found)
              is_strigent_chimera = 5;
          }
          if (end_break_pos)
            delete end_break_pos;
        }
      }
    }
  }

  return is_strigent_chimera;
}

/****************************************
*       PairedEndReadClassifier         *
****************************************/

PairedEndReadClassifier::PairedEndReadClassifier(uint tagId1, uint tagId2,
                         LocateOnGenome *genome, 
                         gkArrays *tags, 
                         Parameters *parameters) {
  this->tags = tags;
  this->genome = genome;
  this->parameters = parameters;
  this->paired_end_chimera = 0;
  classifier1 = new SingleReadClassifier(tagId1,genome,tags,parameters);
  classifier2 = new SingleReadClassifier(tagId2,genome,tags,parameters);
}

PairedEndReadClassifier::~PairedEndReadClassifier() {
  delete classifier1;
  delete classifier2;
}

void PairedEndReadClassifier::classify() {
  classifier1->classify();
  classifier2->classify();

  TagInfo *taginfo1 = classifier1->getTagInfo();
  TagInfo *taginfo2 = classifier2->getTagInfo();

  // reclassement of the chimera(s) in the first read
  uint i = 0;
  while (i < taginfo1->getNbSpliceInter()) {
    if(!chimeraPairedEndCheck(taginfo1->getInfosSpliceInter()[i],classifier2)) {
      taginfo1->removeSpliceInter(i,"Chimera reclassification because paired-end informations does not match the chimera.");
    } else {
      i++;
    }
  }

  // reclassement of the chimera(s) in the second read
  i = 0;
  while (i < taginfo2->getNbSpliceInter()) {
    if(!chimeraPairedEndCheck(taginfo2->getInfosSpliceInter()[i],classifier1)) {
      taginfo2->removeSpliceInter(i,"Chimera reclassification because paired-end informations does not match the chimera.");
    } else {
      i++;
    }
  }
  
  // Looking for paired-end chimera
  if(taginfo1->isSingle() && taginfo2->isSingle() && !taginfo1->hasSpliceInterChr() && !taginfo2->hasSpliceInterChr()) {
    ChrPosition *paire1 = taginfo1->getLocation();
    ChrPosition *paire2 = taginfo2->getLocation();

    if(paire1 != NULL && paire2 != NULL) {

      if(strcmp(paire1->getChrPosition(), paire2->getChrPosition()) != 0) {
        paired_end_chimera = 1;
      } else {
        if(paire1->getStrand() != paire2->getStrand()) {
          long long int relativeDist = labs((long long int)paire1->getRelativePosition() - (long long int)paire2->getRelativePosition());
          if((paire1->getStrand() == 1) && (paire1->getRelativePosition() > paire2->getRelativePosition() && relativeDist > taginfo2->getSupportLength())) {
            paired_end_chimera = 3;
          } else if((paire2->getStrand() == 1) && (paire2->getRelativePosition() > paire1->getRelativePosition() && relativeDist > taginfo1->getSupportLength())) {
            paired_end_chimera = 3;
          } else if(relativeDist > parameters->max_splice_length) {
            paired_end_chimera = 2;
          }
        } else { // Paired-end read should be located on opposite strand according to the protocol
          paired_end_chimera = 4;
        }
      }
    }
  }
}

ostream &PairedEndReadClassifier::samOutput(ostream &os) {
  vector<SamLine*> *sam_lines1 = getFirstTagInfo()->getSamLines();
  vector<SamLine*> *sam_lines2 = getSecondTagInfo()->getSamLines();

  SamLine *primary_line1 = (*sam_lines1)[0];
  SamLine *primary_line2 = (*sam_lines2)[0];


  // If PE read is mapped we print paired-end additional fields
  if(!primary_line2->isSegmentUnmapped()) {
    ostringstream additionalInfos;
    additionalInfos << "loc:" << getSecondTagInfo()->isSingle() << ':' << getSecondTagInfo()->isDuplication() << ':' << getSecondTagInfo()->isMultiple();
    // If there is a chimera between both reads, we print that
    // piece of information
    if(this->hasPairedEndChimera()) {
      additionalInfos << ";chimera:" << *getFirstTagInfo()->getLocation() << ':' << *getSecondTagInfo()->getLocation();
    }
    primary_line1->addOptionalField("XP",additionalInfos.str());
  }


  setPairedEndOptionalFields(*primary_line1,*primary_line2);
  setPairedEndOptionalFields(*primary_line2,*primary_line1);

  writeSamLines(os,*sam_lines1,*primary_line1,*primary_line2,true);
  writeSamLines(os,*sam_lines2,*primary_line2,*primary_line1,false);

  for (vector<SamLine*>::iterator it = sam_lines1->begin() ; it != sam_lines1->end(); ++it)
    delete (*it);
  delete sam_lines1;
  for (vector<SamLine*>::iterator it = sam_lines2->begin() ; it != sam_lines2->end(); ++it)
    delete (*it);
  delete sam_lines2;

  return os;
}

void PairedEndReadClassifier::setPairedEndOptionalFields(SamLine &line, const SamLine &paired_line) {
  ostringstream string_stream;

  line.addOptionalField("R2",paired_line.getSeq());
  string_stream << paired_line.getCigar();
  line.addOptionalField("MC",string_stream.str());
  line.addOptionalField("MQ",paired_line.getMapQ());
}

void PairedEndReadClassifier::writeSamLines(ostream &os, vector<SamLine*> &sam_lines, SamLine &primary_line, SamLine &paired_primary_line, bool is_first_taginfo) {
  // We add informations about paired-end in the SAM lines
  // Setting flag bits for paired-end
  for (vector<SamLine*>::iterator it = sam_lines.begin() ; it != sam_lines.end(); ++it) {
    // Set FLAG
    if(is_first_taginfo) {
      (*it)->setFirstSegmentInTheTemplate();
    } else {
      (*it)->setLastSegmentInTheTemplate();
    }
    (*it)->setTemplateHavingMultipleSegments();
    if(primary_line.isSegmentUnmapped() || paired_primary_line.isSegmentUnmapped()) {
      (*it)->unsetEachSegmentsMapped();
    } else {
      (*it)->setEachSegmentsMapped();
    } 
    if(paired_primary_line.isSegmentUnmapped())
      (*it)->setNextSegmentUnmapped();
    if(paired_primary_line.isSeqReverseComplemented())
      (*it)->setNextSeqReverseComplemented();
    // Set Pnext (field 7)
    (*it)->setPnext(paired_primary_line.getPos());
    // Set Rnext (field 8)
    if((*it)->getRname() == paired_primary_line.getRname()) {
      (*it)->setRnext("=");
    } else {
      (*it)->setRnext(paired_primary_line.getRname());
    }
    // Set TLEN (field 9)
    // if both segment are mapped to the same reference
    if(!primary_line.isSegmentUnmapped() && !paired_primary_line.isSegmentUnmapped()
        && primary_line.getRname() == paired_primary_line.getRname()) {
      int tlen = 0;
      if(primary_line.getPos() < paired_primary_line.getPos()) {
        tlen = (paired_primary_line.getPos() + paired_primary_line.getCigar().getReferenceAlignementLength()) - primary_line.getPos();
        // Confused definition of tlen, should we compute a genomic distance
        // or an exonic distance?
        //tlen = (paired_primary_line.getPos() + paired_primary_line.getCigar().getNbDeletions())
        //  - (primary_line.getPos() + primary_line.getCigar().getNbPositionsRead());
      } else if (primary_line.getPos() > paired_primary_line.getPos()) {
        tlen = paired_primary_line.getPos() - (primary_line.getPos() + primary_line.getCigar().getReferenceAlignementLength());
        //tlen = (paired_primary_line.getPos() + paired_primary_line.getCigar().getNbPositionsRead())
        //  - (primary_line.getPos() + primary_line.getCigar().getNbDeletions());
      }
      (*it)->setTlen(tlen);
      //primary_line1->setTlen(primary_line1->getPos() - (primary_line2->getPos() + primary_line2->getCigar().getReferenceAlignementLength())
      //primary_line2->setTlen(primary_line2->getPos() - (primary_line1->getPos() + primary_line1->getCigar().getReferenceAlignementLength())
    }
    (*it)->writeLine(os);
  }
}

readIterator *PairedEndReadClassifier::setInfos(readIterator *readIt) {
  classifier1->setInfos(readIt);
  classifier2->setInfos(readIt);
  return readIt;
}


void PairedEndReadClassifier::writeOutputs(ostream *snp
      , ostream *bioTagIndel
      , ostream *seqErr 
      , ostream *splice
      , ostream *spliceNoCover
      , ostream *spliceInter
      , ostream *undetermined
      , ostream *repetition
      , ostream *duplication
      , ostream *nothing
      , ostream *normal
      , ostream *almostNormal
      , ostream *multiple
      , ostream *none
      , ostream *bioUndermined
      , ostream *single) {

  classifier1->writeOutputs(snp,
      bioTagIndel,
      seqErr,
      splice,
      spliceNoCover,
      spliceInter,
      undetermined,
      repetition,
      duplication,
      nothing,
      normal,
      almostNormal,
      multiple,
      none,
      bioUndermined,
      single);

  classifier2->writeOutputs(snp,
      bioTagIndel,
      seqErr,
      splice,
      spliceNoCover,
      spliceInter,
      undetermined,
      repetition,
      duplication,
      nothing,
      normal,
      almostNormal,
      multiple,
      none,
      bioUndermined,
      single);
}

void PairedEndReadClassifier::updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable) {
  classifier1->updateStatistics(nb_classes, nb_explainable);
  classifier2->updateStatistics(nb_classes, nb_explainable);
  if(this->hasPairedEndChimera()) {
    (*nb_classes)[posBitInMask(MASK_PAIRED_END_CHIMERA)]++;
  }
}

TagInfo *PairedEndReadClassifier::getFirstTagInfo() {
  return classifier1->getTagInfo();
}

TagInfo *PairedEndReadClassifier::getSecondTagInfo() {
  return classifier2->getTagInfo();
}

bool PairedEndReadClassifier::hasPairedEndChimera() {
  return paired_end_chimera > 0;
}

void PairedEndReadClassifier::writePairedEndChimera(ostream *pairedEndChimera) {
  // If there is a chimera between the paired-end tags
  if(this->hasPairedEndChimera()) {
    *pairedEndChimera << paired_end_chimera << " " << getFirstTagInfo()->getReadNumber() << " " << *getFirstTagInfo() << " " << getSecondTagInfo()->getReadNumber() << " " << *getSecondTagInfo() << endl;
  }
}

/// PRIVATE ///
bool PairedEndReadClassifier::chimeraPairedEndCheck(SpliceInterInfo *chimera,
                              SingleReadClassifier *pairedEndClassifier) {
  TagInfo *pairedEndTag = pairedEndClassifier->getTagInfo();
  ChrPosition P1 = chimera->getChrPosition(); // Position of the chimera before the break
  ChrPosition P2 = chimera->getChromosomeDest(); // Position of the chimera before the break
  ChrPosition *pairePos = NULL;
  
  if (pairedEndTag->isSingle())
    pairePos = pairedEndTag->getLocation(); // Position of the paired-end read only if it is the once (otherwise it is crap!!!)

  bool isAValidChimera = false;

  // If the paired-end read is located and has no chimera
  if (pairePos != NULL && !pairedEndTag->hasSpliceInterChr()) {
    // ------------------------------------
    // 1. We use P1 position of the chimera
    // ------------------------------------
    // a) P1(+) Case
    //                                        ^ (P1)
    //                                 -------|\\\\\|--------- strand : +
    // ======================================================= (ADN)
    // -------------------- (paired-end read)                  strand : -
    if (P1.getStrand() == 1) {
      if (P1.getStrand() != pairePos->getStrand() &&
          strcmp(P1.getChrPosition(), pairePos->getChrPosition()) == 0 &&
          P1.getRelativePosition() >= pairePos->getRelativePosition() &&
          (P1.getRelativePosition() - pairePos->getRelativePosition()) <= parameters->max_splice_length) {
        isAValidChimera = true;
      }
    } else if (!isAValidChimera) {
    // b) P1(-) Case
    //                  (paired-end read) -------------------- strand : +
    // ======================================================= (ADN)
    // -------|\\\\\|---------                                 strand : -
    //              ^ (P1)
      if (P1.getStrand() != pairePos->getStrand() &&
          strcmp(P1.getChrPosition(), pairePos->getChrPosition()) == 0 &&
          P1.getRelativePosition() <= pairePos->getRelativePosition() &&
          (pairePos->getRelativePosition() - P1.getRelativePosition()) <= parameters->max_splice_length) {
        isAValidChimera = true;
      }
    }
    // ------------------------------------
    // 2. We use P2 position of the chimera
    // ------------------------------------
    // a) P2(+) Case
    //              ^ (P2)
    // -------|\\\\\|---------                                 strand : +
    // ======================================================= (ADN)
    //                  (paired-end read) -------------------- strand : -
    if (P2.getStrand() == 1) {
      if (P2.getStrand() != pairePos->getStrand() &&
          strcmp(P2.getChrPosition(), pairePos->getChrPosition()) == 0 &&
          P2.getRelativePosition() <= pairePos->getRelativePosition() &&
          (pairePos->getRelativePosition() - P2.getRelativePosition()) <= parameters->max_splice_length) {
        isAValidChimera = true;
      }
    } else if (!isAValidChimera) {
    // b) P2(-) Case
    // -------------------- (paired-end read)                  strand : +
    // ======================================================= (ADN)
    //                                 -------|\\\\\|--------- strand : +
    //                                        ^ (P2)
      if (P2.getStrand() != pairePos->getStrand() &&
          strcmp(P2.getChrPosition(), pairePos->getChrPosition()) == 0 &&
          P2.getRelativePosition() >= pairePos->getRelativePosition() &&
          (P2.getRelativePosition() - pairePos->getRelativePosition()) <= parameters->max_splice_length) {
        isAValidChimera = true;
      }
    }
  } else {
    // no assumption can be made about the paired-end tag, the chimera is not discarded
    isAValidChimera = true;
  }

  return isAValidChimera;
}
