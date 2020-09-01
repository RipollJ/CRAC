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

#ifndef TAG_INFO_H
#define TAG_INFO_H
#include <iostream>
#include <stdarg.h>
#include "SNPInfo.h"
#include "BioTagIndelInfo.h"
#include "BioUndeterminedInfo.h"
#include "DuplicationInfo.h"
#include "SpliceInfo.h"
#include "SpliceInterInfo.h"
#include "SpliceNoCoverInfo.h"
#include "SpliceIntraInfo.h"
#include "SeqErrInfo.h"
#include "repetitionInfo.h"
#include "undeterminedErrorInfo.h"
#include "genericInfo.h"
#include "../libSSA/chrPosition.h"
#include "../libSSA/locateOnGenome.h"
#include "../types.h"
#include "../Support.h"
#include "../SupportBreak.h"
#include "../Parameters.h"
#include "samLine.h"

using namespace std;

class TagInfo {
 private:
  int code;
  LocateOnGenome *genome;
  byte *nb_each_type;
  GenericInfo ***info_each_type; /* Stores the info on breaks depending on 
                                  * the type it corresponds to */
  uint nb_causes;        /* Number of causes found in the read 
                          * INFO_DUPLICATION and INFO_REPETITION are not stored
                          * are not taken into account
                          */
  byte *types_in_order; /* Stores in which order the add* function have been 
                           called.
                           Each cell stores a type (INFO_SNP, INFO_SEQ_ERR, ...).
                           INFO_DUPLICATION and INFO_REPETITION are not stored
                           This array has nb_causes cells.
                        */
  byte *nb_causes_per_break;    /* Number of causes in each break
                                 * INFO_DUPLICATION and INFO_REPETITION are not stored
                                 * are not taken into account
                                 * This array has at least 
                                 * getNbBreaks() cells */

  Support *support;
  uint current_break;		/* The current break id processed */

  ChrPosition **chrpos[NB_POS_BREAK];

  void addGenericElement(byte type, GenericInfo *gi);

  byte getNbGenericElement(uint type);
  GenericInfo **getGenericInfos(uint type);

  char *read_name, *read_quality;
  uint tagId;


  
 public:
  /**
   * The support is not copied but is deleted when this is deleted
   */
  TagInfo(LocateOnGenome *genome, Support *s);

  ~TagInfo();

  float getAverageHighInside(uint);
  float getAverageLowInside(uint);
  SupportBreak *getBreak(uint, int = 1);
  /**
   * @return the type of the i-th cause 
   */
  byte getCauseType(uint i, int=1);
  ChrPosition *getChrPosEndBreak(uint);
  ChrPosition *getChrPosStartBreak(uint);
  uint getCurrentBreak();
  float getDeviationInsideBreak(uint);
  int getEndPosRepeat();
  gap_size_t getGenomeGapLength(uint);
/*   bool getIsExtended(uint); */
  uint *getLocalisations();
  ChrPosition *getLocation();
  uint getLocationEndBreak(uint);
  uint getLocationStartBreak(uint);
  uint getNbBreaks();
  uint getNbCauses();
  byte getNbCausesInBreak(uint, int = 1);
  uint getNbDuplicate();
  uint getNbMultiple();
  uint getNbSingle();
  uint getPositionOfLocation();
  char *getReadName();
  char *getReadQuality();
  int getStartPosRepeat();
  uint *getSupport();
  Support *getSupportObject();
  uint getSupportLength();
  /* uint getPositionEndBreak(uint); */
  uint getPositionEndBreak(uint, int=1);
  /* uint getPositionStartBreak(uint); */
  uint getPositionStartBreak(uint, int=1);
  float getScoreComputedIntraExon(uint);
  float getScoreComputedInterExon(uint);
  float getScoreInsideBreak(uint);
  float getScoreInsideAverages(uint);
  float getScoreOutsideBreak(uint);
  char *getTag();
  char *getTagChromosome();
  uint getTagBreakLength(uint);
  ulong getTagPosition();
  bool getTagStrand();
  uint getThreshold();

  bool isAlmostNormal();
  bool isDuplicated(uint);
  bool isExplainable();
  bool isNormal();
  bool isMultiple();
  bool isNone();
  bool isRepeated(uint);
  bool isSingle();
  bool isSupportFallingLeft(uint, int=1);
  bool isSupportFallingRight(uint, int=1);
  bool isDuplication();
  bool hasSNP();
  bool hasBioTagIndel();
  bool hasBioUndetermined();
  bool hasSeqErr();
  bool hasRepetition();
  bool hasSplice();
  bool hasSpliceIntraChr();
  bool hasSpliceInterChr();
  bool hasSpliceNoCover();
  bool hasNothing();
  bool hasUndeterminedError();
  int getCode();

  DuplicationInfo **getInfosDuplication();

  RepetitionInfo **getInfosRepetition();
  byte getNbRepetition();
  
  SNPInfo **getInfosSNP();
  byte getNbSNP();
  
  BioIndelInfo **getInfosBioTagIndel();
  byte getNbBioTagIndel();

  BioUndeterminedInfo **getInfosBioUndetermined();
  byte getNbBioUndetermined();

  SeqErrInfo **getInfosSeqErr();
  byte getNbSeqErr();

  UndeterminedErrorInfo **getInfosUndeterminedError();
  byte getNbUndeterminedError();

  SpliceInfo **getInfosSplice();
  byte getNbSplice();

  SpliceIntraInfo **getInfosSpliceIntra();
  byte getNbSpliceIntra();

  SpliceInterInfo **getInfosSpliceInter();
  byte getNbSpliceInter();

  SpliceNoCoverInfo **getInfosSpliceNoCover();
  byte getNbSpliceNoCover();

  /// COMMANDS ///

  void addDuplication();

  void addAlmostNormal();
  void addNormal();
  void addMultiple();
  void addNone();

  /**
   * position_in: position in the tag where the repetition starts (default: -1
   *              unkown)
   * position_out: position in the tag where the repetition ends.
   */
  void addRepetition(int position_in=-1, int position_out=-1);

  void addSingle();

  /**
   * tag_nuc: nucleotide found in the tag, or 0 if unknown (or deletion)
   * genome_nuc: nucleotide found in the genome, or 0 if unknown (or insertion)
   * snpID: identify which kind of SNP (depending on in which context it occurs)
   */
  void addSNP(char tag_nuc=0,
              error_context snpID=FIRST_SUBSTITUTION);

  /**
   * Another version of addSNP, less easy to use but that allows to give all
   * the necessary informations without letting the other addSNP version
   * computing them.
   * @param pos : position of the SNP in the read
   * @param base_read : nucleotide  from the read (or 0 if it doesn't exist)
   * @param base_genome: corresponding nucleotide on the genome 
   *                     (or 0 if it doesn't exist) 
   * @param chrPos : location of the SNP on the genome
   * @param score: score of the SNP.
   */
  void addSNP(uint pos, char base_read, char base_genome, ChrPosition &chrPos,
              float score);

  /**
   * nb_ins: number of insertions
   * nb_del: number of deletions
   */
  void addBioTagIndel(uint nb_ins, uint nb_del);

  /**
   * We know that there is a biological cause, but since
   * it is at the beginning or at the end of the support, it is impossible
   * to determined accurately the biological cause.
   */
  void addBioUndetermined(ChrPosition *, uint, const char *, ...);
  
  /**
   * position: position in the tag where the sequence error is
   * nb_ins: number of insertions
   * nb_del: number of deletions
   */
  void addSeqErr(error_context , uint position, uint nb_ins, uint nb_del,
		 float score,
		 uint genome_seq_length=~0,
		 char *tag_seq=NULL, uint tag_seq_length=~0);

  
  /**
   * position: position where the splice is
   * gap_length: length of the intron
   */
  void addSplice(uint gap_length);
  void addSpliceIntra(uint gap_length);
  void addSpliceInter();
  void addSpliceNoCover(uint gap_length);
  
  void addUndeterminedError(const char *format, ...);

  /**
   * Change a GenericInfo with an other one.
   * The new GenericInfo will be inserted at the same place as the old one,
   * Indeed, it will apply to the same break.
   * This function is useful when we are reclassifying breaks afterwards.
   *
   * @param type the type of the old genericInfo to replace
   * @param newType the new type that will replace the current type
   * @param oldMask: old mask corresponding to the old cause to be replaced
   *                  (see constant whose name starts with MASK_ in const.h)
   * @param newMask: new mask corresponding to the new cause replacing the new one
   * @param i the index of the GenericInfo in info_each_type[type] array.
   * @param newGi the new GenericElement to insert in array info_each_type[newType]
   *
   * @post getNbGenericElement(type) == old getNbGenericElement(type) - 1
   *    && getNbGenericElement(newType) == old getNbGenericElement(newType) + 1
   *    && (old getNbGenericElement(type) == 1 ==> getCode() & type == 0)
   *    && getCode() & newType == newType
   */
  void changeGenericElement(byte type, byte newType, int oldMask, int newMask, 
                            uint i, GenericInfo *newGi);

  /**
   * Same function as above.
   * The difference resides in the parameter given.
   * On that version, you must provide the break ID, the cause lies in, and additionally
   * the number of a cause of that type in this break (starts at 0).
   *
   * @param type the type of the old genericInfo to replace
   * @param newType the new type that will replace the current type
   * @param oldMask: old mask corresponding to the old cause to be replaced
   *                  (see constant whose name starts with MASK_ in const.h)
   * @param newMask: new mask corresponding to the new cause replacing the new one
   * @param break_id The ID of the break where we must change something
   * @param cause_num The number of the cause (starts at 0) we want to change, having
   *                  the same type, in break break_id.
   * @param newGi the new GenericElement to insert in array info_each_type[newType]
   */
  void changeGenericElement(byte type, byte newType, int oldMask, int newMask, 
                            uint break_id, uint cause_num, GenericInfo *newGi);


  /**
   * Remove a splice Inter and add a bioUndetermined info instead.
   *
   * @param i the indice of the spliceInter if Info_each_type array
   * @param info a string containing the reason of the reclassification
   */
  void removeSpliceInter(uint i,const char* info);

  /**
   * SAM output for the given read
   * XU: Is it single located?
   * XD: It it duplicated?
   * XM: Is it multiple?
   * XN: Is it normal? (1: normal, 2: almost, 0: no)
   * XL: location taken to minimise the probability of FP.
   * XP: position in the read of XL
   * XO: Location taken in Support as a single location
   * XQ: position in the read of XO
   * XC: Number of causes in the read
   * XE: Detail of a cause
   * XB: Large details for a cause
   */
  ostream &samOutput(ostream &os);

  /*
   * Return a set of sam line that corresponds to the alignement of the read.
   * Multiple line for one read can occur in two situations:
   * 1. the read has multiple locations on the genome
   * 2. the read has a chimeric alignement
   * The primary line of the read is always returned if first postion
   * of the vector. 
   */
  vector<SamLine*> *getSamLines();

  /**
   * Tells which break we are working on.
   * The given integer is the break number.
   * @post Calling one of addSNP, addGenomeIndel, ... will affect one of those causes
   * to the current break
   */
  void setCurrentBreak(uint);

  void setReadName(char *name);
  void setReadQuality(char *);
  void setReadNumber(uint tagId);
  uint getReadNumber(); 

  friend ostream &operator<<(ostream &os, TagInfo &);

 private:
  pair<ChrPosition, uchar *> retrieveChrPosAndNuc(error_context type,
                                                  uint nb_nuc,
                                                  int &shift);


  /**
   * This method constructs from a single SamLine object a set of SamLines
   * that corresponds to the chimericAlignement of the read if there is one.
   * This method is recursive and called by the "getSamLine()" method of this
   * very Class.
   */
  vector<SamLine*> *getChimericAlignements(SamLine *base_line, uint cigar_index, uint splice_id, uint left_softclip, uint *right_softclip);

  /**
   * Return a prototype of SAMLine that can be used to create
   * multiple alignement for one tagInfo.
   * @loc is the location of a chosen kmer in the read which represents
   * the read location
   * @pos_loc the kmer position in the read
   */
  SamLine *getSamLine(ChrPosition *loc = NULL, ulong pos_loc = 0);


  
};

ostream &operator<<(ostream &os, TagInfo &);


#endif
