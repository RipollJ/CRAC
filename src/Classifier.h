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

#ifndef CLASSIFIER_H
#define CLASSIFIER_H

#include <iostream>
#include <gkArrays.h>
#include "libReadsInfo/tagInfo.h"
#include "Support.h"
#include "const.h"
#include "SupportBreak.h"
#include <config.h>

class Classifier {
protected:
  string output_sam;

public:
  virtual ~Classifier(){}

  /*
   * Classify read(s) depending of the implementation 
   * of Classifier
   */
  virtual void classify()=0;

  /*
   * Set Informations about read(s) using the readIteror
   */
  virtual readIterator *setInfos(readIterator *readIt)=0;

  /**
   * Write the line(s) in SAM file
   */
  virtual ostream &samOutput(ostream &os)=0;


  /**
   * Print specific information about classification in
   * each file
   */
  virtual void writeOutputs(ostream *snp
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
			    , ostream *single)=0;

  /**
   * Update statitics about classification
   */
  virtual void updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable)=0;

};

class SingleReadClassifier : public Classifier {
private:
  Support *suppo;
  TagInfo *taginfo;
  gkArrays *tags;
  LocateOnGenome *genome;
  Parameters *parameters;
  uint tagId;

public:
  /**
   * Constructor of SingleReadClassifier
   *
   * @param tagId the id of the tag in the gkArrays
   * @param genome the LocateOnGenome object
   * @param tags the gkArrays object
   * @param parameters the parameters object
   */
  SingleReadClassifier(uint tagId, LocateOnGenome *genome, gkArrays *tags, Parameters *parameters);
  ~SingleReadClassifier();

  /**
   * Classify the read contained if the classifier
   */
  virtual void classify();

  /**
   * Set Informations about the read using the iterator.
   * This method also move the iterator one step further.
   */
  virtual readIterator *setInfos(readIterator *readIt);

  /**
   * Print the SAM line of the read in the file in arguement
   *
   * @param os the stream where the SAM line will be written
   */
  virtual ostream &samOutput(ostream &os);

  /**
   * Print informations about the classification of the read
   * in the related stream
   */
  virtual void writeOutputs(ostream *snp
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
			    , ostream *single);

  /**
   * Update the statistics about classification with the current classified read
   *
   * @param nb_classes the array containing statistics about each type of classification
   * @param nb_explainable the number of read explainable
   */
  virtual void updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable);

  /**
   * @return the tag id of the read handled by the classifier
   */
  uint getTagId();

  /**
   * @return the tagInfo that has been created by the classifier.
   */
  TagInfo *getTagInfo();

private:
  /**
   * Return a substring of the tag of length <length> and which
   * ends at the position corresponding to the end of the break.
   * The position retrieved can be shifted using value <shift>.
   * length > params->max_bases_retrieved => result == NULL
   */
 char *getTagAtEndBreak(SupportBreak *suppoBreak, char *tag,
			uint length, uint shift=0);

  /**
   * Try to find a SNP even is some information is lacking (eg. the break
   * is not complete).
   * This is done by naively retrieving the modified nucleotide on the genome
   * and to check if we succeed to locate (a) k-mer(s) on the genome when
   * correcting the faulty nucleotide.
   * @post if a (or several) SNP(s) is/are found, the TagInfo is updated
   * accordingly.
   */
  void perform_deep_snp_search(SupportBreak *, uint break_num);

  /**
   * Tests if the supportBreak that have been classified as a chimera
   * is a stringent chimera according to stringent classification tests.
   * This function stops when one of the test has failed and return a 
   * specific flag for this cause.
   * This does not mean that the chimera has only failed this tests, but
   * that this test is the first one that has failed.
   *
   * @return 0 if the chimera is stringent
   *         1 if the chimera has a break_length too small
   *         2 if the chimera comes from a break with too many merges
   *         3 if the chimera is in a repeated zone
   *         4 if the chimera has a high support variation
   *         5 if the kmer verification test doesn't match the chimera
   *         6 if no single loc have been found around the break
   */
  uint isStringentChimera(SupportBreak *, uint break_num);
};

class PairedEndReadClassifier : public Classifier {
private:
  gkArrays *tags;
  LocateOnGenome *genome;
  Parameters *parameters;
  SingleReadClassifier *classifier1;
  SingleReadClassifier *classifier2;
  // 0 : nothing
  // 1 : Class 1 chimera
  // 2 : Class 2 chimera
  // 3 : Class 3 chimera
  // 4 : Class 4 chimera
  uint paired_end_chimera;

public:
  /**
   * Constructor of PairedEndReadClassifier
   *
   * @param tagId1 the id of the first tag of the pair in the gkArrays
   * @param tagId2 the id of the second tag of the pair in the gkArrays
   * @param genome the LocateOnGenome object
   * @param tags the gkArrays object
   * @param parameters the parameters object
   */
  PairedEndReadClassifier(uint tagId1, uint tagId2, LocateOnGenome *genome, gkArrays *tags, Parameters *parameters);
  ~PairedEndReadClassifier();

  /**
   * Classify paired-end reads contained in the classifier
   */
  virtual void classify();

  /**
   * Set Informations about the paird-end reads using the iterator.
   * This method also move the iterator two step further.
   */
  virtual readIterator *setInfos(readIterator *readIt);

  /**
   * Print two SAM lines corresponding to the paired-end reads classified by the
   * Classifier in the output file in parameter
   *
   * @param os the stream where the SAM lines will be written
   */
  virtual ostream &samOutput(ostream &os);

  /**
   * Print informations about the classification of paired-end reads
   * in the related stream
   */
  virtual void writeOutputs(ostream *snp
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
			    , ostream *single);
  /**
   * Update the statistics about classification with the information about paired-end reads classified
   *
   * @param nb_classes the array containing statistics about each type of classification
   * @param nb_explainable the number of read explainable
   */
  virtual void updateStatistics(uint (*nb_classes)[NB_MASKS+1], uint *nb_explainable);

  /**
   * If the paired-end tags contains a chimera in the non sequenced region
   * between the tags informations about it in the pairedEndChimera file.
   * Otherwise nothing is done.
   */
  void writePairedEndChimera(ostream *pairedEndChimera);

  /**
   * Return true if the paired-end tags contains a chimera in the non sequenced
   * region between the tags.
   */
  bool hasPairedEndChimera();

  /**
   * @return the tagInfo that has been created by the classifier.
   */
  TagInfo *getFirstTagInfo();
  TagInfo *getSecondTagInfo();

private:
  /**
   * Check if a chimera is valid according to the informations
   * about its paired read location.
   *
   * @param chimera the chimera to analyze
   * @param pairedEndClassifier the paired read classfier
   * @return true if chimera is valid, false if not.
   */
  bool chimeraPairedEndCheck(SpliceInterInfo *chimera, SingleReadClassifier *pairedEndClassifier);

  /*
   * Set @line optional fields deticated to store mate informations
   */
  void setPairedEndOptionalFields(SamLine &line, const SamLine &paired_line);

  /**
   * Write a vector of samlines to the output stream.  According to the primary
   * line and the paired-end read primary line the method update some flags to
   * be SAM consistent.
   */
  void writeSamLines(ostream &os, vector<SamLine*> &sam_lines, SamLine &primary_line, SamLine &paired_primary_line, bool is_first_taginfo);
};
#endif
