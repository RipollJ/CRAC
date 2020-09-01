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

#include "classifyTags.h"
#include "utils.h"
#include "libReadsInfo/tagInfo.h"
#include "Support.h"
#include "SupportBreak.h"
#include <config.h>
#include <pthread.h>
#include <cstdlib>
#include <string>
#include <sstream>

#ifdef HAVE_LIBPROGRESSBAR
#include <sys/ioctl.h>
#include <libProgressBar/progressBar.h>
#endif

using namespace std;

ClassifyTags::ClassifyTags(LocateOnGenome *genome,
			   gkArrays *tags,
			   Parameters *param):
  genome(genome),tags(tags),parameters(param)
{}

void ClassifyTags::classify(ostream *snp
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
			    , ostream *single
			    , ostream *sam
          , ostream *pairedEndChimera) {
  Classifier *classifier;
  pthread_t threads[parameters->nb_threads];
  classify_thread_t params[parameters->nb_threads];
  readIterator *readIt = tags->getReads()->begin();
  Classifier **store_tagsInfo[parameters->nb_threads-1];
  sem_t remaining_seats[parameters->nb_threads-1];
  sem_t processed_tags[parameters->nb_threads-1];
  uint *positions = NULL;
  uint nb_tagsInfo;
  bool isPairedEnd = tags->getReads()->isPairedEnd();


#ifdef HAVE_LIBPROGRESSBAR
  DoccY::ProgressBar *PB = NULL;
  unsigned int maxval = 0;
  unsigned int cpt = 0;
  unsigned int reset_cpt = 0;
  if (parameters->show_progressbar) {
    unsigned int cols = 80;
#ifdef TIOCGSIZE
      struct ttysize ts;
      cols = ioctl(STDIN_FILENO, TIOCGSIZE, &ts) ? 80 : ts.ts_cols;    
#elif defined(TIOCGWINSZ)
      struct winsize ts;
      cols = ioctl(STDIN_FILENO, TIOCGWINSZ, &ts) ? 80 : ts.ws_col;
#endif /* TIOCGSIZE */

    maxval = cols ? cols : 100;
    reset_cpt = tags->getNbTags() / maxval;
    PB = new DoccY::ProgressBar("Processing Reads", maxval, cols, cerr, false);
    PB->ShowPercent();
    PB->ShowTime();
    PB->update();
  }
#endif

  for (uint i=0; i <= NB_MASKS; i++)
    nb_classes[i] = 0;

  nb_explainable=0;

  if (parameters->nb_threads >  1) {
    nb_tagsInfo = parameters->nb_tags_info_stored;

    positions = new uint[parameters->nb_threads - 1];

    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      positions[j] = 0;
      params[j].classify = this;
      sem_init(&remaining_seats[j],0,nb_tagsInfo-1);
      sem_init(&processed_tags[j],0,0);
      params[j].remaining_seats = &remaining_seats[j];
      params[j].processed_tags = &processed_tags[j];
      params[j].total_seats = nb_tagsInfo;
      if(isPairedEnd) {
        store_tagsInfo[j] = (Classifier**)calloc(nb_tagsInfo, sizeof(PairedEndReadClassifier*));
      } else {
        store_tagsInfo[j] = (Classifier**)calloc(nb_tagsInfo, sizeof(SingleReadClassifier*));
      }
      params[j].seats = store_tagsInfo[j];
      params[j].thread_num = j;
      params[j].nb_threads = parameters->nb_threads-1;
    }
  } else {
    nb_tagsInfo = 1;
  }

  if (parameters->nb_threads > 1) {
    // Launching threads
    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      pthread_create(&threads[j], NULL, fillInfos, &params[j]);
    }
  } 

  uint i = 0;
  while (!readIt->isFinished()) {
    if (parameters->nb_threads > 1) {
      uint j;
      if(isPairedEnd && i > 0) {
        j = i/2 % (parameters->nb_threads - 1);
      }
      else {
        j = i % (parameters->nb_threads - 1);
      }
      sem_wait(&processed_tags[j]);
      classifier = store_tagsInfo[j][positions[j]];
      sem_post(&remaining_seats[j]);
      positions[j] = (positions[j] + 1) % nb_tagsInfo;
    } else {
      if(isPairedEnd) {
        classifier = new PairedEndReadClassifier(i,i+1,genome,tags,parameters);
      } else {
        classifier = new SingleReadClassifier(i,genome,tags,parameters);
      }
      classifier->classify();
    }

    classifier->setInfos(readIt);
    classifier->writeOutputs(snp,
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
    classifier->samOutput(*sam);
    if(isPairedEnd && pairedEndChimera != NULL)
      ((PairedEndReadClassifier*)classifier)->writePairedEndChimera(pairedEndChimera);
    classifier->updateStatistics(&nb_classes,&nb_explainable);
    delete classifier;

    // If we are using paired-end tags, increments loop ind by 2
    if(isPairedEnd) {
      i += 2;
#ifdef HAVE_LIBPROGRESSBAR
      cpt += 2;
#endif
    } else {
      i++;
#ifdef HAVE_LIBPROGRESSBAR
      cpt++;
#endif
    }
#ifdef HAVE_LIBPROGRESSBAR
    if ((cpt >= reset_cpt) && parameters->show_progressbar) {
      PB->Step();
      cpt = 0;
    }
#endif
  }
#ifdef HAVE_LIBPROGRESSBAR
  if (parameters->show_progressbar) {
    PB->SetVal(maxval);
    PB->update(false);
    cerr << endl;
    delete PB;
  }
#endif

  // deleting the iterator
  delete readIt;
  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    free(store_tagsInfo[j]);
    sem_destroy(&remaining_seats[j]);
  }
  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    pthread_join(threads[j], NULL);
  }

  if (positions)
    delete [] positions;
}

LocateOnGenome *ClassifyTags::getGenomeIndex() {
  return genome;
}

Parameters *ClassifyTags::getParameters() {
  return parameters;
}

gkArrays *ClassifyTags::getReadIndex() {
  return tags;
}

uint ClassifyTags::getNbNone() {
  return nb_classes[posBitInMask(MASK_NONE)];
}

uint ClassifyTags::getNbMultiple() {
  return nb_classes[posBitInMask(MASK_MULTIPLE)];
}

uint ClassifyTags::getNbSingle() {
  return nb_classes[posBitInMask(MASK_SINGLE)];
}

uint ClassifyTags::getNbAlmostNormal() {
  return nb_classes[posBitInMask(MASK_ALMOST_NORMAL)];
}

uint ClassifyTags::getNbExplainable() {
  return nb_explainable;
}

uint ClassifyTags::getNbNormal() {
  return nb_classes[posBitInMask(MASK_NORMAL)];
}

uint ClassifyTags::getNbDuplication() {
  return nb_classes[posBitInMask(MASK_DUPLICATION)];
}

uint ClassifyTags::getNbSNP() {
  return nb_classes[posBitInMask(MASK_SNP)];
}

uint ClassifyTags::getNbBioTagIndel() {
  return nb_classes[posBitInMask(MASK_BIOLOGICAL_TAG_INDEL)];
}

uint ClassifyTags::getNbBioUndetermined() {
  return nb_classes[posBitInMask(MASK_BIOLOGICAL_UNDETERMINED)];
}

uint ClassifyTags::getNbSeqErr() {
  return nb_classes[posBitInMask(MASK_SEQ_ERR)];
}


uint ClassifyTags::getNbRepetition() {
  return nb_classes[posBitInMask(MASK_REPETITION)];
}

uint ClassifyTags::getNbSplice() {
  return nb_classes[posBitInMask(MASK_SPLICE)];
}

uint ClassifyTags::getNbSpliceIntra() {
  return nb_classes[posBitInMask(MASK_INTRA_TRANSPLICING)];
}

uint ClassifyTags::getNbSpliceInter() {
  return nb_classes[posBitInMask(MASK_INTER_TRANSPLICING)];
}

uint ClassifyTags::getNbSpliceNoCover() {
  return nb_classes[posBitInMask(MASK_SPLICE_NO_COVER)];
}

uint ClassifyTags::getNbNothing() {
  return nb_classes[NB_MASKS];
}

uint ClassifyTags::getNbUndetermined() {
  return nb_classes[posBitInMask(MASK_UNDETERMINED_ERROR)];
}

uint ClassifyTags::getNbPairedEndChimera() {
  return nb_classes[posBitInMask(MASK_PAIRED_END_CHIMERA)];
}

void ClassifyTags::getHeaders(ostream *snp
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
			    , ostream *single
          , ostream *pairedEndChimera) {

  stringstream headerTemplate;
  headerTemplate << "#CRAC version " << PACKAGE_VERSION << endl;
  headerTemplate << "#Created on " <<__DATE__ << endl << "#" << endl;

  string chrPosPatternInfo = "#A pattern [ chr|strand,pos ] is an occurrence of location on the index (genome or reference)";
  
  // For SNV
  if (snp){
    *snp << headerTemplate.str() << chrPosPatternInfo << endl;
    *snp << "#A tag \"single\" is noted if there is one possibility of SNV, a tag \"duplicate\" is noted otherwise" << endl;
    *snp << "#read_id tag_snv loc_snv_on_genome pos_snv_on_read snv crac_score single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For indels
  if (bioTagIndel){
    *bioTagIndel << headerTemplate.str() << chrPosPatternInfo << endl;
    *bioTagIndel << "#A tag \"single\" is noted if there is one possibility of indel, a tag \"duplicate\" is noted otherwise" << endl;
    *bioTagIndel << "#read_id tag_indel loc_indel_on_genome pos_start_indel_on_read nb_of_insertion nb_of_deletion single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For error
  if (seqErr){
    *seqErr << headerTemplate.str() << chrPosPatternInfo << endl;
    *seqErr << "#A tag \"single\" is noted if there is one possibility of error, a tag \"duplicate\" is noted otherwise" << endl;
    *seqErr << "#Two different formats for error: "<< endl;
    *seqErr << "#read_id tag_error loc_error_on_genome pos_ponctual_error_on_read error crac_score single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
    *seqErr << "#or" << endl;
    *seqErr << "#read_id tag_error loc_error_on_genome pos_start_indel_error_on_read nb_of_insertion nb_of_deletion single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For splice and splice_no_cover
  if (splice){
    *splice << headerTemplate.str() << chrPosPatternInfo << endl;
    *splice << "#A tag \"single\" is noted if there is one possibility of splice, a tag \"duplicate\" is noted otherwise" << endl;
    *splice << "#read_id tag_splice loc_end_first_exon_on_genome pos_junction_on_read splice_length single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }

  if (spliceNoCover) {
    *spliceNoCover << headerTemplate.str() << chrPosPatternInfo << endl;
    *spliceNoCover << "#A tag \"single\" is noted if there is one possibility of splice, a tag \"duplicate\" is noted otherwise" << endl;
    *spliceNoCover << "#read_id tag_splice loc_end_first_exon_on_genome pos_junction_on_read splice_length single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For chimera
  if (spliceInter){
    *spliceInter << headerTemplate.str() << chrPosPatternInfo << endl;
    *spliceInter << "#A tag \"single\" is noted if there is one possibility of chimera, a tag \"duplicate\" is noted otherwise" << endl;
    *spliceInter << "#read_id tag_chimera loc_end_first_exon_on_genome loc_start_second_exon_on_genome pos_junction_on_read single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }

  if (pairedEndChimera){
    *pairedEndChimera << headerTemplate.str() << chrPosPatternInfo << endl;
    *pairedEndChimera << "#A paire of paired-end tags classified as \"single\" are noted if there is one possibility of chimera between them" << endl;
    *pairedEndChimera << "#chimera_class read_id1 single_loc_on_genome1(private) pos_single_loc_on_read1(private) read1 p_support1 p_loc1 read_id2 single_loc_on_genome2(private) pos_single_loc_on_read2(private) read2 p_support2 p_loc2" << endl;
  }

  // For undetermined cause and bioUndermined
  if (undetermined){
    *undetermined << headerTemplate.str() << chrPosPatternInfo << endl;
    *undetermined << "#read_id undetermined_cause_features single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  if (bioUndermined) {
    *bioUndermined << headerTemplate.str() << chrPosPatternInfo << endl;
    *bioUndermined << "#read_id bioUndetermined_cause_features single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For repetition
  if (repetition){
    *repetition << headerTemplate.str() << chrPosPatternInfo << endl;
    *repetition << "#read_id pos_start_repeat_on_read pos_end_repeat_on_read single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For duplication and multiple
  if (duplication){
    *duplication << headerTemplate.str() << chrPosPatternInfo << endl;
    *duplication << "#read_id occurrence_loc_on_genome(private) pos_loc_on_read(private) read p_support p_loc" << endl;
  }
  if (multiple) {
    *multiple << headerTemplate.str() << chrPosPatternInfo << endl;
    *multiple << "#read_id one_occurrence_loc_on_genome(private) pos_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For nothing
  if (nothing){
    *nothing << headerTemplate.str();
    *nothing << "#read_id read p_support p_loc" << endl;
  }
  // For normal and almostNormal
  if (normal){
    *normal << headerTemplate.str() << chrPosPatternInfo << endl;
    *normal << "#read_id single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  if (almostNormal) {
    *almostNormal << headerTemplate.str() << chrPosPatternInfo << endl;
    *almostNormal << "#read_id single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
  // For none
  if (none){
    *none << headerTemplate.str() << chrPosPatternInfo << endl;
    *none << "#read_id [loc_on_genome(private) pos_loc_on_read(private)] read p_support p_loc" << endl;
  }
  // For single
  if (single){
    *single << headerTemplate.str() << chrPosPatternInfo << endl;
    *single << "#read_id single_loc_on_genome(private) pos_single_loc_on_read(private) read p_support p_loc" << endl;
  }
}

/// PRIVATE ///

void *fillInfos(void *args) {
  classify_thread_t params = ((classify_thread_t *)args)[0];
  uint id = params.thread_num;
  uintSA nb_reads = params.classify->getReadIndex()->getNbTags();
  bool isPairedEnd = params.classify->getReadIndex()->getReads()->isPairedEnd();
  uint current_pos = 0;
  uint k = id;
  uint increment = params.nb_threads;
  if(isPairedEnd) {
    k = 2*k;
    increment = 2*increment;
  }
  for (; k < nb_reads; k += increment) {
    // We need one seat for the current read
    if (isPairedEnd) {
      params.seats[current_pos] = new PairedEndReadClassifier(k,k+1,params.classify->getGenomeIndex(),params.classify->getReadIndex(),params.classify->getParameters());
    } else {
      params.seats[current_pos] = new SingleReadClassifier(k,params.classify->getGenomeIndex(),params.classify->getReadIndex(),params.classify->getParameters());
    }
    params.seats[current_pos]->classify();
    //params.seats[current_pos]->generateSAM();
    sem_post(params.processed_tags);
    sem_wait(params.remaining_seats);
    current_pos = (current_pos + 1) % params.total_seats;
  }
  pthread_exit(NULL);
  return NULL;
}
