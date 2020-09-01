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

#include <iostream>
#include <fstream>
#include <ostream>
#include <string>
#include <cstdlib>
#include <algorithm>
#include "types.h"
#include "utils.h"
#include <gkArrays.h>
#include <gkArraysTypes.h>
#include "classifyTags.h"
#include "libSSA/locateOnGenome.h"
#include "locateTags.h"
#include "libSSA/cracIndex.h"
#include "Parameters.h"
#include "locateServer.h"
#include <cstdio>
#include <unistd.h>
#include <getopt.h>
#include <sys/stat.h>
#ifdef HAVE_CONFIG_H
#  include <config.h>
#endif
#include <gzstream.h>

using namespace std;

#define SNP 8
#define BIOTAG 9
#define SEQERR 10
#define SPLICE 11
#define SPLICE_NOCOVER 12
#define CHIMERA 13
#define UNDETERMINED 14
#define REPETITION 15
#define DUPLICATION 16
#define NOTHING 17
#define NORMAL 18
#define ALMOSTNORMAL 19
#define MULTIPLE 20
#define NONE 21
#define BIOUNDETERMINED 22
#define SINGLE 23

#define MAX_SPLICE 24
#define MAX_BIO_INDEL 25
// #define MAX_SEQERR_INDEL 26 
#define MAX_DUPLICATION 27 /* Maximal number of occurrences for a factor
					  to be considered as duplicated */
#define MIN_DUPLICATION 28 /* Minimal number of occurrences for a factor
					  to be considered as duplicated */
#define PERCENT_MIN_SINGLE_REPETITION 29 /* For classifying as repetition, minimal percentage of factors located once */
#define PERCENT_MIN_MULTIPLE_CST 30 /* To be classified as single, a support
                                  * must have at least PERCENT_MIN_SINGLE%
                                  * single-located factors 
				  */
#define PERCENT_MIN_DUPLICATION 31  /* To be classified as duplicate, a support must 
					 have at least this percentage of factors 
					 located at most MAX_LOCALISATION_DUPLICATION on the genome */
#define MIN_LOC_REPETITION 32	/* Minimal number of occurrences of a factor, for being a repetition */
#define PERCENT_VARIATION_COVERAGE 33 /* 50% */
#define P_VALUE_VARIATION_COVERAGE 34 /* P-Value which is the relation of two average
					  * to determine whether we have a sequence error 
					  * or a biological explanation. */
#define MAX_BIOLOGICAL_DEVIATION 35 /* Maximal deviation allowed for the support
					  * inside the break, to consider it as a biological event*/

#define MAX_RANDOMLY_MATCHED 36 /* Maximal number of bases we may match, before the
				       * expected position (e.g. for genome ins/del,
				       * there is 0.25 probability to match a factor of
				       * the tag one position, before the one expected)*/

#define MAX_EXTENSION 38 /* A break is extended right and left by, at most,
				 * MAX_EXTENSION_LENGTH bases. */
#define MIN_BEFORE_BREAK 39 /* Minimal number of bases we must have to consider a
				    new break */ 
#define MAX_RETRIEVED 40	/* Maximal number of bases retrieve for information
				 *  when discovering  a sequence error. */

#define SUPPORT_SCORE_WINDOW 41

#define MINIMUM_SUPPORT_NO_COVER 42

#define MINIMUM_BREAK_LENGTH 43

#define MAXIMUM_SUPPORT_OUT_NO_COVER 44

#define MAX_AMBIGUOUS_AVG_HIGH 45
#define NB_BASES_SUBSTRACTED_LOW_AVG 46
#define NB_POS_CHECK_ONES 47
#define MIN_PERC_ONES_INSIDE 48
#define MIN_RATIO_SUPPORT_FALL_CST 49
#define NB_TAGS_INFO_CST 50
#define NB_THREADS_CST 51

#define ALL_OUTPUT 52

#define MAX_NB_OVERLAPPING_BREAKS_CST 53
#define PERCENT_MIN_SINGLE_CST 54 /* To be classified as single, a support
                                  * must have at least PERCENT_MIN_SINGLE%
                                  * single-located factors 
				  */
#define SCORE_INTRA_AMBIGUOUS_CST 55
#define SCORE_INTER_AMBIGUOUS_CST 56

#define SAM 57

#define DEEP_SNP_SEARCH_CST 58
#define NUMBER_NUCLEOTIDES_SNP_COMPARISON_CST 59
#define DETAILED_SAM_CST 60
#define SERVER_CST 61
#define INPUT_FIFO_CST 62
//63 for '?'
#define OUTPUT_FIFO_CST 64

#define MAX_VERYNICE_MERGE_CST 65
#define MAX_NB_LOCATED_OCCS_CST 66  // -n

#define READS_LENGTH 67
#define IS_STRANDED_CST 68

#define DGE 69
#define STRINGENT_CHIMERA_CST 70
#define NO_AMBIGUITY_CST 71
#define PAIRED_END_CHIMERA 72
#define GZ_OUTPUT 73
#define SUMMARY 74
#define PROGRESSBAR 75
#define TREAT_MULTIPLE_CST 76

void print_help() {
  printf(PACKAGE_STRING);
  printf("\tCompiled on %s.\n\n",__DATE__);

  printf("   -h, --help           <none>          print this help and exit\n");
  printf("   -f, --full-help      <none>          print a complete help and exit\n");
  printf("   -v                   <none>          print version and exit\n\n");
 
  printf("Mandatory arguments\n");  
  printf("   -i                   <FILE>          set genome index file (without the extension filename)\n");
  printf("   -r                   <FILE> [FILE2]  set read file. Specify FILE2 in case of paired-end reads\n");
  printf("   -k                   <INT>           set k-mer length\n");
  printf("   -o, --sam            <FILE>          set SAM output filename or print on STDOUT with \"-o -\" argument\n\n");
  

  printf("Optional arguments\n");
  printf("  * Protocol\n");
  printf("   --stranded           <none>          set the read mapping with for a strand specific library (DEFAULT non-strand specific)\n\n");
  printf("  * Efficiency\n");
  printf("   --nb-threads         <INT>           set the number of worker threads (DEFAULT %d)\n", NB_THREADS);
  printf("   --reads-length, -m    <INT>           set read length in case of all reads have the same length to optimize\n\
                                        CPU and memory times\n"); 
  printf("   --treat-multiple     <none>          write alignments with multiple locations (>max-duplication and <=max-locs) in the SAM file rather than only one occurrence\n"); 
  printf("   --max-locs           <INT>           set the maximum number of locations on the reference index (DEFAULT %d)\n\n", MAX_NB_LOCATED_OCCS);
  printf("  * Accuracy\n");
  printf("   --no-ambiguity       <none>          discard biological events (splice, snv, indel, chimera) which have several matches on the reference index\n\n");
}

void print_fullhelp() {
  print_help();
  printf("\nOptional output arguments\n");
  printf("   --all                              <FILE>     set output base filename for all homemade file formats (see man page for more details)\n\
                                                 in addition to SAM output \n");
  printf("   --gz                               <none>     all output files specified after this argument are gzipped (included SAM output file)\n");
  printf("   --summary                          <FILE>     set output summary file with some statistics about mapping and classification\n"); 
#ifdef HAVE_LIBPROGRESSBAR
  printf("   --show-progressbar                 <none>     show a progress bar for the process times on STDERR\n\n"); 
#else
  printf("\n"); 
#endif

  // printf("  * Mapping\n");
  // printf("   --single                           <FILE>     set output single file\n"); 
  // printf("   --duplicate                        <FILE>     set output duplication file\n"); 
  // printf("   --multiple                         <FILE>     set output multiple file\n"); 
  // printf("   --none                             <FILE>     set output none file\n"); 
  // printf("   --normal                           <FILE>     set output normal file\n"); 
  // printf("   --almost-normal                    <FILE>     set output almost normal file\n\n"); 

  // printf("  * Biological causes\n");
  // printf("   --snv                              <FILE>     set output SNV file\n"); 
  // printf("   --indel                            <FILE>     set output short indel file\n");
  // printf("   --splice                           <FILE>     set output splice junction file\n"); 
  // printf("   --weak-splice                      <FILE>     set output coverless splice junction file\n"); 
  // printf("   --chimera                          <FILE>     set output chimera junction file\n"); 
  // printf("   --paired-end-chimera               <FILE>     set output for paired-end chimera file\n");
  // printf("   --biological                       <FILE>     set output bio-undetermined file\n\n");  

  // printf("  * Sequence errors\n");
  // printf("   --errors                           <FILE>     set output sequence errors file\n\n"); 

  // printf("  * Repetition\n");
  // printf("   --repeat                           <FILE>     set output repetition file\n\n"); 
   
  // printf("  * Other causes\n");
  // printf("   --undetermined                     <FILE>     set output undetermined file\n"); 
  // printf("   --nothing                          <FILE>     set output nothing file\n\n"); 

  printf("Optional process for specific research\n");
  printf("   --deep-snv                         <none>     will search hard to find SNPs\n");
  printf("   --stringent-chimera                <none>     will search chimeras with more accuracy (but less sensitivity)\n\n");

  printf("Optional process launcher (once must be selected)\n");
  printf("  * Exact matching tool\n");
  printf("   --emt                              <none>     launch CRAC-emt for exact mapping of short reads\n\n");
  printf("  * Server tool (for debugging) \n");
  printf("   --server                           <none>     launch CRAC server,the output arguments will\n\
                                                 not be taken into account\n");
  printf("   --input-name-server                <STRING>   DEFAULT %s\n", INPUT_FIFO);
  printf("   --output-name-server               <STRING>   DEFAULT %s\n\n", OUTPUT_FIFO);
  printf("Additional settings for users\n");
  printf("  * Sam output file\n");  
  printf("   --detailed-sam                     <none>     more informations are added in SAM output file\n\n"); 
 
  printf("  * Mapping\n");
  printf("   --min-percent-single-loc           <FLOAT>    DEFAULT %.2f\n", PERCENT_MIN_SINGLE);
  printf("   --min-duplication                  <INT>      DEFAULT %d\n", MIN_LOCALISATION_DUPLICATION);
  printf("   --max-duplication                  <INT>      DEFAULT %d\n", MAX_LOCALISATION_DUPLICATION);
  printf("   --min-percent-duplication-loc      <FLOAT>    DEFAULT %.2f\n", PERCENT_MIN_DUPLICATE);
  printf("   --min-percent-multiple-loc         <FLOAT>    DEFAULT %.2f\n", PERCENT_MIN_MULTIPLE);
  printf("   --min-repetition                   <INT>      DEFAULT %d\n", MIN_OCC_REPETITION);  
  printf("   --min-percent-repetition-loc       <FLOAT>    DEFAULT %.2f\n\n", PERCENT_MIN_UNIQUE_REPETITION);
  
 
  printf("  * Biological causes\n");
  printf("   --max-splice-length                <INT>      DEFAULT %d\n", MAX_SPLICE_LENGTH);
  printf("   --max-bio-indel                    <INT>      DEFAULT %d\n", MAX_BIO_INS_DEL);
  printf("   --max-bases-retrieved              <INT>      DEFAULT %d\n\n", MAX_BASES_RETRIEVED);

  printf("  * Undetermined\n");
  printf("   --min-support-no-cover             <FLOAT>    DEFAULT %.2f\n\n", MIN_SUPPORT_NO_COVER);


  printf("Additional settings for advanced users\n");
  printf("  * Break verification and fusion (merging mirage breaks)\n");
  printf("   --min-break-length                 <FLOAT> DEFAULT %.2f\n", MIN_BREAK_LENGTH);
  printf("   --max-bases-randomly-matched       <INT>   DEFAULT %d\n", MAX_BASES_RANDOMLY_MATCHED);
  printf("   --max-extension-length             <INT>   DEFAULT %d\n\n", MAX_EXTENSION_LENGTH);

  printf("  * Threading\n");
  printf("   --nb-tags-info-stored              <INT>   DEFAULT %d\n\n", NB_TAGS_INFO_STORED);

  printf("  * Deep SNV search option\n");
  printf("   --nb-nucleotides-snv-comparison    <INT>   DEFAULT %d\n\n", NUMBER_NUCLEOTIDES_SNP_COMPARISON);
  
  // printf("   --percent-variation-coverage       <FLOAT> DEFAULT %.2f\n", PERCENT_SUPPORT_VARIATION_ALMOST_NORMAL);
  // printf("   --p-value-variation-coverage       <FLOAT> DEFAULT %.2f\n", P_VALUE_VARIATION_BIOLOGICAL);
  // printf("   --score-intra-ambiguous            <FLOAT> DEFAULT %.2f\n", SCORE_INTRA_AMBIGUOUS);
  // printf("   --score-inter-ambiguous            <FLOAT> DEFAULT %.2f\n", SCORE_INTER_AMBIGUOUS);
  // printf("   --min-bases-before-break           <INT>   DEFAULT %d\n", MIN_BASES_BEFORE_BREAK);
  // printf("   --support-score-window-length      <INT>   DEFAULT %d\n", SUPPORT_SCORE_WINDOW_LENGTH);
  // printf("   --max-support-out-no-cover         <FLOAT> DEFAULT %.2f\n", MAX_SUPPORT_OUT_NO_COVER);
  // printf("   --max-ambiguous-high-average       <FLOAT> DEFAULT %.2f\n", MAX_AMBIGUOUS_AVERAGE_HIGH);
  // printf("   --max-verynice-merge               <INT>   DEFAULT %d\n", MAX_VERYNICE_MERGE);
  // printf("   --nb-bases-substracted-low         <INT>   DEFAULT %d\n", NB_BASES_SUBSTRACTED_LOW_AVERAGE);
  // printf("   --nb-positions-check-ones          <INT>   DEFAULT %d\n", NB_POSITIONS_CHECK_ONES);
  // printf("   --min-perc-ones-inside             <FLOAT> DEFAULT %.2f\n", MIN_PERCENTAGE_ONES_INSIDE);
  // printf("   --min-ratio-support-fall           <FLOAT> DEFAULT %.2f\n", MIN_RATIO_SUPPORT_FALL);
  // printf("   --max-nb-overlapping-breaks        <INT>   DEFAULT %d\n\n", MAX_NB_OVERLAPPING_BREAKS);

}

void print_version() {
  printf("CRAC version " PACKAGE_VERSION "\n\n");
  printf("/******************************************************************************\n");
  printf("*                                                                             *\n");
  printf("*  Copyright © 2010-2013 -- IRB/INSERM                                        *\n");
  printf("*                           (Institut de Recherches en Biothérapie /          *\n");
  printf("*                           Institut National de la Santé et de la Recherche  *\n");
  printf("*                           Médicale)                                         *\n");
  printf("*                           LIFL/INRIA                                        *\n");
  printf("*                           (Laboratoire d'Informatique Fondamentale de       *\n");
  printf("*                           Lille / Institut National de Recherche en         *\n");
  printf("*                           Informatique et Automatique)                      *\n");
  printf("*                           LIRMM/CNRS                                        *\n");
  printf("*                           (Laboratoire d'Informatique, de Robotique et de   *\n");
  printf("*                           Microélectronique de Montpellier /                *\n");
  printf("*                           Centre National de la Recherche Scientifique)     *\n");
  printf("*                           LITIS                                             *\n");
  printf("*                           (Laboratoire d'Informatique, du Traitement de     *\n");
  printf("*                           l'Information et des Systèmes).                   *\n");
  printf("*                                                                             *\n");
  printf("*                                                                             *\n");
  printf("*  Auteurs/Authors: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *\n");
  printf("*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *\n");
  printf("*                   Thérèse COMMES   <commesd@univ-montp2.fr>                 *\n");
  printf("*                   Éric RIVALS      <eric.rivals@lirmm.fr>                   *\n");
  printf("*                                                                             *\n");
  printf("*  Programmeurs                                                               *\n");
  printf("*      /Progammers: Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>              *\n");
  printf("*                   Mikaël SALSON    <mikael.salson@lifl.fr>                  *\n");
  printf("*                   Jérôme AUDOUX    <jerome.audoux@univ-montp2.fr            *\n");
  printf("*  with additional contribution for the packaging of:	                      *\n");
  printf("*                   Alban MANCHERON  <alban.mancheron@lirmm.fr>               *\n");
  printf("*                                                                             *\n");
  printf("*  Contact:         CRAC list        <crac-bugs@lists.gforge.inria.fr>        *\n");
  printf("*                                                                             *\n");
  printf("/******************************************************************************\n");
}

void writeSamHeader(ostream *sam
          , LocateOnGenome *genome
          , int argc
          , char **argv
          , string sam_filename) {
  *sam << "@HD\tVN:1.4\tSO:unsorted" << endl;
  // Display information on sequence names and lengths
  for (uint i=0; i < genome->getGenomeInfo()->getNbChr(); i++) {
    *sam << "@SQ\tSN:" << genome->getGenomeInfo()->getChrName(i) 
         << "\tLN:" << genome->getGenomeInfo()->getChrLength(i) << endl;
  }
  struct tm tm;
  char tmp[120];
  time_t t;
  time(&t);
  localtime_r(&t,&tm); 
  sprintf(tmp,"%i-%02i-%02iT%02i:%02i:%02i", tm.tm_year + 1900, tm.tm_mon+1, tm.tm_mday, tm.tm_hour, tm.tm_min, tm.tm_sec);
  *sam << "@RG\tID:" << rand() << "\tSM:" << sam_filename << "\tDT:" << tmp << endl;
  *sam << "@PG\tID:crac\tPN:crac\tVN:" << PACKAGE_VERSION;
  if (argc) {
    *sam << "\tCL:";
    for (int i=0; i < argc; i++)
      *sam << argv[i] << " ";
    *sam << endl;
  }
}

int main(int argc, char **argv) {
  fstream file;
  int nb_streams;
  ostream *streams[19] = {NULL,NULL,NULL
			  ,NULL,NULL,NULL
			  ,NULL,NULL,NULL
			  ,NULL,NULL,NULL
			  ,NULL,NULL, NULL
                          ,NULL,NULL,NULL
			  ,NULL};
  const char *extensions[19] = {
    "snp", "indel", "error", "splice", 
    "coverlessSplice", "chimera", "undetermined", "repetition", 
    "duplication", "nothing", "normal", "almostNormal", "multiple",
    "none", "bioUndetermined", "single", "sam", "pairedEndChimera", 
    "summary"
  };
  nb_streams = 19;

  bool gzOutput = false;
  gkArrays *indexTags = NULL;
  ulong pattern_length = 0;
  bool is_stranded = false;
  ClassifyTags *ct;
  LocateOnGenome *log = NULL;
  Parameters *params = new Parameters();
  params->show_progressbar = false;
  uint threshold = ~0;
  char *reads = 0;
  char *reads2 = 0;
  char *indexGenome = 0;
  char *fileBasename = 0;
  char *filename = 0;
  uint max_length = 0;
  int nb_streams_to_init = 0;
  bool use_server = false;
  bool use_emt = false;
  bool use_paired_end_chimera = false;
  bool use_all_outputs = false;
  bool find_output_file = false;
  char *fifo_input_name = NULL;
  char *fifo_output_name = NULL;
  string sam_filename;
  ostream *summary = NULL;

  int option_index = 0;
  int opt;
  
  static const char *short_options =  "hi:r:o:k:m:fv";
  
  static struct option long_options[] = {
    { "help",    no_argument,       NULL, 'h' },
    { "full-help",    no_argument,       NULL, 'f' },
    { "version",    no_argument,       NULL, 'v' },
    { "all", required_argument, NULL, ALL_OUTPUT},
    { "reads-length", required_argument, NULL, READS_LENGTH},
    { "stranded", no_argument, NULL, IS_STRANDED_CST},
    { "sam",    required_argument, NULL, SAM },
    { "snv",    required_argument, NULL, SNP },
    { "indel",    required_argument, NULL, BIOTAG },
    { "errors",    required_argument, NULL, SEQERR },
    { "splice",    required_argument, NULL, SPLICE },
    { "weak-splice",    required_argument, NULL, SPLICE_NOCOVER },
    { "chimera",    required_argument, NULL, CHIMERA },
    { "undetermined",    required_argument, NULL, UNDETERMINED },
    { "repeat",    required_argument, NULL, REPETITION },
    { "duplicate",    required_argument, NULL, DUPLICATION },
    { "nothing",    required_argument, NULL, NOTHING },
    { "normal",    required_argument, NULL, NORMAL },
    { "almost-normal",    required_argument, NULL, ALMOSTNORMAL },
    { "multiple",    required_argument, NULL, MULTIPLE},
    { "none",    required_argument, NULL, NONE },
    { "biological",    required_argument, NULL, BIOUNDETERMINED },
    { "single",    required_argument, NULL, SINGLE },
    { "max-locs", required_argument, NULL, MAX_NB_LOCATED_OCCS_CST},
    { "max-splice-length",    required_argument, NULL, MAX_SPLICE },
    { "max-bio-indel",    required_argument, NULL, MAX_BIO_INDEL },
    { "min-duplication",    required_argument, NULL, MIN_DUPLICATION },
    { "max-duplication",    required_argument, NULL, MAX_DUPLICATION },
    { "min-percent-multiple-loc",    required_argument, NULL, PERCENT_MIN_MULTIPLE_CST },
    { "min-percent-duplication-loc",    required_argument, NULL, PERCENT_MIN_DUPLICATION },
    { "min-percent-single-loc", required_argument, NULL, PERCENT_MIN_SINGLE_CST},
    { "min-percent-repetition-loc",    required_argument, NULL, PERCENT_MIN_SINGLE_REPETITION },
    { "max-verynice-merge",    required_argument, NULL, MAX_VERYNICE_MERGE_CST },
    { "min-repetition",    required_argument, NULL, MIN_LOC_REPETITION },
    { "percent-variation-coverage",    required_argument, NULL, PERCENT_VARIATION_COVERAGE },
    { "p-value-variation-coverage",    required_argument, NULL, P_VALUE_VARIATION_COVERAGE },
    { "score-intra-ambiguous", required_argument, NULL, SCORE_INTRA_AMBIGUOUS_CST},
    { "score-inter-ambiguous", required_argument, NULL, SCORE_INTER_AMBIGUOUS_CST},
    { "max-bio-deviation",    required_argument, NULL, MAX_BIOLOGICAL_DEVIATION },
    { "max-bases-randomly-matched",    required_argument, NULL, MAX_RANDOMLY_MATCHED },
    { "min-bases-before-break",    required_argument, NULL, MIN_BEFORE_BREAK },
    { "max-bases-retrieved",    required_argument, NULL, MAX_RETRIEVED },
    { "max-extension-length",    required_argument, NULL, MAX_EXTENSION },
    { "support-score-window-length",    required_argument, NULL, SUPPORT_SCORE_WINDOW },
    { "min-support-no-cover",    required_argument, NULL, MINIMUM_SUPPORT_NO_COVER },
    { "min-break-length",    required_argument, NULL, MINIMUM_BREAK_LENGTH },
    { "max-support-out-no-cover",    required_argument, NULL, MAXIMUM_SUPPORT_OUT_NO_COVER },
    { "max-ambiguous-high-average",    required_argument, NULL, MAX_AMBIGUOUS_AVG_HIGH },
    { "nb-bases-substracted-low",    required_argument, NULL, NB_BASES_SUBSTRACTED_LOW_AVG },    
    { "nb-positions-check-ones",    required_argument, NULL, NB_POS_CHECK_ONES },
    { "min-perc-ones-inside",    required_argument, NULL, MIN_PERC_ONES_INSIDE },    
    { "min-ratio-support-fall",    required_argument, NULL, MIN_RATIO_SUPPORT_FALL_CST },
    { "nb-tags-info-stored", required_argument, NULL, NB_TAGS_INFO_CST},
    { "nb-threads", required_argument, NULL, NB_THREADS_CST},
    { "max-nb-overlapping-breaks", required_argument, NULL, 
      MAX_NB_OVERLAPPING_BREAKS_CST},
    { "deep-snv", no_argument, NULL, DEEP_SNP_SEARCH_CST},
    { "nb-nucleotides-snv-comparison", required_argument, NULL, 
      NUMBER_NUCLEOTIDES_SNP_COMPARISON_CST},
    { "detailed-sam", no_argument, NULL, DETAILED_SAM_CST},
    { "server", no_argument, NULL, SERVER_CST},
    { "input-name-server", required_argument, NULL, INPUT_FIFO_CST},
    { "output-name-server", required_argument, NULL, OUTPUT_FIFO_CST},
    { "emt", no_argument, NULL, DGE },
    { "stringent-chimera", no_argument, NULL, STRINGENT_CHIMERA_CST },
    { "no-ambiguity", no_argument, NULL, NO_AMBIGUITY_CST },
    { "paired-end-chimera", required_argument, NULL, PAIRED_END_CHIMERA },
    { "gz", no_argument, NULL, GZ_OUTPUT},
    { "summary",    required_argument, NULL, SUMMARY },
#ifdef HAVE_LIBPROGRESSBAR
    { "show-progressbar", no_argument, NULL, PROGRESSBAR },
#endif
    { "treat-multiple", no_argument, NULL, TREAT_MULTIPLE_CST },
    { 0, 0, 0, 0 } //terminator
  };
  
  while ((opt = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1){
    switch(opt) {
      
    case 'h':
      print_help();
      exit(0);
    case 'f' : 
      print_fullhelp();
      exit(0);
    case 'v' : 
      print_version();
      exit(0);
    case 'i':
      indexGenome = optarg;
      break;
    case 'r':
      reads = optarg;   
      break;
    case 'k':
      threshold = atoi(optarg); 
      params->threshold = threshold;
      break;
    case 'o':
      sam_filename = (const char *)optarg;
      if (strcmp(sam_filename.c_str(),"-") != 0){ 
	streams[16] = create_stream(gzOutput, optarg); 
      }else{
	streams[16] = &cout;
      }
      break;
    case 'm':
      pattern_length = atol(optarg);
      break;
    case ALL_OUTPUT:
      fileBasename = optarg;
      use_all_outputs = true;
      break;
    case READS_LENGTH:
      pattern_length = atol(optarg);
      break;
    case IS_STRANDED_CST:
      is_stranded = true;
      break;
    case MAX_NB_LOCATED_OCCS_CST:
      params->max_nb_located_occs = atol(optarg);
      break;
    case SNP:
      streams[0] = create_stream(gzOutput, optarg);
      break;
    case BIOTAG:
      streams[1] = create_stream(gzOutput, optarg);
      break;
    case SEQERR:
      streams[2] = create_stream(gzOutput, optarg);
      break;
    case SPLICE:
      streams[3] = create_stream(gzOutput, optarg);
      break;
    case SPLICE_NOCOVER:
      streams[4] = create_stream(gzOutput, optarg);
      break;
    case CHIMERA:
      streams[5] = create_stream(gzOutput, optarg);
      break;
    case UNDETERMINED:
      streams[6] = create_stream(gzOutput, optarg);
      break;
    case REPETITION:
      streams[7] = create_stream(gzOutput, optarg);
      break;
    case DUPLICATION:
      streams[8] = create_stream(gzOutput, optarg);
      break;
    case NOTHING:
      streams[9] = create_stream(gzOutput, optarg);
      break;
    case NORMAL:
      streams[10] = create_stream(gzOutput, optarg);
      break;
    case ALMOSTNORMAL:
      streams[11] = create_stream(gzOutput, optarg);
      break;
    case MULTIPLE:
      streams[12] = create_stream(gzOutput, optarg);
      break;
    case NONE:
      streams[13] = create_stream(gzOutput, optarg);
      break;
    case BIOUNDETERMINED:
      streams[14] = create_stream(gzOutput, optarg);
      break;
    case SINGLE:
      streams[15] = create_stream(gzOutput, optarg);
      break;
    case SAM:
      streams[16] = create_stream(gzOutput, optarg);
      sam_filename = (const char *)optarg;
      break;
    case PAIRED_END_CHIMERA:
      streams[17] = create_stream(gzOutput, optarg);
      use_paired_end_chimera = true;
      break;
    case SUMMARY:
      streams[18] = create_stream(gzOutput, optarg);
      summary = streams[18]; 
      break;
    case MAX_SPLICE:
      params->max_splice_length = atoi(optarg);
      break;
    case MAX_BIO_INDEL:
      params->max_bio_ins_del = atoi(optarg);
      break;
    case MAX_DUPLICATION:
      params->max_localisation_duplication = atoi(optarg);
      break;
    case MIN_DUPLICATION:
      params->min_localisation_duplication = atoi(optarg);
      break;
    case PERCENT_MIN_SINGLE_REPETITION:
      params->percent_min_unique_repetition = atof(optarg);
      break;
    case PERCENT_MIN_SINGLE_CST: 
      params->percent_min_single = atof(optarg);
      break;
    case PERCENT_MIN_MULTIPLE_CST:
      params->percent_min_multiple = atof(optarg);
      break;
    case PERCENT_MIN_DUPLICATION:
      params->percent_min_duplicate = atof(optarg);
      break;
    case MIN_LOC_REPETITION:
      params->min_occ_repetition = atoi(optarg);
      break;
    case PERCENT_VARIATION_COVERAGE:
      params->percent_support_variation_almost_normal = atof(optarg);
      break;
    case P_VALUE_VARIATION_COVERAGE:
      params->p_value_variation_biological = atof(optarg);
      break;
    case SCORE_INTRA_AMBIGUOUS_CST:
      params->score_intra_ambiguous = atof(optarg);
      break;
    case SCORE_INTER_AMBIGUOUS_CST:
      params->score_inter_ambiguous = atof(optarg);
      break;
    case MAX_RANDOMLY_MATCHED:
      params->max_bases_randomly_matched = atoi(optarg);
      break;
    case MAX_VERYNICE_MERGE_CST:
      params->max_verynice_merge = atoi(optarg);
      break;
    case MIN_BEFORE_BREAK:
      params->min_bases_before_break = atoi(optarg);
      break;
    case MAX_RETRIEVED:
      params->max_bases_retrieved = atoi(optarg);
      break;
    case MAX_EXTENSION:
      params->max_extension_length = atoi(optarg);
      break;
    case SUPPORT_SCORE_WINDOW:
      params->support_score_window_length = atoi(optarg);
      break;
    case MINIMUM_SUPPORT_NO_COVER:
      params->min_support_no_cover = atof(optarg);
      break;
    case MINIMUM_BREAK_LENGTH:
      params->min_break_length = atoi(optarg);
      break;
    case MAXIMUM_SUPPORT_OUT_NO_COVER:
      params->max_support_out_no_cover = atof(optarg);
      break;
    case MAX_AMBIGUOUS_AVG_HIGH:
      params->max_ambiguous_average_high = atof(optarg);
      break;
    case NB_BASES_SUBSTRACTED_LOW_AVG:
      params->nb_bases_substracted_low_average = atoi(optarg);
      break;
    case NB_POS_CHECK_ONES:
      params->nb_positions_check_ones = atoi(optarg);
      break;
    case MIN_PERC_ONES_INSIDE:
      params->min_perc_ones_inside = atof(optarg);
      break;
    case MIN_RATIO_SUPPORT_FALL_CST:
      params->min_ratio_support_fall = atof(optarg);
      break;
    case NB_TAGS_INFO_CST:
      params->nb_tags_info_stored = atoi(optarg);
      break;
    case NB_THREADS_CST:
      params->nb_threads = atoi(optarg);
      break;
    case MAX_NB_OVERLAPPING_BREAKS_CST:
      params->max_nb_overlapping_breaks = atof(optarg);
      break;
    case DEEP_SNP_SEARCH_CST:
      params->deep_snp_search = true;
      break;
    case NUMBER_NUCLEOTIDES_SNP_COMPARISON_CST:
      params->number_nucleotides_snp_comparison = atoi(optarg);
      break;
    case DETAILED_SAM_CST:
      params->detailed_sam = true;
      break;
    case SERVER_CST:
      use_server=true;
      break;
    case STRINGENT_CHIMERA_CST:
      params->stringent_chimera = true;
      break;
    case INPUT_FIFO_CST:
      fifo_input_name = new char[strlen(optarg)+1];
      strcpy(fifo_input_name, optarg);
      break;
    case OUTPUT_FIFO_CST:
      fifo_output_name = new char[strlen(optarg)+1];
      strcpy(fifo_output_name, optarg);
      break;
    case DGE:
      use_emt=true;
      break;
    case GZ_OUTPUT:
      gzOutput = true;
      break;
    case NO_AMBIGUITY_CST:
      params->no_ambiguity=true;
      break;
    case PROGRESSBAR:
      params->show_progressbar = true;
      break;
    case TREAT_MULTIPLE_CST:
      params->treat_multiple=true;
      break;
    case '?':
      print_help();
      exit(1);
    default:
      print_help();
      exit(1);
    }
  }

  // Init parameters with values we have fixed with options
  params->init();

  // little hack to get the second file argument of -r option when doing paired-end
  if (optind < argc)
  {
      reads2 = argv[optind++];
  }

  // Init output files if we use --all option
  if(use_all_outputs) {
    nb_streams_to_init = nb_streams;
    // If we are not using paired-end reads we dont create the last stream which is
    // dedicated for paired-end reads only
    if(reads2 == NULL) {
      nb_streams_to_init--;
    }
    max_length = 0;
    for (int i = 0; i < nb_streams_to_init; i++) {
      uint current_length = strlen(extensions[i]);
      if (current_length > max_length)
        max_length = current_length;
    }
    filename = new char[strlen(fileBasename)+max_length+2];
    for (int i=0; i < nb_streams_to_init; i++) {
      sprintf(filename, "%s.%s", fileBasename, extensions[i]);
      streams[i] = create_stream(gzOutput, filename);
      if (i==16)
	sam_filename = filename;
    }
    delete [] filename;
  }

  if (indexGenome == NULL  || reads == 0 || 
      threshold == (uint)~0 || (streams[16] == NULL && !use_server)){
    print_help();
    cerr << endl << "Missing at least a mandatory argument among: -i <index> -r <reads> -k <k_length> -o <output_SAM>" << endl << endl;    
    exit(1);
  }else{
    for (int i=0; i < nb_streams && !find_output_file ; i++) {
      if (streams[i] != NULL){
	find_output_file = true;
      }
    }
    if (!find_output_file && !use_server){
      print_help();
      cerr << endl << "Obviously, you must specify at least one output argument" << endl << endl;
      exit(1);
    }
    // check if the there is paired-end input file when using paired-end options
    if (reads2 == NULL && use_paired_end_chimera) {
       print_help();
       cerr << endl << "Please specify a second read input file if you want to use a specific option for paired-end tags." << endl;
    }
  }

  // Checking fifos for server
  if (use_server) {
    struct stat buf;
    // If no name was specified, use default ones.
    if (! fifo_input_name)
      fifo_input_name = (char *)INPUT_FIFO;
    if (! fifo_output_name)
      fifo_output_name = (char *)OUTPUT_FIFO;

    if (stat(fifo_input_name,&buf) == 0 || stat(fifo_output_name,&buf) == 0) {
      cerr << fifo_input_name << " or " << fifo_output_name << " already exist."
	   << endl << "Please specity another name or delete the existing files." << endl;
      exit(4);
    }
  }

  // Launching exact-matching for reads mapping
  if (use_emt){
    if (reads2 != NULL)
      cerr << "Warning: CRAC--emt version does not support paired-ends protocol, only the first input file "<<reads <<" is considered" << endl;
    if (summary != NULL)
      *summary << "Running CRAC-emt version " <<PACKAGE_VERSION <<" optimized for exact mapping reads " << endl << endl;
    // Indexation of the genome
    if (threshold > 0)
      cerr << "Warning: Running CRAC--emt with k=" << threshold << " which means that it only maps the first located kmer of length "<< threshold <<" for each read" << endl;
    else
      cerr << "Warning: Running CRAC--emt with k=0 to map the entire reads" << endl;
    log = new LocateOnGenome((uchar *)indexGenome); 
    log->setNbLocations(params->max_nb_located_occs);
    log->startChrono();
    LocateTags *lt; 
    lt = new LocateTags(log, params, reads, threshold);
    // put headers in the sam file
    writeSamHeader(streams[16],log,argc,argv,sam_filename);
    // launch the locate process
    lt->locate(streams[16]);
    
    if (summary != NULL) {
      *summary << "Time to locate all tags: "<< log->getElapsedTime() << " s" << endl;       
      *summary << endl << "----------------------------------" << endl;
      *summary << "           Some STATISTICS            " << endl;
      *summary << endl << "----------------------------------" << endl;
      *summary << "Total number of reads analyzed: " << lt->getNbTags() << endl << endl;
      *summary << "Single: " << lt->getNbSingle() << " (" << lt->getNbSingle()*100./lt->getNbTags() << "%)" << endl;
      *summary << "Multiple: " << lt->getNbMultiple() << " (" << lt->getNbMultiple()*100./lt->getNbTags() << "%)" << endl;
      *summary << "None: " << lt->getNbNone() << " (" << lt->getNbNone()*100./lt->getNbTags() << "%)" << endl;
    }
    delete lt;
  }
  // Launching CRAC for RNA-Seq mapping 
  else{
    // Indexation of reads
    float start = getChrono();
    // If we have pairedEnd reads we create gkArrays using both reads files
    try {
      if (reads2 != 0) {
        if (summary != NULL)
          *summary << "Running CRAC using paired-end files." << endl;
        indexTags = new gkArrays(reads, reads2, threshold, false, pattern_length, is_stranded, params->nb_threads, params->show_progressbar);
      } else {
        indexTags = new gkArrays(reads, threshold, false, pattern_length, is_stranded, params->nb_threads, params->show_progressbar);
      }
    } catch (bad_alloc e) {
      cerr << "You do not have enough memory to store all the reads. Aborting" << endl;
      exit(42);
    }

    if (summary != NULL){
      if (indexTags->isLarge()){
	*summary << "Running CRAC version " <<PACKAGE_VERSION <<" optimized for large data set ";
	if (pattern_length){
	  *summary << "and fixed read length " ;
	}else{
	  *summary << "and variable read length ";
	}
      }else{
	*summary << "Running CRAC version " <<PACKAGE_VERSION <<" optimized for small data set < 4Gbp (eg. less than 50M reads of 75 bp) ";
	if (pattern_length){
	  *summary << "and fixed read length " ;
	}else{
	  *summary << "and variable read length " ;
	}
      }
      if (is_stranded)
	*summary << "for a strand specific library." << endl;
      else
	*summary << "for a non-strand specific library." << endl;
      
      *summary << endl <<  "Time to index all reads: " << (getChrono()-start)/1000000. << " s" << endl;
    }      

    // Indexation of the genome
    log = new LocateOnGenome((uchar *)indexGenome); 
    log->setNbLocations(params->max_nb_located_occs);
    
    // classify the reads
    ct = new ClassifyTags(log, indexTags, params);
    if (use_server) {
      locateServer(ct, fifo_input_name, fifo_output_name);
    } else {
      log->startChrono();
      // put headers in the sam file
      writeSamHeader(streams[16],log,argc,argv,sam_filename);
      // put headers with feature for each files
      ct->getHeaders(streams[0], streams[1], streams[2], streams[3],
		     streams[4], streams[5], streams[6], streams[7],
		     streams[8], streams[9], streams[10], streams[11],
		     streams[12], streams[13], streams[14], streams[15],
		     streams[17]);
      // and classify
      ct->classify(streams[0], streams[1], streams[2], streams[3],
		   streams[4], streams[5], streams[6], streams[7],
		   streams[8], streams[9], streams[10], streams[11],
		   streams[12], streams[13], streams[14], streams[15],
		   streams[16], streams[17] );
     
      if (summary != NULL){ 
	*summary << "Time to classify all reads: "<< log->getElapsedTime() << " s" << endl; 
	
	// making statistics on STDOUT
	*summary << endl << "----------------------------------" << endl;
	*summary << "           Some STATISTICS            " << endl;
	*summary << endl << "----------------------------------" << endl;
	*summary << "Total number of reads analyzed: " << indexTags->getNbTags() << endl << endl;
	*summary << "Single: " << ct->getNbSingle() << " (" << ct->getNbSingle()*100./indexTags->getNbTags() << "%)" << endl;
	*summary << "Multiple: " << ct->getNbMultiple() << " (" << ct->getNbMultiple()*100./indexTags->getNbTags() << "%)" << endl;
	*summary << "None: " << ct->getNbNone() << " (" << ct->getNbNone()*100./indexTags->getNbTags() << "%)" << endl;
	*summary << "Duplication: " << ct->getNbDuplication() << " (" << ct->getNbDuplication()*100./indexTags->getNbTags() << "%)" << endl;
	
	*summary << endl << "Warning: the sum of the four previous categories may not be equal to 100%." << endl << " This is normal: reads are considered by chunks." << endl << " In a given read, a chunk may appear multiple times, while another just appears once." << endl ; 
	*summary << endl << "----------------------------------" << endl;
	
	*summary << "Explainable: " << ct->getNbExplainable() << " (" << ct->getNbExplainable()*100./indexTags->getNbTags() << "%)" << endl << endl;
	
	*summary << "Repetition: " << ct->getNbRepetition() << " (" << ct->getNbRepetition()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Normal: " << ct->getNbNormal() << " (" << ct->getNbNormal()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Almost-Normal: " << ct->getNbAlmostNormal() << " (" << ct->getNbAlmostNormal()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Sequence-Errors: " << ct->getNbSeqErr() << " (" << ct->getNbSeqErr()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "SNV: " << ct->getNbSNP() << " (" << ct->getNbSNP()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Short-Indel: " << ct->getNbBioTagIndel() << " (" << ct->getNbBioTagIndel()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Splice: " << ct->getNbSplice() << " (" << ct->getNbSplice()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Weak-Splice: " << ct->getNbSpliceNoCover() << " (" << ct->getNbSpliceNoCover()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Chimera: " << ct->getNbSpliceInter() << " (" << ct->getNbSpliceInter()*100./indexTags->getNbTags() << "%)" << endl ;
  if(reads2 != NULL) {
    *summary << "Paired-end Chimera: " << ct->getNbPairedEndChimera() << " (" << ct->getNbPairedEndChimera()*100./indexTags->getNbTags() << "%)" << endl ;
  }
	*summary << "Bio-Undetermined: " << ct->getNbBioUndetermined() << " (" << ct->getNbBioUndetermined()*100./indexTags->getNbTags() << "%)" << endl ;
	*summary << "Undetermined: " << ct->getNbUndetermined() << " (" << ct->getNbUndetermined()*100./indexTags->getNbTags() << "%)" << endl ;
	// *summary << "Nothing: " << ct->getNbNothing() << " (" << ct->getNbNothing()*100./indexTags->getNbTags() << "%)" << endl ;
      }
    }  
    delete ct;
  }
  
  for (int i = 0; i < nb_streams; i++) {
    if (streams[i] != NULL && streams[i] != &cout)
      delete streams[i]; 
  }
  delete indexTags;
  delete params;
  delete log;
  return 0;
}
