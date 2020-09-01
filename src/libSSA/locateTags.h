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

#ifndef LOCATETAGS_H
#define LOCATETAGS_H

#include "locateOnGenome.h"
#include "../Parameters.h"
#include <ostream>
#include <string>
#include <readsReader.h>
#include <semaphore.h>

using namespace gkarrays;

typedef struct read_s {
  string name; /* Name of the read in the input file */
  string seq;  /* DNA sequence of the read */
  string qual; /* Quality of the DNA sequence */
} read_t;

class LocateTags{
 private:
  LocateOnGenome *genome; /* Genome object to locate reads on reference genome */
  Parameters *parameters; /* Parameters object */
  char *tags_file;
  uint nb_single;         /* Number of single locations */
  uint nb_multiple;       /* Number of multiple locations */
  uint nb_none;           /* Number of None locations */
  uint nb_tags;           /* Number of tags in the file */
  uint threshold;         /* Threshold to use */
  sem_t stats_mutex;      /* Mutex when thread updates stats values (nb_single, nb_multiple...) */
  
 public: 
 
  /**
   * Constructor of locateTags
   *
   * @param genome the genome object used for getting locations
   * @param param the parameters object
   * @param tags_file the name of the FAST/FASTQ file containing reads to process
   * @param sam output stream to print the sam file
   * @param threshold the size of the prefix of the read we are going to use for mapping
   *                  if the size equal 0 locateTags will map the whole sequence of the read.
   */
  LocateTags(LocateOnGenome *genome, Parameters *param, char *tags_file, uint threshold = 0);
  
  /**
   * Destructor of locateTags
   */
  ~LocateTags(){};

  /**
   * Launch location process, and multiple threads if any.
   */
  void locate(ostream *sam);
  
  /**
   * locate the tag in the genome index and generate the coresponding sam line.
   *
   * @param read the read to process
   * @return a string containing the sam line of the mapping of this read
   */
  string *locateTag(read_t *read);
  
  /**
   * @return the number of tags that have been processed
   */
  uint getNbTags();

  /**
   * @return the number of single locations
   */
  uint getNbSingle();
  
  /**
   * @return the number of multiple locations
   */
  uint getNbMultiple();
  
  /**
   * @return the number of none locations
   */
  uint getNbNone();

 private:
  /**
   * Set a read_t object with info from iterato
   *
   * @readIt the iterator used to fill data of the read object
   * @read the read to fill
   */
  void setReadInfo(readIterator *readIt, read_t *read);

  /**
   * The function that is call for each thread
   */
  friend void *fillLocations(void *);

};  

void *fillLocations(void *);

typedef struct locate_thread_s {
  LocateTags *locate;
  sem_t *remaining_seats;  /* Number of remaining cells containing reads to process */
  sem_t *processed_reads;  /* Number of reads that have been processed and ready to write in the ouput */
  string **seats;          /* Seats where the samLines of the processed reads are stored */
  read_t **reads;          /* Seats where the reads to process are stored */
  uint total_seats;        /* Total number of seats */
  bool keepGoing;          /* Boolean to stop the tread */
} locate_thread_t;

#endif
