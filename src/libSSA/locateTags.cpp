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

#include <cstdlib>
#include "locateTags.h"
#include "utils.h"
#include <config.h>
#include <sstream>
#include <algorithm>
#include "../libReadsInfo/samLine.h"

using namespace std;

LocateTags::LocateTags(LocateOnGenome *genome, Parameters *param, char *tags_file, uint threshold):
  genome(genome),parameters(param),tags_file(tags_file),threshold(threshold)
{}

void LocateTags::locate(ostream *sam) {
  // Init stat variables
  nb_single = 0;
  nb_multiple = 0;
  nb_none = 0;
  nb_tags = 0;
  sem_init(&stats_mutex,0,1);

  pthread_t threads[parameters->nb_threads];
  locate_thread_t params[parameters->nb_threads];
  string **store_samLines[parameters->nb_threads-1];
  read_t **store_reads[parameters->nb_threads-1];
  readsReader reads(tags_file,threshold);
  readIterator *readItForBuffer = reads.begin();
  readIterator *readIt = reads.begin();
  sem_t remaining_seats[parameters->nb_threads-1];
  sem_t processed_reads[parameters->nb_threads-1];
  uint *positions = NULL;
  uint nb_samLines;
  string *samLine;

  if (parameters->nb_threads >  1) {
    nb_samLines = parameters->nb_tags_info_stored;

    positions = new uint[parameters->nb_threads - 1];

    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      positions[j] = 0;
      params[j].locate = this;
      sem_init(&processed_reads[j],0,0);
      params[j].processed_reads = &processed_reads[j];
      params[j].total_seats = nb_samLines;
      store_reads[j] = (read_t **)malloc(nb_samLines*sizeof(read_t *));
      params[j].reads = store_reads[j];
      store_samLines[j] = (string **)malloc(nb_samLines*sizeof(string *));
      params[j].seats = store_samLines[j];
      params[j].keepGoing = true;
    }
  } else {
    nb_samLines = 1;
  }

  if (parameters->nb_threads > 1) {
    // add reads data in threads
    uint i = 0;
    uint j = 0;
    uint readCpt = 0;
    while (!readItForBuffer->isFinished() && i < nb_samLines) {
      j = 0;
      while (j < (parameters->nb_threads -1) && !readItForBuffer->isFinished()) {
        read_t *r = new read_t();
        setReadInfo(readItForBuffer,r);
        store_reads[j][i] = r;
        ++(*readItForBuffer);
        readCpt++;
        j++;
      }
      i++;
    }
    // Init the semaphore according to the number of reads stored
    // in each thread
    for (uint k = 0; k < parameters->nb_threads -1; k++) {
      if (k >= j) {
        sem_init(&remaining_seats[k],0,i-1);
        params[k].remaining_seats = &remaining_seats[k];
      } else {
        sem_init(&remaining_seats[k],0,i);
        params[k].remaining_seats = &remaining_seats[k];
      }
    }

    // Launching threads
    for (uint j = 0; j < parameters->nb_threads - 1; j++) {
      pthread_create(&threads[j], NULL, fillLocations, &params[j]);
    }
  } 

  uint i = 0;
  while (!readIt->isFinished()) {
    if (parameters->nb_threads > 1) {
      // get the number of the thread we are waiting for
      uint j = i % (parameters->nb_threads - 1);
      // wait until the read has been processed
      sem_wait(&processed_reads[j]);
      // get the sam line resulting of the read location
      samLine = store_samLines[j][positions[j]];
      // delete the read object that has been used to get the location
      delete store_reads[j][positions[j]];
      // Add a new read to the thread buffer
      if (!readItForBuffer->isFinished()) {
        read_t *r = new read_t();
        params[j].reads[positions[j]] = r;
        setReadInfo(readItForBuffer,params[j].reads[positions[j]]);
        ++(*readItForBuffer);
        // Tell the tread that it can locate this new read
        sem_post(&remaining_seats[j]);
      }
      // Update the position of the next read to get from
      // this thread
      positions[j] = (positions[j] + 1) % nb_samLines;
    } else {
      // If we dont use multi-threading we directly compute
      // sam line for the current read
      read_t *read = new read_t();
      setReadInfo(readIt,read);
      samLine = locateTag(read);
      delete read;
    }
    // Write the sam line in the file
    *sam << *samLine;
    delete samLine;
    ++(*readIt);
    i++;
    nb_tags++;
  }

  // unlock threads
  for (uint j = 0; j < parameters->nb_threads -1; j++) {
    params[j].keepGoing = false;
    sem_post(&remaining_seats[j]);
  }
  
  // deleting iterators
  delete readIt;
  delete readItForBuffer;

  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    free(store_samLines[j]);
    free(params[j].reads);
    sem_destroy(&remaining_seats[j]);
    sem_destroy(&processed_reads[j]);
  }

  for (uint j = 0; j < parameters->nb_threads - 1; j++) {
    pthread_join(threads[j], NULL);
  }

  if (positions)
    delete [] positions;

}

string *LocateTags::locateTag(read_t *read) {
  ostringstream os;
  pair<ChrPosition **,uint>locs;
  uint pos_locs;
  string tag_s = read->seq; 
  ulong tag_length = (ulong)tag_s.size();
  SamLine base_line;
  // minuscule nucleotide to capital
  transform(tag_s.begin(), tag_s.end(), tag_s.begin(), static_cast<int(*)(int)>(toupper));

  uchar *tag = (uchar *)tag_s.c_str();
  string *samLines = new string();
  
  if (threshold == 0){
    locs = genome->getLocations((uchar *)&tag[0], tag_length, tag_length, pos_locs);
  }else{
    locs = genome->getLocations((uchar *)&tag[0], threshold, tag_length, pos_locs);
  }

  // write generic SAM features
  base_line.setQname(read->name);
  base_line.setSeq(read->seq);
  base_line.setQual(read->qual);
  
  // write sam_line, if no loc
  if (locs.second == 0) {
    base_line.setSegmentUnmapped();
    base_line.setMapQ(255);
    base_line.addOptionalField("NH",locs.second);
    base_line.writeLine(os);
    // set mutex for none stats
    sem_wait(&stats_mutex);
    nb_none++;
    sem_post(&stats_mutex);
    // end mutex zone
  }
  // write sam_line for each loc (or a single loc)
  else{    
    for(ulong i = 0; i < locs.second; i++) {
      if (i==0 || parameters->treat_multiple){
	SamLine current_line = base_line;
	string rname = string(locs.first[i]->getChrPosition());
	current_line.setRname(rname);
	current_line.setPos(locs.first[i]->getRelativePosition());
	current_line.setMapQ(254);
	
	Cigar cigar;
	if (threshold > 0){
	  uint soft_size = tag_length-threshold-pos_locs;
	  if (locs.first[i]->getStrand() == 1){
	    if (pos_locs > 0)
	      cigar.append(pos_locs,'S');
	    cigar.append(threshold,'M');
	    if (soft_size  > 0)
	      cigar.append(soft_size,'S');
	  }else{
	    if (soft_size > 0)
	      cigar.append(soft_size,'S');
	    cigar.append(threshold,'M');
	    if (pos_locs > 0)
	      cigar.append(pos_locs,'S');
	  }
	}else
	  cigar.append(tag_length,'M');
	current_line.setCigar(cigar);
	if (locs.first[i]->getStrand() == -1){
	  current_line.setSeqReverseComplemented();
	  current_line.reverseComplementeSeq();
	  current_line.reverseQual();
	}
	//multiple_alignement.writeLine(cout);
	// If this is not the primary line
	if(i>0) {
	  current_line.setSecondaryAlignement();
	}
	current_line.addOptionalField("NH",locs.second);
	// set single and multiple and update stats
	// use a mutex to protect globa stats update
	sem_wait(&stats_mutex);
	if (locs.second == 1)
	  nb_single++;
	else 
	  nb_multiple++;
	sem_post(&stats_mutex);
	// end of mutex zone
	current_line.writeLine(os);
      }
      delete locs.first[i];
    }
    free(locs.first);
  }
  *samLines = os.str();  
  //delete [] tag;
  return samLines;
}

uint LocateTags::getNbTags(){
  return nb_tags;
}

uint LocateTags::getNbSingle(){
  return nb_single;
}

uint LocateTags::getNbMultiple(){
  return nb_multiple;
}

uint LocateTags::getNbNone(){
  return nb_none;
}

// private methods

void LocateTags::setReadInfo(readIterator *readIt,read_t *read) {
  if(readIt->getName() != NULL) {
    read->name = readIt->getName();
  } else { read->name = "*"; }
  read->seq = readIt->getSequence();
  if(readIt->getQuality() != NULL) {
    read->qual = readIt->getQuality();
  } else { read->qual = "*";}
}

void *fillLocations(void *args) {
  locate_thread_t params = ((locate_thread_t *)args)[0];
  LocateTags *locate = params.locate;
  uint current_pos = 0;
  int semValue;
  sem_getvalue(params.remaining_seats,&semValue);
  while(((locate_thread_t *)args)[0].keepGoing || semValue > 0) {
    // We need one seat for the current read
    sem_wait(params.remaining_seats);
    sem_getvalue(params.remaining_seats,&semValue);

    if(((locate_thread_t *)args)[0].keepGoing && semValue >= 0) {
      params.seats[current_pos] = locate->locateTag(params.reads[current_pos]);
      sem_post(params.processed_reads);
      // delete the read
      //delete params.reads[current_pos];
      current_pos = (current_pos + 1) % params.total_seats;
    }
  }
  pthread_exit(NULL);
  return NULL;
}
