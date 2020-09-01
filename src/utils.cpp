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

#include "utils.h"

#include <iostream>
#include <string>
#include <stdlib.h>
#include <sstream>
#include <gzstream.h>
#include <cstring>

using namespace std;


uintSA DNAtoInt(char *dna, uintSA dna_length){
  uintSA dna_int = 0;
  for (uintSA i=0 ; i< dna_length ; i++){
    dna_int <<= 2;      
    dna_int |= convNuc(dna[i]);
  }
  return dna_int;
}

uint convNuc(char nuc){
   switch (nuc){
   case 'a' : case 'A' : return 0 ;
   case 'c' : case 'C' : return 1 ;
   case 'g' : case 'G' : return 2 ;
   case 't' : case 'T' : return 3 ;
   default : cerr << "invalid nucleotide "<< nuc <<endl; exit(4);
   }  
   return 0;
}

char intToNuc(uint c) {
  switch(c) {
  case 0: return 'A';
  case 1: return 'C';
  case 2: return 'G';
  case 3: return 'T';
  default: cerr << "Invalid number (0-3): " << c <<endl; exit(5);
  }
}

void intToDNA(uint64 code, uint dna_length, char *dna) {
  uint64 mask = 3;
  for (uint i=0; i < dna_length; i++) {
    dna[dna_length-i-1] = intToNuc(code & mask);
    code >>=2;
  }
}

uint64 factorsToInt(uintSA pos, uchar *dna, uint length) {
  uintSA posFirstFactor = 4 - (pos % 4);
  uintSA numFirstFactor = pos / 4;
  uintSA posLastFactor = (length - posFirstFactor) % 4;
  uint64 factor = dna[numFirstFactor] & ((1 << (posFirstFactor*2)) - 1);

  for(uintSA i = 1 ; i <= (length-posFirstFactor)/4 ; i++){
    factor <<= 8;
    factor |= dna[numFirstFactor+i];
  }
  factor <<= posLastFactor*2;

  factor |= dna[numFirstFactor+(length-posFirstFactor)/4+1] >> (8-2*posLastFactor) & ((1 << (posLastFactor*2)) - 1);
  
  return factor;
}

int comparUint(const void *a1, const void* a2){
  if ((* (uint *) a1) < (* (uint *) a2)){
    return -1;
  }else{
    return ((* (uint *) a1) == (* (uint *) a2)) ? 0 : 1; 
  }
}

uint posBitInMask(uint mask) {
  uint i = 1;
  uint pos = 0;
  if (mask == 0)
    return sizeof(uint)*8;
  while ((mask & i) == 0) {
    i <<= 1;
    pos++;
  }
  return pos;
}

float getChrono() {
	struct timeval time;
	time_t sec;
	suseconds_t usec;

	if (gettimeofday(&time,NULL) != -1) {
		sec = time.tv_sec;
		usec = time.tv_usec;
		return (float)sec*1000000+usec;
	}
	return 0;
}

string intToString(uint i) {
  ostringstream oss;
  oss << i;
  return oss.str();
}

ostream *create_stream(bool gzipped, char *name) {
  if (gzipped) {
    uint len = strlen(name);
    char *gzName = new char[len+4];
    ostream *result;
    sprintf(gzName, "%s.gz", name);
    result = new ogzstream((const char *)gzName, ofstream::out);
    delete [] gzName;
    return result;
  }
  return new ofstream((const char *)name, ofstream::out);
}
