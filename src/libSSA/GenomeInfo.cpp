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

#include "GenomeInfo.h"
#include <cassert>
#include <fstream>

  /**
   * Array containing the sequence names
   */
  uchar **chr_names;
  /**
   * Number of sequences
   */
  ulong nb_chr;
  /**
   * Array containing the length (in nt) of each sequence
   */
  ulong *length_chr;
  /**
   * Array containing the cumulative length of each sequence.
   * First cell's value is 0
   */
  ulong *total_length;


  /**
   * Read the file whose name is stored in <basename>.conf.
   * Retrieves the information from the file and store them 
   * in the class.
   * @param basename: name of the file without the .conf extension
   */
GenomeInfo::GenomeInfo(char *basename) {
  char *genomeConf = new char[strlen(basename)+6];
  fstream file;
  sprintf(genomeConf, "%s.conf", basename);
     
  file.open(genomeConf, ios::in);
  if (! file.is_open()) {
    cerr << "Cannot open chromosomes file: " << genomeConf << endl;
    exit(2);
  }

  delete [] genomeConf;

  // Store the number of sequences
  file >> skipws >> nb_chr;

  // Allocate memory for the arrays
  chr_names = new uchar*[nb_chr];
  total_length = new ulong[nb_chr+1];
  length_chr = new ulong[nb_chr];
  total_length[0] = 0;

  int i = 0;
  file.width(MAX_LENGTH_CHR_NAME);
  while (! file.eof()) {
    chr_names[i] = new uchar[MAX_LENGTH_CHR_NAME];
    file >> skipws >> chr_names[i];
    file >> skipws >> length_chr[i];
    total_length[i+1] = total_length[i]+ length_chr[i];
    i++;
  }
  file.close();
}

GenomeInfo::GenomeInfo(uchar **chr_names, ulong nb_chr, ulong *length_chr):
  chr_names(chr_names), nb_chr(nb_chr), length_chr(length_chr) {
  total_length = new ulong[nb_chr+1];
  total_length[0] = 0;
  for (ulong i = 0; i < nb_chr; i++) {
    total_length[i+1] = total_length[i] + length_chr[i];
  }
}

GenomeInfo::~GenomeInfo() {
  for (uint i = 0; i < nb_chr; i++)
    delete [] chr_names[i];
  delete [] chr_names;
  delete [] total_length;
  delete [] length_chr;
}

uchar *GenomeInfo::getChrName(ulong i) {
  assert(i >= 0 && i < getNbChr());
  
  return chr_names[i];
}

pair<const uchar *,const ulong> GenomeInfo::getChrNameFromPosition(ulong pos) {
  pair<ulong, ulong> tmp_pair = getNumFromPosition(pos);
  uchar *chr_name;

  if (tmp_pair.first >= getNbChr()) {
    chr_name = NULL;
  } else {
    chr_name = getChrName(tmp_pair.first);
  }
  return (pair<const uchar *, const ulong>(chr_name, tmp_pair.second));
}

ulong GenomeInfo::getChrLength(ulong i) {
  assert(i >= 0 && i < getNbChr());
  
  return length_chr[i];
}

ulong GenomeInfo::getGenomeLength() {
  return total_length[getNbChr()];
}

ulong GenomeInfo::getGenomePosition(uchar *chr, ulong position) {
  return getGenomePosition(getNumFromName(chr), position);
}

ulong GenomeInfo::getGenomePosition(ulong id, ulong position) {
  if (id < getNbChr() && position < length_chr[id]) {
    return total_length[id] + position;
  }
  return ~0;
}

ulong GenomeInfo::getNbChr() {
  return nb_chr;
}

ulong GenomeInfo::getNumFromName(uchar *name) {
  ulong num_chr = 0;
  while (num_chr < getNbChr()
         && strcmp((char *)chr_names[num_chr], (char *)name) != 0){
    num_chr++;
  }
  if (num_chr >= getNbChr()) {
    return getNbChr();
  }
  return num_chr;
}

pair<const ulong, const ulong> GenomeInfo::getNumFromPosition(ulong pos) {
  ulong num_chr = 0;
  
  while (num_chr < getNbChr() && pos >= total_length[num_chr+1]) {
    num_chr++;
  }
  return pair<const ulong, const ulong>(num_chr, pos - total_length[num_chr]);
}

bool GenomeInfo::isValidPosition(ulong pos) {
  return pos < getGenomeLength();
}


void GenomeInfo::save(char *name, const char *ext) {
  char *confName = new char[strlen(name) + strlen(ext)];
  ofstream conf;
  sprintf(confName, "%s%s", name, ext);
  conf.open(confName, ios_base::out | ios_base::binary);
  if (! conf.is_open()) {
    cerr << "Can't write file " << confName << endl;
    exit(3);
  }
  delete [] confName;
  
  // Writing configuration file
  conf << *this;
  conf.close();

}

ostream &operator<<(ostream &os, GenomeInfo &g) {
  os << g.getNbChr();
  for (uint i=0; i < g.getNbChr(); i++) {
    os << endl << g.getChrName(i) 
       << endl << g.getChrLength(i);
  }
  return os;
}
