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

#include "Cigar.h"
#include <cassert>

Cigar::Cigar():nb_M(0), nb_I(0), nb_D(0), nb_S(0), nb_P(0), nb_H(0), nb_EQ(0), nb_X(0), nb_N(0), nb_chimeras(0){}


void Cigar::append(uint nb, char op) {
  assert(op == 'M' || op == 'I' || op == 'D' || op == 'S' || op == '='
         || op == 'P' || op == 'H' || op == 'X' || op == 'N');
  cigar_type ct = {nb, op};
  append(ct);
}

void Cigar::append(cigar_type cigar) {
  update_counters(cigar);
  this->cigar.push_back(cigar);
}

void Cigar::filter() {
  char old_type = 0;
  uint old_nb = 0;
  vector<cigar_type>::iterator it = cigar.begin();
  while (it < cigar.end()) {
    // Two identical consecutive types
    if (it->type == old_type 
        && (it->type != 'N' || (it->nb > 0 && old_nb > 0))) {
      it->nb += old_nb;
      it = cigar.erase(it-1);
    } else if (it->type != 'N' && it->nb == 0) {
      it = cigar.erase(it)-1;
    } 

    old_type = it->type;
    old_nb = it->nb;
    it++;
  }
}

void Cigar::prepend(uint nb, char op) {
  assert(op == 'M' || op == 'I' || op == 'D' || op == 'S' || op == '='
         || op == 'P' || op == 'H' || op == 'X' || op == 'N');
  cigar_type ct = {nb, op};
  prepend(ct);
}

void Cigar::prepend(cigar_type cigar) {
  update_counters(cigar);
  this->cigar.insert(this->cigar.begin(), cigar);
}

  /* Queries */

uint Cigar::count() const {
  return cigar.size();
}

const cigar_type& Cigar::get(uint i) const{
  return cigar[i];
}

uint Cigar::getReferenceAlignementLength() const {
  return nb_M + nb_D + nb_N + nb_EQ + nb_X;
}

uint Cigar::getNbDeletions() const {
  return nb_D + nb_N;
}

uint Cigar::getNbPositionsRead() const {
  return nb_M + nb_EQ + nb_X + nb_I + nb_S;
}

uint Cigar::getEditDistance() const {
  return nb_X + nb_I + nb_D + nb_N;
}

uint Cigar::getNbChimeras() const {
  return nb_chimeras;
}

void Cigar::reverse() {
  vector<cigar_type> reverse_cigar;
  for(int i = cigar.size(); i > 0; i--) {
    reverse_cigar.push_back(cigar[i-1]);
  }
  cigar = reverse_cigar;
}

cigar_type &Cigar::operator[](uint i) {
  return cigar[i];
}


ostream &operator<<(ostream &os, const Cigar &c) {
  if (c.count() > 0) {
    for (uint i=0; i < c.count(); i++) {
      os << c.get(i);
    }
  } else {
    os << "*";
  }
  return os;
}

void Cigar::update_counters(const cigar_type &cigar) {
  if(cigar.type == 'M') {
    nb_M += cigar.nb;
  } else if(cigar.type == '=') {
    nb_EQ += cigar.nb;
  } else if(cigar.type == 'X') {
    nb_X += cigar.nb;
  } else if(cigar.type == 'I') {
    nb_I += cigar.nb;
  } else if(cigar.type == 'D') {
    nb_D += cigar.nb;
  } else if(cigar.type == 'N') {
    if(cigar.nb == 0) {
      nb_chimeras += 1;
    } else {
      nb_N += cigar.nb;
    }
  } else if(cigar.type == 'S') {
    nb_S += cigar.nb;
  } else if(cigar.type == 'H') {
    nb_H += cigar.nb;
  }
}
