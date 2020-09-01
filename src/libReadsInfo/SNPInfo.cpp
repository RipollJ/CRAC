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

#include "SNPInfo.h"
#include "libSSA/utils.h"

SNPInfo::SNPInfo(ChrPosition chrPos, uint position,
		 float score, bool is_duplicated, 
		 char actual_nuc, char expected_nuc):
  chrPos(chrPos), score(score),
  is_duplicated(is_duplicated),
  actual_nuc(actual_nuc),
  expected_nuc(expected_nuc) 
{
  setPosition(position);
  if (actual_nuc == 0) {
    if (expected_nuc == 0) {
      type = SNP_UNKNOWN;
    } else {
      type = SNP_DELETION;
    }
  } else {
    if (expected_nuc == 0) {
      type = SNP_INSERTION;
    } else {
      type = SNP_SUBSTITUTION;
    }
  }
}

SNPInfo::~SNPInfo() {
}

char SNPInfo::getActualNucleotide(int strand) {
  if (getSNPType() == SNP_UNKNOWN || getSNPType() == SNP_DELETION)
    return '?';
  return (strand == 1) ? actual_nuc : complementDNA(actual_nuc);;
}

char SNPInfo::getExpectedNucleotide(int strand) {
  if (getSNPType() == SNP_UNKNOWN || getSNPType() == SNP_INSERTION)
    return '?';
  return (strand == 1) ? expected_nuc : complementDNA(expected_nuc);
}

ChrPosition &SNPInfo::getChrPosition() {
  return chrPos;
}

float SNPInfo::getScore() {
  return score;
}

snp_type SNPInfo::getSNPType() {
  return type;
}

bool SNPInfo::isDuplicated() {
  return is_duplicated;
}

void SNPInfo::output(ostream &os) {
  os << *this;
}

void SNPInfo::samOutput(ostream &os, int strand) {
  os << "SNP:" << getScore() << ":" << getPosition(strand) << ":" << getChrPosition()
     << ":" << getExpectedNucleotide(strand) << ":" << getActualNucleotide(strand);
}

cigar_type SNPInfo::getCigarInfo() {
  cigar_type cigar;
  cigar.nb = 1;
  cigar.type = 0;               // Just initializing so that gcc is happy
  if (getSNPType() == SNP_SUBSTITUTION) {
    cigar.type = 'X';
  } else if (getSNPType() == SNP_DELETION) {
    cigar.type = 'D';
  } else if (getSNPType() == SNP_INSERTION) {
    cigar.type = 'I';
  }
  return cigar;
}

ostream &operator<<(ostream &os, SNPInfo& i) {
  if (i.isDuplicated()){
    os << "duplicate ";  
  }else{
    os << "single ";
  }
  os << i.getChrPosition() 
     << " pos_SNV="<<i.getPosition() 
     << " ";
  if (i.getSNPType() == SNP_UNKNOWN)
    os << "?->?";
  else {
    if (i.getSNPType() == SNP_INSERTION)
      os << "?";
    else 
      os << i.getExpectedNucleotide() ;
    
    os << "->";
    
    if (i.getSNPType() == SNP_DELETION)
      os << "?";
    else 
      os << i.getActualNucleotide();
    os << " score=" << i.getScore() ; 
	  
  }
  return os;
}