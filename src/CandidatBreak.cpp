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

#include "CandidatBreak.h"
#include "SupportBreak.h"
#include <math.h>

using namespace std;

CandidatBreak::CandidatBreak(uint pos_start, uint loc_start, int strand_start, 
			     uint pos_end, uint loc_end, int strand_end,
			     SupportBreak *sb, Parameters *p 
			     ):pos_start_break(pos_start), loc_start_break(loc_start), strand_start_break(strand_start), pos_end_break(pos_end), loc_end_break(loc_end), strand_end_break(strand_end), supportBreak(sb), parameters(p) {
  is_single_consistent = false;
  single_distance = ~0;
  is_duplicated = false;
} 

CandidatBreak::CandidatBreak(const CandidatBreak &c):pos_start_break(c.pos_start_break), loc_start_break(c.loc_start_break), strand_start_break(c.strand_start_break), pos_end_break(c.pos_end_break), loc_end_break(c.loc_end_break), strand_end_break(c.strand_end_break), supportBreak(c.supportBreak), parameters(c.parameters) {
  is_single_consistent = c.is_single_consistent;
  is_duplicated = c.is_duplicated;
  single_distance = c.single_distance;
}

CandidatBreak::~CandidatBreak(){}


void CandidatBreak::checkSingleCorrespondence(uint pos_single, uint loc_single, int strand_single){
  // case of loc_single does not exist
  if (strand_single == 0){
    is_single_consistent = false;
    single_distance = ~0;
  }else{
    ulong chr_single = getGenome()->getIdChromosome(loc_single);
    // if pos_single is before break, checking of strand and chr
    if (strand_single == strand_start_break 
	&& getPosEndBreak() > pos_single
	&&  chr_single == getChrIdStartBreak()){
      // compute the distance between loc_single and loc_start_break
      if ((strand_single == 1 && pos_single <= getPosStartBreak()) 
	  || (strand_single == -1 && pos_single >= getPosStartBreak()))
	single_distance = (long long int) getLocStartBreak() - loc_single;
      else
	single_distance = (long long int) loc_single - getLocStartBreak();
      // is_single consistent iff single_distance >= 0 && single_distance <= parameters->max_splice_length
      if (single_distance >= 0 && single_distance <= parameters->max_splice_length) 
	is_single_consistent = true;
      else{
	is_single_consistent = false;
	single_distance = ~0;
      }
    }
    // if pos_single is after break, checking of strand and chr
    else if (strand_single == strand_end_break 
	     && getPosStartBreak() < pos_single
	     && chr_single == getChrIdEndBreak()){
      // compute the distance between loc_single and loc_start_break
      if ((strand_single == 1 && pos_single >= getPosEndBreak()) 
	  || (strand_single == -1 && pos_single <= getPosEndBreak()))
	single_distance = (long long int) loc_single - getLocEndBreak();
      else
	single_distance = (long long int) getLocEndBreak() - loc_single;
      // is_single consistent iff single_distance >= 0
      if (single_distance >= 0) 
	is_single_consistent = true;
      else{
	is_single_consistent = false;
	single_distance = ~0;
      }
    }
    // other cases: not single consistency
    else{
      is_single_consistent = false;
      single_distance = ~0;
    }
  }
}


const uchar *CandidatBreak::getChrEndBreak(){
  return getGenome()->getChromosome(loc_end_break);
}


ChrPosition *CandidatBreak::getChrPositionEndBreak(){
  return getGenome()->getChrPos(loc_end_break,strand_end_break);
}


ChrPosition *CandidatBreak::getChrPositionStartBreak(){
  return getGenome()->getChrPos(loc_start_break,strand_start_break);
}

ulong CandidatBreak::getChrIdEndBreak(){
  return getGenome()->getIdChromosome(loc_end_break);
}

ulong CandidatBreak::getChrIdStartBreak(){
  return getGenome()->getIdChromosome(loc_start_break);
}

const uchar *CandidatBreak::getChrStartBreak(){
    return getGenome()->getChromosome(loc_start_break);
}

LocateOnGenome *CandidatBreak::getGenome(){
  return supportBreak->getGenome();
}

gap_size_t CandidatBreak::getGenomeGapLength(){
  if (strand_start_break != strand_end_break
      || (getChrIdStartBreak() != getChrIdEndBreak()))
    return min(labs(loc_end_break - loc_start_break),
	       labs(loc_start_break - loc_end_break)) - 1;
  if (strand_end_break == -1)
    return (gap_size_t)loc_start_break - (gap_size_t)loc_end_break-1;
  return (gap_size_t)loc_end_break - (gap_size_t)loc_start_break-1;
}

uint CandidatBreak::getLength(){
  return supportBreak->getLength();
}

uint CandidatBreak::getLocEndBreak(){
  return loc_end_break;
}

uint CandidatBreak::getLocStartBreak(){
  return loc_start_break;
}

uint CandidatBreak::getPosEndBreak(){
  return pos_end_break;
}

uint CandidatBreak::getPosStartBreak(){
  return pos_start_break;
}

uint CandidatBreak::getReadBreakLength(){
  return (pos_end_break - pos_start_break + 1);
}

long long int CandidatBreak::getSingleDistance(){
  return single_distance;
}

int CandidatBreak::getStrandEndBreak(){
  return strand_end_break;
}

int CandidatBreak::getStrandStartBreak(){
  return strand_start_break;
}

bool CandidatBreak::hasNoEndBreak(){
  return (pos_end_break == (getLength() - 1)); 
}

bool CandidatBreak::hasNoStartBreak(){
  return (pos_start_break == 0);
}

bool CandidatBreak::isDuplicated(){
  return is_duplicated;
}

bool CandidatBreak::isNiceCandidat(){
  return (isSameChr() && isSameStrand()
	  && (
	      (getStrandStartBreak() == 1 && (getLocStartBreak() + supportBreak->getThreshold() - 1) < getLocEndBreak())
	      || (getStrandStartBreak() == -1 && (getLocEndBreak() + supportBreak->getThreshold() - 1) < getLocStartBreak())
	      )
	  && (getGenomeGapLength() <= parameters->max_splice_length)
	  );
}

bool CandidatBreak::isSameChr(){
  if (getChrIdStartBreak() != getChrIdEndBreak())
    return false;
  else
    return true;
}

bool CandidatBreak::isSameStrand(){
  if (strand_start_break != strand_end_break)
    return false;
  else
    return true;
}

bool CandidatBreak::isSingleConsistent(){
  return is_single_consistent;
}
 
void CandidatBreak::setDuplicated(bool flag){
  is_duplicated = flag;
}

void CandidatBreak::setLocEndBreak(uint loc_end){
  loc_end_break = loc_end;  
}

void CandidatBreak::setLocStartBreak(uint loc_start){
  loc_start_break = loc_start;
}

void CandidatBreak::setPosEndBreak(uint pos_end){
  pos_end_break = pos_end;
}

void CandidatBreak::setPosStartBreak(uint pos_start){
  pos_start_break = pos_start;
}

void CandidatBreak::setSingleConsistent(bool flag){
  is_single_consistent = flag;
}

void CandidatBreak::setStrandEndBreak(int flag){
  strand_end_break = flag;
}

void CandidatBreak::setStrandStartBreak(int flag){
  strand_start_break = flag;
}
