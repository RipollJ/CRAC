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

#include "samLine.h"
#include "../libSSA/utils.h"
#include <sstream>
#include <cstring>


SamLine::SamLine() {
  // init SAM fields with default values
  rname = "*";
  qname = "*";
  rnext = "*";
  seq = "*";
  qual = "*";
  flag = 0;
  pos = 0;
  mapQ = 255;
  //mapQ = 0;
  pnext = 0;
  tlen = 0;
  //is_mapped = true;
}

// Query template NAME
void SamLine::setQname(const string& qname) {
  this->qname = qname;
}
const string& SamLine::getQname() const {
  return qname;
}

// bitwize FLAG
void SamLine::setFlag(uint flag){
  this->flag = flag;
}
uint SamLine::getFlag() const {
  return flag;
}

// generic flag bit methods
bool SamLine::isFlagBitSet(uint bit) const {
  return flag & bit;
}
void SamLine::setFlagBit(uint bit) {
  flag = flag | bit;
}
void SamLine::unsetFlagBit(uint bit) {
  flag = flag^bit;
}

// Flag operations
// bit 1
bool SamLine::isTemplateHavingMultipleSegments() const {
  return isFlagBitSet(1);
}
void SamLine::setTemplateHavingMultipleSegments() {
  setFlagBit(1);
}
void SamLine::unsetTemplateHavingMultipleSegments() {
  unsetFlagBit(1);
}

// bit 2
bool SamLine::isEachSegmentsMapped() const {
  return isFlagBitSet(2);
}
void SamLine::setEachSegmentsMapped() {
  setFlagBit(2);
}
void SamLine::unsetEachSegmentsMapped() {
  unsetFlagBit(2);
}

// bit 4
bool SamLine::isSegmentUnmapped() const {
  return isFlagBitSet(4);
}
void SamLine::setSegmentUnmapped() {
  setFlagBit(4);
}
void SamLine::unsetSegmentUnmapped() {
  unsetFlagBit(4);
}

// bit 8
bool SamLine::isNextSegmentUnmapped() const {
  return isFlagBitSet(8);
}
void SamLine::setNextSegmentUnmapped() {
  setFlagBit(8);
}
void SamLine::unsetNextSegmentUnmapped() {
  unsetFlagBit(8);
}

// bit 16
bool SamLine::isSeqReverseComplemented() const {
  return isFlagBitSet(16);
}
void SamLine::setSeqReverseComplemented() {
  setFlagBit(16);
}
void SamLine::unsetSeqReverseComplemented() {
  unsetFlagBit(16);
}

// bit 32
bool SamLine::isNextSeqReverseComplemented() const {
  return isFlagBitSet(32);
}
void SamLine::setNextSeqReverseComplemented() {
  setFlagBit(32);
}
void SamLine::unsetNextSeqReverseComplemented() {
  unsetFlagBit(32);
}

// bit 64
bool SamLine::isFirstSegmentInTheTemplate() const {
  return isFlagBitSet(64);
}
void SamLine::setFirstSegmentInTheTemplate() {
  setFlagBit(64);
}
void SamLine::unsetFirstSegmentInTheTemplate() {
  unsetFlagBit(64);
}

// bit 128
bool SamLine::isLastSegmentInTheTemplate() const {
  return isFlagBitSet(128);
}
void SamLine::setLastSegmentInTheTemplate() {
  setFlagBit(128);
}
void SamLine::unsetLastSegmentInTheTemplate() {
  unsetFlagBit(128);
}

// bit 256
bool SamLine::isSecondaryAlignement() const {
  return isFlagBitSet(256);
}
void SamLine::setSecondaryAlignement() {
  setFlagBit(256);
}
void SamLine::unsetSecondaryAlignement() {
  unsetFlagBit(256);
}

// bit 512
bool SamLine::isFailingQualityControl() const {
  return isFlagBitSet(512);
}
void SamLine::setFailingQualityControl() {
  setFlagBit(512);
}
void SamLine::unsetFailingQualityControl() {
  unsetFlagBit(512);
}

// bit 1024
bool SamLine::isPCRDuplicated() const {
  return isFlagBitSet(1024);
}
void SamLine::setPCRDuplicated() {
  setFlagBit(1024);
}
void SamLine::unsetPCRDuplicated() {
  unsetFlagBit(1024);
}

// bit 2048
bool SamLine::isChimericAlignement() const {
  return isFlagBitSet(2048);
}
void SamLine::setChimericAlignement() {
  setFlagBit(2048);
}
void SamLine::unsetChimericAlignement() {
  unsetFlagBit(2048);
}

// Reference sequence NAME
void SamLine::setRname(const string& rname) {
  this->rname = rname;
}
const string& SamLine::getRname() const {
  return rname;
}

// 1-based leftmost mapping POSition
void SamLine::setPos(uint pos) {
  this->pos = pos;
}
uint SamLine::getPos() const {
  return pos;
}

// MAPping Quality
void SamLine::setMapQ(uint mapQ) {
  this->mapQ = mapQ;
}
uint SamLine::getMapQ() const {
  return mapQ;
}

// CIGAR string
void SamLine::setCigar(const Cigar &cigar) {
  this->cigar = cigar;
}
const Cigar& SamLine::getCigar() const {
  return cigar;
}

// Ref. name of the mate/next segment
void SamLine::setRnext(const string& rnext) {
  this->rnext = rnext;
}
const string& SamLine::getRnext() const{
  return rnext;
}

// Position of the mate/next segment
void SamLine::setPnext(uint pnext) {
  this->pnext = pnext;
}
uint SamLine::getPnext() const{
  return pnext;
}

// observed Template LENgth
void SamLine::setTlen(int tlen) {
  this->tlen = tlen;
}
int SamLine::getTlen() const{
  return tlen;
}

// segment SEQuence
void SamLine::setSeq(const string& seq) {
  this->seq = seq;
}
const string& SamLine::getSeq() const {
  return seq;
}
void SamLine::reverseComplementeSeq() {
  string new_seq;
  for (int i = seq.length() - 1; i >= 0; i--) {
    new_seq.push_back(complementDNA(seq[i]));
  }
  setSeq(new_seq);
}

// ASCII of Phred-scaled base QUALity+33
void SamLine::setQual(const string& qual) {
  this->qual = qual;
}
const string& SamLine::getQual() const {
  return qual;
}
void SamLine::reverseQual() {
  string new_qual;
  for (int i = qual.length() - 1; i >= 0; i--) {
    new_qual.push_back(qual[i]);
  }
  setQual(new_qual);
}

// Optional fields
bool SamLine::addOptionalField(const char *tag, int val) {
  ostringstream convert;
  convert << "i:" << val;
  optionalFields[tag] = convert.str();
  return true;
}
bool SamLine::addOptionalField(const char *tag, const string& val) {
  string new_val = val;
  optionalFields[tag] = new_val.insert(0,"Z:"); 
  return true;
}
bool SamLine::addOptionalField(const char *tag, const char* val) {
  string new_val = val;
  optionalFields[tag] = new_val.insert(0,"Z:"); 
  return true;
}
const string& SamLine::getOptionalField(const char *tag) {
  return optionalFields[tag];
}
bool SamLine::isOptionalFieldDefined(const char *tag) {
  return optionalFields.count(tag) > 0;
}

void SamLine::removeAllUserOptionalFields() {
  map<string,string>::iterator it=optionalFields.begin();
  while(it != optionalFields.end()) {
    if(it->first[0] == 'X' || it->first[0] == 'Y' || it->first[0] == 'Z')
      optionalFields.erase(it++);
    else
      ++it;
  }
}

// write the SAM line in the file
ostream &SamLine::writeLine(ostream &os){
  os << getQname()  << "\t"
     << getFlag()   << "\t"
     << getRname()  << "\t"
     << getPos()    << "\t"
     << getMapQ()   << "\t"
     << getCigar()  << "\t"
     << getRnext()  << "\t"
     << getPnext()  << "\t"
     << getTlen()   << "\t"
     << getSeq()    << "\t"
     << getQual();
  for (map<string,string>::iterator it=optionalFields.begin(); it!=optionalFields.end(); ++it)
        os << "\t" << it->first << ":" << it->second;
  os << endl;
  return os; 
}

