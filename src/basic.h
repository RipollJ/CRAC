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

#ifndef BASICSINCLUDED
#define BASICSINCLUDED

#include <sys/types.h>
#include <sys/resource.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define mask31 0x0000001F

#define max(x,y) ((x)>(y)?(x):(y))
#define min(x,y) ((x)<(y)?(x):(y))




/*numero de bits del entero de la maquina*/
#define W 32
/* W-1 */
#define Wminusone 31
/*numero de bits del entero de la maquina*/
#define WW 64
/*bits para hacer la mascara para contar mas rapido*/
#define bitsM 8
/*bytes que hacen una palabra */
#define BW 4
#define size_uchar 256


/* bits needed to represent a number between 0 and n */
inline uint bits (uint n){
  uint b = 0;
  while (n) { b++; n >>= 1; }
  return b;
}

/* reads bit p from e */
#define bitget(e,p) ((((e)[(p)/W] >> ((p)%W))) & 1)
/* sets bit p in e */
#define bitset(e,p) ((e)[(p)/W] |= (1<<((p)%W)))
/* cleans bit p in e */
#define bitclean(e,p) ((e)[(p)/W] &= ~(1<<((p)%W)))

/* numero de enteros necesarios para representar e elementos de largo n */
#define enteros(e,n) ((e)*(n))/W+(((e)*(n))%W > 0)

inline ulong GetField(ulong *A, register  ulong len, register ulong index) {
  register ulong i=index*len/W, j=index*len-W*i, result;
  if (j+len <= W)
    result = (A[i] << (W-j-len)) >> (W-len);
  else {
    result = A[i] >> j;
    result = result | (A[i+1] << (WW-j-len)) >> (W-len);
  }
  return result;
}

inline void SetField(ulong *A,register ulong len, register ulong index,register  ulong x) {
   ulong i=index*len/W, j=index*len-i*W;
   ulong mask = ((j+len) < W ? ~0u << (j+len) : 0) | ((W-j) < W ? ~0u >> (W-j) : 0);
   A[i] = (A[i] & mask) | x << j;
   if (j+len>W) {
      mask = ((~0u) << (len+j-W));
      A[i+1] = (A[i+1] & mask)| x >> (W-j);
   }
}

inline unsigned GetFieldW32(unsigned *A,register unsigned index) {
  return A[index];
}

inline void SetField32(unsigned *A, register unsigned index,register unsigned x) {
  A[index]=x;
}

inline unsigned GetFieldW16(unsigned *A,register unsigned index) {
  register unsigned i=index/2, j=(index&1)<<4, result;
  result = (A[i] << (16-j)) >> (16);
  return result;
}

inline unsigned GetFieldW4(unsigned *A,register unsigned index) {
  register unsigned i=index/8, j=(index&0x7)<<2;
  /*register unsigned i=index/8, j=index*4-32*i; */
  return (A[i] << (28-j)) >> (28);
}



#endif

