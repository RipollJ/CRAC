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

#include "bitrankw32int.h"
#include "assert.h"
#include "math.h"
#include <sys/types.h>
#include <cstdlib>

/////////////
//Rank(B,i)// 
/////////////
//_factor = 0  => s=W*lgn
//_factor = P  => s=W*P
//Is interesting to notice
//factor=2 => overhead 50%
//factor=3 => overhead 33%
//factor=4 => overhead 25%
//factor=20=> overhead 5%

BitRankW32Int::BitRankW32Int( ulong *bitarray, ulong _n, bool owner, ulong _factor){
  data=bitarray;
  this->owner = owner;
  this->n=_n;
	ulong lgn=bits(n-1);
	this->factor=_factor;
	if (_factor==0) this->factor=lgn;
	else this->factor=_factor;
	b=32;
	s=b*this->factor;
  integers = n/W+1;
 	BuildRank();
}

BitRankW32Int::~BitRankW32Int() {
	delete [] Rs;
	if (owner) free(data);
}

//Metodo que realiza la busqueda d
void BitRankW32Int::BuildRank(){
	ulong num_sblock = n/s;
	Rs = new ulong[num_sblock+1];// +1 pues sumo la pos cero
	for(ulong i=0;i<num_sblock+1;i++)
	  Rs[i]=0;
	ulong j;
	Rs[0]=0;
	for (j=1;j<=num_sblock;j++) {
		Rs[j]=Rs[j-1];
		Rs[j]+=BuildRankSub((j-1)*factor,factor);
	}
}

ulong BitRankW32Int::BuildRankSub(ulong ini,ulong bloques){
	ulong rank=0,aux;
	for(ulong i=ini;i<ini+bloques;i++) {
		if (i < integers) {
			aux=data[i];
			rank+=popcount32(aux);
		}
	}
	return rank; //retorna el numero de 1's del intervalo

}

ulong BitRankW32Int::rank(ulong i) {
  ++i;
  ulong resp=Rs[i/s];
  ulong aux=(i/s)*factor;
  for (ulong a=aux;a<i/W;a++)
    resp+=popcount32(data[a]);
  resp+=popcount32(data[i/W]  & ((1<<(i & mask31))-1));
  return resp;
}

bool BitRankW32Int::IsBitSet(ulong i) {
  return (1u << (i % W)) & data[i/W];
}

int BitRankW32Int::save(FILE *f) {
  if (f == NULL) return 20;
  if (fwrite (&n,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (&factor,sizeof(ulong),1,f) != 1) return 21;
  if (fwrite (data,sizeof(ulong),n/W+1,f) != n/W+1) return 21;
  if (fwrite (Rs,sizeof(ulong),n/s+1,f) != n/s+1) return 21;
  return 0;
}

int BitRankW32Int::load(FILE *f) {
  if (f == NULL) return 23;
  if (fread (&n,sizeof(ulong),1,f) != 1) return 25;
  b=32; // b is a word
  if (fread (&factor,sizeof(ulong),1,f) != 1) return 25;
  s=b*factor;
  ulong aux=(n+1)%W;
  if (aux != 0)
    integers = (n+1)/W+1;
  else
    integers = (n+1)/W;
  data= (ulong *)calloc(n/W+1, sizeof(ulong));
  if (!data) return 1;
  if (fread (data,sizeof(ulong),n/W+1,f) != n/W+1) return 25;
  this->owner = true;
  Rs= new ulong[n/s+1];
  if (!Rs) return 1;
  if (fread (Rs,sizeof(ulong),n/s+1,f) != n/s+1) return 25;
  return 0;
}

BitRankW32Int::BitRankW32Int(FILE *f, int *error) {
  *error = BitRankW32Int::load(f);
}

ulong BitRankW32Int::SpaceRequirementInBits() {
  return (owner?n:0)+(n/s)*sizeof(ulong)*8 +sizeof(BitRankW32Int)*8; 
}

ulong BitRankW32Int::prev(ulong start) {
      // returns the position of the previous 1 bit before and including start.
      // tuned to 32 bit machine

      ulong i = start >> 5;
      int offset = (start % W);
      ulong answer = start;
      ulong val = data[i] << (Wminusone-offset);

      if (!val) { val = data[--i]; answer -= 1+offset; }

      while (!val) { val = data[--i]; answer -= W; }

      if (!(val & 0xFFFF0000)) { val <<= 16; answer -= 16; }
      if (!(val & 0xFF000000)) { val <<= 8; answer -= 8; }

      while (!(val & 0x80000000)) { val <<= 1; answer--; }
      return answer;
}

ulong BitRankW32Int::select(ulong x) {
  // returns i such that x=rank(i) && rank(i-1)<x or n if that i not exist
  // first binary search over first level rank structure
  // then sequential search using popcount over a int
  // then sequential search using popcount over a char
  // then sequential search bit a bit

  //binary search over first level rank structure
  ulong l=0, r=n/s;
  ulong mid=(l+r)/2;
  ulong rankmid = Rs[mid];
  while (l<=r) {
    if (rankmid<x)
      l = mid+1;
    else
      r = mid-1;
    mid = (l+r)/2;
    rankmid = Rs[mid];
  }
  //sequential search using popcount over a int
  ulong left;
  left=mid*factor;
  x-=rankmid;
        ulong j=data[left];
        ulong ones = popcount32(j);
        while (ones < x) {
    x-=ones;left++;
    if (left > integers) return n;
          j = data[left];
      ones = popcount32(j);
        }
  //sequential search using popcount over a char
  left=left*b;
  rankmid = popcount8(j);
  if (rankmid < x) {
    j=j>>8;
    x-=rankmid;
    left+=8;
    rankmid = popcount8(j);
    if (rankmid < x) {
      j=j>>8;
      x-=rankmid;
      left+=8;
      rankmid = popcount8(j);
      if (rankmid < x) {
        j=j>>8;
        x-=rankmid;
        left+=8;
      }
    }
  }

  // then sequential search bit a bit
        while (x>0) {
    if  (j&1) x--;
    j=j>>1;
    left++;
  }
  return left-1;
}
