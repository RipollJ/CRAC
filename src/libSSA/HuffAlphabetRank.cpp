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

#include "HuffAlphabetRank.h"
#include <cstdlib>

   THuffAlphabetRank::THuffAlphabetRank(uchar *s, ulong n, TCodeEntry *codetable, ulong level, ulong factor) {
      left = NULL;
      right = NULL;
      bitrank = NULL;
//       parent = NULL;
      ch = s[0];
      this->codetable = codetable;
      ulong *Binbits = (ulong *)calloc(n/W+1, sizeof(ulong));;
      ulong sum=0,i;

      for (i=0;i<n;i++)
         if ((codetable[s[i]].code & (1u << level)) )  {
           SetField(Binbits,1,i,true);
           sum++;
         } else 
           SetField(Binbits,1,i,false);
      leaf=true;
      for (i=1;i<n;i++)
        if (s[i] != ch)
          leaf = false;
      if (leaf) {
        free(Binbits);
        delete [] s;
        return;
      }
      
      uchar *sfirst=NULL, *ssecond=NULL;
      if (n-sum > 0) {
        sfirst = new uchar[n-sum];
      }
      if (sum > 0) {
        ssecond = new uchar[sum];
      }
      ulong j=0,k=0;
      for (i=0;i<n;i++)
        if (GetField(Binbits,1,i)) ssecond[k++] = s[i];
        else sfirst[j++] = s[i];

      delete [] s;
      bitrank = new BitRankW32Int(Binbits,n,true,factor);
      if (j > 0) {
        left = new THuffAlphabetRank(sfirst,j,codetable,level+1,factor);
      } else if (sfirst) {
        delete [] sfirst;
      }
      if (k>0) {
        right = new THuffAlphabetRank(ssecond,k,codetable,level+1,factor);
      } else if (ssecond) {
        delete [] ssecond;
      }
//       if (level == 0)
//         init_leaves();
   }

// void THuffAlphabetRank::init_leaves() {
//   leaves = new THuffAlphabetRank[size_uchar];
//   for (uchar c = 0; c < size_uchar; c++) {
//     leaves[c] = getLeaf(c);
//   }
// }

THuffAlphabetRank *THuffAlphabetRank::getLeaf(uchar c) {
  if (codetable[c].count == 0) return NULL;
  ulong level = 0;
  ulong code = codetable[c].code;
  THuffAlphabetRank *currentNode = this;
  while (! currentNode->leaf) {
    if ((code & (1u<<level)) == 0) {
      currentNode = currentNode->left;
    } else {
      currentNode = currentNode->right;
    }
    ++level;
  }
  return currentNode;
}

   ulong THuffAlphabetRank::rank(uchar c, ulong i) { // returns the number of characters c before and including position i
      THuffAlphabetRank *temp=this;
      if (codetable[c].count == 0) return 0;
      ulong level = 0;
      ulong code = codetable[c].code;
      while (!temp->leaf) {
        if ((code & (1u<<level)) == 0) {
          i = i-temp->bitrank->rank(i);
          temp = temp->left;
        }
        else {
          i = temp->bitrank->rank(i)-1;
          temp = temp->right;
        }
        ++level;
      }
      return i+1;
   }

// ulong THuffAlphabetRank::select(uchar c, ulong i) {
//   THuffAlphabetRank *node = leaves[c];

//   if (node->parent == NULL)
//     return i;
  
//   bool bit = node->parent->right == node;

//   while(node->parent != NULL) {
//     node = node->parent;
//     node->bitrank->select
//   }
// }


   bool THuffAlphabetRank::IsCharAtPos(uchar c, ulong i) {
      THuffAlphabetRank *temp=this;
      if (codetable[c].count == 0) return false;
      ulong level = 0;
      ulong code = codetable[c].code;
      while (!temp->leaf) {
         if ((code & (1u<<level))==0) {
            if (temp->bitrank->IsBitSet(i)) return false;
            i = i-temp->bitrank->rank(i);
            temp = temp->left;
         }
         else {
            if (!temp->bitrank->IsBitSet(i)) return false;
            i = temp->bitrank->rank(i)-1;
            temp = temp->right;
         }
         ++level;
      }
      return true;
   }
   uchar THuffAlphabetRank::charAtPos(ulong i) {
      THuffAlphabetRank *temp=this;
      while (!temp->leaf) {
        if (temp->bitrank->IsBitSet(i)) {
          i = temp->bitrank->rank(i)-1;
          temp = temp->right;
        } else {
          i = i-temp->bitrank->rank(i);
          temp = temp->left;
        }
      }
      return temp->ch;
   }
   uchar THuffAlphabetRank::charAtPos2(ulong i, ulong *rank) {
     THuffAlphabetRank *temp=this;
     while (!temp->leaf) {
       if (temp->bitrank->IsBitSet(i)) {
         i = temp->bitrank->rank(i)-1;
         temp = temp->right;
       } else {
         i = i-temp->bitrank->rank(i);
         temp = temp->left;
       }
     }
     (*rank)=i;
     return temp->ch;
   }

   ulong THuffAlphabetRank::SpaceRequirementInBits() {
      ulong bits=sizeof(THuffAlphabetRank)*8;
      if (left!=NULL) bits+=left->SpaceRequirementInBits();
      if (right!=NULL) bits+=right->SpaceRequirementInBits();
      if (bitrank != NULL)
         bits += bitrank->SpaceRequirementInBits();
      return bits;
   }
   bool THuffAlphabetRank::Test(uchar *s, ulong n) {
      // testing that the code works correctly
      uchar C[size_uchar];
      ulong i,j;
      bool correct=true;
      for (j=0;j<size_uchar;j++)
         C[j] = 0;
      for (i=0;i<n;i++) {
         C[s[i]]++;
         if (C[s[i]]!=rank(s[i],i)) {
            correct = false;
            printf("%u (%c): %d<>%u\n",i,s[i],C[s[i]],rank(s[i],i));
         }
      }
      return correct;
   }
   THuffAlphabetRank::~THuffAlphabetRank() {
      if (left!=NULL) delete left;
      if (right!=NULL) delete right;
      if (bitrank!=NULL)
         delete bitrank;
   }

   int THuffAlphabetRank::save(FILE *f) {
     uchar sons;
     if (f == NULL) return 20;
     if (fwrite (&leaf,sizeof(uchar),1,f) != 1) return 21;
     if (fwrite (&ch,sizeof(uchar),1,f) != 1) return 21;
     if (!leaf) {
       if (bitrank->save(f) != 0) return 21;
       sons=0;
       if (left != NULL) sons=sons | 1;
       if (right != NULL) sons=sons | 2;
       if (fwrite (&sons,sizeof(uchar),1,f) != 1) return 21;
       if (left != NULL) 
         if (left->save(f)!=0) return 21;
       if (right != NULL)
         if (right->save(f)!=0) return 21;
     }
     return 0;
   }

   int THuffAlphabetRank::load(FILE *f, TCodeEntry *_codetable) {
     int error;
     uchar sons;
     left = NULL;
     right = NULL;
//      parent = NULL;
     bitrank = NULL;

     if (f == NULL) return 23;
     if (fread (&leaf,sizeof(uchar),1,f) != 1) return 25;
     if (fread (&ch,sizeof(uchar),1,f) != 1) return 25;
     codetable=_codetable;
     if (!leaf){ 
       bitrank= new BitRankW32Int(f,&error); 
       if (error != 0) return error;
       if (fread (&sons,sizeof(uchar),1,f) != 1) return 25;
       if ((sons&1) > 0 ) {
         left = new THuffAlphabetRank(f, _codetable, &error);
//          left->parent = this;
         if (error != 0) return error;
       }
       if ((sons&2) > 0 ) {
         right = new THuffAlphabetRank(f, _codetable, &error);
//          right->parent = this;
         if (error != 0) return error;
       }
     }
     return 0;
   }
   THuffAlphabetRank::THuffAlphabetRank(FILE *f, TCodeEntry *_codetable, int *error) {
     *error = load(f, _codetable);
//      init_leaves();
   }
