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

#ifndef CANDIDAT_BREAK_H
#define CANDIDAT_BREAK_H

#include <utility>
#include "Parameters.h"
#include "libSSA/locateOnGenome.h"
#include "libSSA/chrPosition.h"
#include "types.h"

class SupportBreak;

class CandidatBreak {

  uint pos_start_break;		/* Position in the support of the start of
				   the break.*/
  uint loc_start_break;           /* the last absolute location before the 
				   break */

  int strand_start_break;       /* the strand of the start_break position */

  uint pos_end_break;		/* Position in the support of the end of 
				   the break. */
  
  uint loc_end_break;   /* the first absolute location after the 
			   break */
  

  int strand_end_break;         /* the strand  of the end_break position*/

  SupportBreak* supportBreak; 	    /* The supportBreak we are working on. */

  Parameters* parameters;	/* Parameters chosen by the user or chosen 
				   by default by the great developers. */

  bool is_duplicated;             /* true iff there are at least two possible causes for a same break */

  bool is_single_consistent;             /* true if the candidat locations follow the single location */

  long long int single_distance;           /* distance between (loc_start_break || loc_end_break) and loc_single in case of is_single_consistent */

 public:

  CandidatBreak(uint pos_start, uint loc_start, int strand_start, 
		uint pos_end, uint loc_end, int strand_end,
		SupportBreak *sb, Parameters *p);


  CandidatBreak(const CandidatBreak &c); 

  ~CandidatBreak();

  /**
   * Check the correspondence between Candidat and the single Location.
   * - same strand
   * - same chr
   * - loc_single <= loc_start_break if loc_single is before break
   *   or loc_single => loc_end_break if loc_single is after break  
   */
  void checkSingleCorrespondence(uint pos_single, uint loc_single, int strand_single);

  /**
   * @return chr for the pos_end_break
   */
  const uchar *getChrEndBreak();

  /**
   * @return a ChrPosition for the pos_end_break
   */
  ChrPosition *getChrPositionEndBreak();

  /**
   * @return a ChrPosition for the pos_start_break
   */
  ChrPosition *getChrPositionStartBreak();

  /**
   * @return chr id for the pos_end_break
   */
  ulong getChrIdEndBreak();

  /**
   * @return chr id for the pos_start_break
   */
  ulong getChrIdStartBreak();

  /**
   * @return chr for the pos_start_break
   */
  const uchar *getChrStartBreak();
  
  /**
   * @return the reference genome.
   */
  LocateOnGenome *getGenome();

  /**
   * @return the distance between loc_start_break and loc_end_break
   *         on the reference genome.
   */
  gap_size_t getGenomeGapLength();
  
  /**
   * @return the support length: 
   *         supportBreak->getLength()
   */
  uint getLength();

  /**
   * @return loc_end_break
   */
  uint getLocEndBreak();

  /**
   * @return loc_start_break
   */
  uint getLocStartBreak();

  /**
   * @return pos_end_break
   */
  uint getPosEndBreak();

  /**
   * @return pos_start_break
   */
  uint getPosStartBreak();

  /**
   * @return the difference between pos_start_break and pos_end_break
   *         on the read (the break length).
   */
  uint getReadBreakLength();

  /**
   * @return single_distance;
   */
  long long int getSingleDistance();

  /**
   * @return strand_end_break
   */
  int getStrandEndBreak();

  /**
   * @return strand_start_break
   */
  int getStrandStartBreak();

  /**
   * @return true if pos_end_break == (getReadLength() -1)
   */
  bool hasNoEndBreak();

  /**
   * @return true if pos_start_break == 0
   */
  bool hasNoStartBreak();
  
  /**
   * @return true if the candidat is duplicated (two similar candidats 
   *         for a same SupportBreak),
   *         false otherwise.
   */
  bool isDuplicated();

  /**
   * @return true is the candidat is a SNV, an indel or a splice,
   *         false is the candidat looks like a chimera. 
   */
  bool isNiceCandidat();

  /**
   * @return true is loc_start_break and loc_end_break are on the same chr,
   *         false otherwise.
   */
  bool isSameChr();
  
  /**
   * @return true is loc_start_break and loc_end_break are on the same strand,
   *         false otherwise.
   */
  bool isSameStrand();

  /**
   * @return true is the candidat is single consistent 
   *                         (see checkSingleCorrespondence procedure)
   *         false otherwise.
   */
  bool isSingleConsistent();
 

  /**
   * isDuplicated = flag.
   */
  void setDuplicated(bool flag);
  
  /**
   * loc_end_break = loc_end.
   */
  void setLocEndBreak(uint loc_end);
  
  /**
   * loc_start_break = loc_start.
   */
  void setLocStartBreak(uint loc_start);

  /**
   * pos_end_break = pos_end.
   */
  void setPosEndBreak(uint pos_end);
  
  /**
   * pos_start_break = pos_start.
   */
  void setPosStartBreak(uint pos_start);
  
  /**
   * isSingleConsistent = flag.
   */
  void setSingleConsistent(bool flag);

  /**
   * strand_end_break = flag.
   */
  void setStrandEndBreak(int flag);

  /**
   * strand_start_break = flag.
   */
  void setStrandStartBreak(int flag);

};
#endif
