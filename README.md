# CRAC
Forked from [INRIA forge](http://crac.gforge.inria.fr/) for bioconda recipes

-------------------------------------------------------------

crac - Package

CRAC is a tool to analyze High Throughput Sequencing (HTS) data in
comparison to a reference genome. It is intended for transcriptomic
and genomic sequencing reads. More precisely, with transcriptomic
reads as input, it predicts point mutations, indels, splice junction,
and chimeric RNAs (ie, non colinear splice junctions). CRAC can also
output positions and nature of sequence error that it detects in the
reads. CRAC uses a genome index. This index must be computed before
running the read analysis. For this sake, use the command "crac-index"
on your genome files. You can then process the reads using the command
crac. See the man page of CRAC (help file) by typing "man crac". CRAC
requires large amount of main memory on your computer. For processing
against the Human genome, say 50 million reads of 100 nucleotide each,
CRAC requires about 40 gigabytes of main memory. Check whether the
system of your computing server is equipped with sufficient amount of
memory before launching an analysis.  

crac - Source code

CRAC is freely available under a GPL-compliant CeCILL license. The
distribution comes as an archive (a gzipped tarball file .tar.gz)
containing the source code.

crac - Installation

1) Unpack the archive
2) Enter the directory crac-version-number
3) Type ./configure
4) If everything went fine, run make
5) You may want to check everything is ok by running make check
6) Finally, you can install the software on your system by running make install

If the configure step failed, that may be due to a missing library. In particular the zlib is needed. On a Debian or Debian-like system, you'll need to install zlib1g, zlib1g-dev.

After the compilation succeded you can find a more complete documentation 
inside the doc directory and a man page is also available.

-------------------------------------------------------------------------

Copyright © 2010-2013 --
- IRB/INSERM  (Institut de Recherches en Biothérapie / Institut National de la Santé et de la Recherche Médicale)
- LIFL/INRIA  (Laboratoire d'Informatique Fondamentale de Lille / Institut National de Recherche en Informatique et Automatique)
- LIRMM/CNRS   (Laboratoire d'Informatique, de Robotique et de Microélectronique de Montpellier / Centre National de la Recherche Scientifique)
- LITIS        (Laboratoire d'Informatique, du Traitement de l'Information et des Systèmes).

Auteurs/Authors:  
 - Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>
 - Mikaël SALSON    <mikael.salson@lifl.fr>
 - Thérèse COMMES   <commes@univ-montp2.fr>
 - Éric RIVALS      <eric.rivals@lirmm.fr>

Programmeurs/Progammers: 
 - Nicolas PHILIPPE <nicolas.philippe@lirmm.fr>
 - Mikaël SALSON    <mikael.salson@lifl.fr>        

with additional contribution for the packaging of:	                      
 - Alban MANCHERON  <alban.mancheron@lirmm.fr>               
                                                                             
Contact: CRAC list   <crac-bugs@lists.gforge.inria.fr> 
  
 -------------------------------------------------------------------------

 Ce fichier fait partie du programme CRAC.

 Ce logiciel est régi  par la licence CeCILL  soumise au droit français et
 respectant les principes  de diffusion des logiciels libres.  Vous pouvez
 utiliser, modifier et/ou redistribuer ce programme sous les conditions de
 la licence CeCILL  telle que diffusée par le CEA,  le CNRS et l'INRIA sur
 le site "http://www.cecill.info".

 En contrepartie de l'accessibilité au code source et des droits de copie,
 de modification et de redistribution accordés par cette licence, il n'est
 offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
 seule une responsabilité  restreinte pèse  sur l'auteur du programme,  le
 titulaire des droits patrimoniaux et les concédants successifs.

 À  cet égard  l'attention de  l'utilisateur est  attirée sur  les risques
 associés  au chargement,  à  l'utilisation,  à  la modification  et/ou au
 développement  et à la reproduction du  logiciel par  l'utilisateur étant
 donné  sa spécificité  de logiciel libre,  qui peut le rendre  complexe à
 manipuler et qui le réserve donc à des développeurs et des professionnels
 avertis  possédant  des  connaissances  informatiques  approfondies.  Les
 utilisateurs  sont donc  invités  à  charger  et  tester  l'adéquation du
 logiciel  à leurs besoins  dans des conditions  permettant  d'assurer  la
 sécurité de leurs systêmes et ou de leurs données et,  plus généralement,
 à l'utiliser et l'exploiter dans les mêmes conditions de sécurité.

 Le fait  que vous puissiez accéder  à cet en-tête signifie  que vous avez
 pris connaissance  de la licence CeCILL,  et que vous en avez accepté les
 termes.

 -------------------------------------------------------------------------

 This File is part of the CRAC program.

 This software is governed by the CeCILL license under French law and
 abiding by the rules of distribution of free software. You can use,
 modify and/ or redistribute the software under the terms of the CeCILL
 license as circulated by CEA, CNRS and INRIA at the following URL
 "http://www.cecill.info".

 As a counterpart to the access to the source code and rights to copy,
 modify and redistribute granted by the license, users are provided only
 with a limited warranty and the software's author, the holder of the
 economic rights, and the successive licensors have only limited
 liability.

 In this respect, the user's attention is drawn to the risks associated
 with loading, using, modifying and/or developing or reproducing the
 software by the user in light of its specific status of free software,
 that may mean that it is complicated to manipulate, and that also
 therefore means that it is reserved for developers and experienced
 professionals having in-depth computer knowledge. Users are therefore
 encouraged to load and test the software's suitability as regards their
 requirements in conditions enabling the security of their systems and/or
 data to be ensured and, more generally, to use and operate it in the same
 conditions as regards security.

 The fact that you are presently reading this means that you have had
 knowledge of the CeCILL license and that you accept its terms.

-------------------------------------------------------------------------

Les commentaires sont les bienvenus/Comments are welcome.
