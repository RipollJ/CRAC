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

#include "Parameters.h"
#include "const.h"

using namespace std;

Parameters::Parameters(uint threshold,
		       uint max_nb_located_occs,
		       uint max_splice_length,
		       uint max_bio_ins_del,
		       uint max_localisation_duplication,
		       uint min_localisation_duplication,
		       float percent_min_unique_repetition,
		       float percent_min_duplicate,
		       uint min_occ_repetition,
		       float percent_support_variation_almost_normal,
		       float p_value_variation_biological,
                       float score_intra_ambiguous,
                       float score_inter_ambiguous,
		       uint max_bases_randomly_matched,
		       uint min_bases_before_break,
		       uint max_bases_retrieved,
		       int max_extension_length,
		       uint support_score_window_length,
                       uint nb_bases_substracted_low_average,
                       uint nb_positions_check_ones,
                       uint nb_tags_info_stored,
                       uint nb_threads,
		       float min_support_no_cover,
		       float max_support_out_no_cover,
		       float max_ambiguous_average_high,
                       float min_perc_ones_inside,
                       float min_ratio_support_fall,
                       float max_nb_overlapping_breaks,
                       float percent_min_single,
                       float percent_min_multiple,
                       bool deep_snp_search,
                       uint number_nucleotides_snp_comparison,
                       bool detailed_sam,
		       uint max_verynice_merge,
		       bool stringent_chimera,
		       uint stringent_chimera_min_break_length,
		       uint stringent_chimera_max_number_of_merges,
		       bool no_ambiguity,
		       bool treat_multiple,
		       bool show_progressbar
		       ):
  threshold(threshold)
  ,max_nb_located_occs(max_nb_located_occs)
  ,max_splice_length(max_splice_length)
  ,max_bio_ins_del(max_bio_ins_del)
  ,max_localisation_duplication(max_localisation_duplication)
  ,min_localisation_duplication(min_localisation_duplication)
  ,percent_min_unique_repetition(percent_min_unique_repetition)
  ,percent_min_duplicate(percent_min_duplicate)
  ,min_occ_repetition(min_occ_repetition)
  ,percent_support_variation_almost_normal(percent_support_variation_almost_normal)
  ,p_value_variation_biological(p_value_variation_biological)
  ,score_intra_ambiguous(score_intra_ambiguous)
  ,score_inter_ambiguous(score_inter_ambiguous)
  ,max_bases_randomly_matched(max_bases_randomly_matched)
  ,min_bases_before_break(min_bases_before_break)
  ,max_bases_retrieved(max_bases_retrieved)
  ,max_extension_length(max_extension_length)
  ,support_score_window_length(support_score_window_length)
  ,nb_bases_substracted_low_average(nb_bases_substracted_low_average)
  ,nb_positions_check_ones(nb_positions_check_ones)
  ,nb_tags_info_stored(nb_tags_info_stored)
  ,nb_threads(nb_threads)
  ,min_support_no_cover(min_support_no_cover)
  ,max_support_out_no_cover(max_support_out_no_cover)
  ,max_ambiguous_average_high(max_ambiguous_average_high)
  ,min_perc_ones_inside(min_perc_ones_inside)
  ,min_ratio_support_fall(min_ratio_support_fall)
  ,max_nb_overlapping_breaks(max_nb_overlapping_breaks)
  ,percent_min_single(percent_min_single)
  ,percent_min_multiple(percent_min_multiple)
  ,deep_snp_search(deep_snp_search)
  ,number_nucleotides_snp_comparison(number_nucleotides_snp_comparison)
  ,detailed_sam(detailed_sam)
  ,max_verynice_merge(max_verynice_merge)
  ,stringent_chimera(stringent_chimera)
  ,stringent_chimera_max_number_of_merges(stringent_chimera_max_number_of_merges)
  ,no_ambiguity(no_ambiguity)
  ,treat_multiple(treat_multiple)
  ,show_progressbar(show_progressbar)
{
  if (max_nb_located_occs < min_occ_repetition){
    cerr << "max_nb_located_occs must be greater than min_loc_repetition" << endl << endl;
    exit(2);
  }
}

void Parameters::init() {
  min_break_length = (uint) (MIN_BREAK_LENGTH*threshold);
  stringent_chimera_min_break_length = (uint) (STRINGENT_CHIMERA_MIN_BREAK_LENGTH*threshold);
}
