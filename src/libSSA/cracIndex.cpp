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

#include <iostream>
#include <fstream>
#include <getopt.h>
#include <list>
#include <cstring>
#include <cstdlib>
#include <ctime>
using namespace std;


extern "C" {
#include <sys/types.h>
#include <unistd.h>
}

#include <blockwise_sa.h>
#include <SSA.h>
#include <config.h>
#include <seqan/sequence.h>
#include <seqan/file.h>

#include <cracIndex.h>
#include <GenomeInfo.h>
#include <utils.h>

typedef String<Dna, Packed<> > TStr;

using namespace seqan;

int verbose = 0;

void usage(char *prog) {
  cerr << endl << prog << " version " << PACKAGE_VERSION;
  cerr << "\tCompiled on " <<__DATE__ << endl << endl;
  cerr << "Usage : "<<prog << " [options] <command> <output file> <input file>+"
       << endl << endl 
       << "  command must be one of:" << endl
       << "    index: create an index on the specified input file(s)." 
       << endl << endl
       << "    get:   get a (multi)FASTA file containing the original sequences."
       << endl << endl
       << "  options can be (for the index command only):" << endl << endl
       << "  -b <bucket_size>\t the size of the bucket for the index construction"
       << endl 
       << "                  \t (default " << THIS_BUCKET_SIZE << ")" << endl
       << "  -d <diff-cover> \t parameter for the index construction (default "
       << COVER << ")"
       << endl 
       << "  -v              \t verbose mode"
       << endl << endl 
       << "  Examples: " << endl
       << "\tIndexing:" << endl
       << "\t\t" << prog << " index myIndex sequence1.fa sequence2.fa sequence3.fa"
       << endl
       << "\t\t\tYou can specify FASTA or MultiFASTA file(s)." << endl
       << "\t\t\tIn this example, two files will be created:" << endl
       << "\t\t\t- myIndex.ssa (index storing the compressed sequences)" << endl
       << "\t\t\t- myIndex.conf (information on sequence names and length)" << endl
       << "\tExtracting:" << endl
       << "\t\t" << prog << " get sequences.fa myIndex" << endl
       << "\t\t\tSequences indexed in myIndex will be extracted" << endl
       << "\t\t\tto the sequences.fa file" << endl;
  exit(1);
}


void index(char **input_files, int nb_input_files, char *output_file, 
           uint cover, uint bucket_size) {
  unsigned long long n = 0;
  void *index;
  ifstream file;
  ofstream conf;
  TStr *text = new TStr();
  uchar *bwt;
  int i = 0;
  list<pair<uint, String<char> > > sequences;

  for (i = 0; i < nb_input_files; i++) {
    if (verbose) 
      cerr << "Reading file " << input_files[i] << endl;

    file.open(input_files[i], ios_base::in | ios_base::binary);
    if (file.is_open()) {
      String<Dna5> current_seq;
      String<char> fasta_tag;
      uint current_length;
      do {
        readMeta(file, fasta_tag, Fasta());
        read(file, current_seq, Fasta());
        current_length = length(current_seq);
        if (verbose && current_length)
          cerr << "\tSequence of length " << current_length << endl;
        if (current_length) {
          // Search sequence name in fasta_tag
          Iterator<String<char> >::Type it = begin(fasta_tag);
          while (it != end(fasta_tag) && value(it) != ' ')
            it++;
          fasta_tag = prefix(fasta_tag, it);

          if (length(fasta_tag) > MAX_LENGTH_CHR_NAME) {
            cerr << "Warning the name of the sequence is too long and will be cut:" << endl
                 << "\t" << fasta_tag << endl;
            fasta_tag = prefix(fasta_tag, MAX_LENGTH_CHR_NAME-1);
          }
          
          sequences.push_back(pair<uint, String<char> >(current_length,
                                                        fasta_tag));
          
          uint this_seq_length = length(current_seq);
          for (uint i = 0; i < this_seq_length; i++) {
            // Replace N with a random nucleotide... that's not great. I know.
            // But the text is encoded on two bits.
            if (current_seq[i] == 'N') {
              current_seq[i] = randomDNA();
            }
          }

          n += current_length;
          *text += current_seq;
        }
      } while (current_length > 0);
      file.close();
    } else {
      cerr<<"Unable to open the input file " << input_files[i] <<  endl;
      exit(2);
    }
  }


  // Creating configuration file
  char *confName = new char[strlen(output_file)+6];
  sprintf(confName, "%s.conf", output_file);
  conf.open(confName, ios_base::out | ios_base::binary);
  if (! conf.is_open()) {
    cerr << "Can't write configuration file " << confName << endl;
    exit(3);
  }
  delete [] confName;
  
  // Writing configuration file
  list<pair<uint, String<char> > >::iterator it = sequences.begin();
  conf << sequences.size();
  while (it != sequences.end()) {
    conf << endl << (*it).second 
         << endl << (*it).first;
    it++;
  }
  conf.close();
  
  if (n > ((unsigned long long) 1 << 32) - 1) {
    cerr << "Your genome sequence is too large. I cannot index it. Sorry" << endl;
    exit(4);
  }

  if (verbose) 
    cerr << "Initialising KBSA" <<  endl;

  BlockwiseSA<TStr> *kSA = new KarkkainenBlockwiseSA<TStr>(*text, bucket_size, cover);
  
  if (verbose)
    cerr << "Building BWT. Do you know what the BWT is?" << endl
         << "You can have a look at http://enwp.org/Burrows-Wheeler%20transform" << endl;

  // Build the BWT of chars.
  bwt = new unsigned char[n+1];
  for (uint i = 0; i <= n; i++) {
    uint current = kSA->nextSuffix();
    if (current == 0)
      bwt[i] = 255;
    else {
      bwt[i] = convert<char>((*text)[current - 1]);
    }
  }
  delete kSA;
  delete text;

  if (verbose)
    cerr << "Building the BWT-based index." << endl
         << "Do you know what an FM-index is?" << endl
         << "You can have a look at http://enwp.org/FM-index" << endl;

  // build_index has been modified so that it takes a BWT as parameter
  char *params = new char[13];
  sprintf(params, "terminator %c",255);
  build_index((unsigned char *)bwt, n+1, params, &index);
  if (verbose)
    cerr << "Saving the index in the output file" << endl;
  save_index(index, output_file);
  delete [] params;
  free_index(index);
}


void get(char *input, char *output) {
  void *index;
  unsigned int n;
  uchar *chars;
  uint total_pos;
  uint current_pos;
  ofstream fasta;
  GenomeInfo *gInfo;

  if (verbose)  
    cerr << "Reading informations on sequences" << endl;
  gInfo = new GenomeInfo(input);

  fasta.open(output, ios_base::out | ios_base::binary);
  if (! fasta.is_open()) {
    cerr << "Can't write FASTA file " << output << endl;
    exit(3);
  }

  if (verbose)
    cerr << "Loading index" << endl;
  load_index(input,&index);
  get_length(index, &n);
  if (verbose)
    cerr << "Extracting the sequences from the index" << endl;
  extract(index, 0, n-1, &chars, &n);

  uint i = 0;
  uchar substr[NB_CHAR_FASTA_LINE+1];
  substr[NB_CHAR_FASTA_LINE]=0;
  total_pos = 0;
  while (i < gInfo->getNbChr()) {
    current_pos = 0;
    // Output FASTA header for the sequence
    if (verbose)
      cerr << "Retrieving " << gInfo->getChrName(i) << endl;
    fasta << ">" << gInfo->getChrName(i) << endl;
    // Reading each sequence
    while (current_pos + NB_CHAR_FASTA_LINE < gInfo->getChrLength(i)) {
      // Writing each line separately
      strncpy((char *)substr, (char *)&chars[total_pos], NB_CHAR_FASTA_LINE);
      fasta << substr << endl;
      total_pos += NB_CHAR_FASTA_LINE;
      current_pos += NB_CHAR_FASTA_LINE;
    }
    // Read the remaining (if any)
    if (current_pos < gInfo->getChrLength(i)) {
      strncpy((char *) substr, (char *) & chars[total_pos], gInfo->getChrLength(i) - current_pos);
      substr[gInfo->getChrLength(i) - current_pos] = 0;
      fasta << substr << endl;
      total_pos += gInfo->getChrLength(i) - current_pos;
    }
    i++;
  }
  fasta.close();
  free(chars);

  free_index(index);

  delete gInfo;
}


int main(int argc, char **argv) {
  uint bucket_size = THIS_BUCKET_SIZE;
  uint cover = COVER;
  char *output_file = NULL;
  int command=-1;                  // O: index, 1: get
  int option_index;
  int opt;

  static const char *short_options = "d:b:hv";
  static struct option long_options[] = 
    {{"help", no_argument, NULL, 'h'},
     {"diff-cover", required_argument, NULL, 'd'},
     {"bucket-size", required_argument, NULL, 'b'},
     { 0, 0, 0, 0 } //terminator
    };

  while ((opt = getopt_long(argc, argv, short_options, long_options, &option_index)) != -1){
    switch(opt) {
    case 'd' : 
      cover = atoi(optarg);
      break;
    case 'b' : 
      bucket_size = atoi(optarg);
      break;
    case 'v':
      verbose = 1;
      break;
    case 'h' : 
    default: 
      usage(argv[0]);
    }
  }

  if (argc - optind < 3) {
    usage(argv[0]);
  }

  // Get the command argument
  if (! strcmp(argv[optind], "index")) {
    command = 0;
  } else if (! strcmp(argv[optind], "get")) {
    command = 1;
  } else {
    cerr << "Unknown command " << argv[optind] << endl;
    usage(argv[0]);
  }
  optind++;

  // Get the output file
  output_file = argv[optind++];

  // Input files are what is remaining...

  if (command == 0) {
    // index
    index(&argv[optind], argc - optind, output_file, cover, bucket_size);
  } else {
    // get
    get(argv[optind], output_file);
  }

    exit(0);
}
