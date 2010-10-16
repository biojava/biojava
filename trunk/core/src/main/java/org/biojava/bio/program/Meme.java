/*
 *                    BioJava development code
 *
 * This code may be freely distributed and modified under the
 * terms of the GNU Lesser General Public Licence.  This should
 * be distributed with the code.  If you do not have a copy,
 * see:
 *
 *      http://www.gnu.org/copyleft/lesser.html
 *
 * Copyright for this code is held jointly by the individual
 * authors.  These should be listed in @author doc comments.
 *
 * For more information on the BioJava project and its aims,
 * or to join the biojava-l mailing list, visit the home page
 * at:
 *
 *      http://www.biojava.org/
 *
 */


package org.biojava.bio.program;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StreamTokenizer;
import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dp.SimpleWeightMatrix;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.ChangeVetoException;

/**
 * The results of a meme run.
 *
 * @author Matthew Pocock
 */
public class Meme {
  private List motifs;
  private List seqIDs;

  {
    motifs = new ArrayList();
    seqIDs = new ArrayList();
  }

  public List getMotifs() {
    return motifs;
  }

  public List getSeqIDs() {
    return seqIDs;
  }

  public Meme(InputStream is, SymbolTokenization symParser)
         throws IOException, IllegalSymbolException, IllegalAlphabetException {
    StreamTokenizer st = new StreamTokenizer(
      new BufferedReader(new InputStreamReader(is)));
    st.eolIsSignificant(true);
    st.wordChars('*', '*');
    st.parseNumbers();

    SymbolList sym = null;

   ALPHABET:
    while( true ) {
      int nt = st.nextToken();
      if (nt == StreamTokenizer.TT_EOF) {
          return;
      } else if (nt == StreamTokenizer.TT_WORD) {
          if(st.sval.startsWith("ALPHABET")) {
            while(st.nextToken() != StreamTokenizer.TT_WORD) {}
            sym = new SimpleSymbolList(symParser, st.sval);
            break ALPHABET;
          }
      }
    }

    while(st.nextToken() != StreamTokenizer.TT_EOL) {}
    while(st.nextToken() != StreamTokenizer.TT_EOL) {}

   SEQLIST:
    while( true ) {
      if(st.nextToken() == StreamTokenizer.TT_WORD) {
          if(st.sval != null && st.sval.startsWith("*"))
            break SEQLIST;

          //need this cause lines sometimes wrap!?
          if(! st.sval.startsWith("Length"))
           seqIDs.add(st.sval.intern());
      }
    }

   OUTER:
    while( true ) {
      int width = 0;

     FINDMOTIF:
      while( true ) {
        int nt = st.nextToken();
        if (nt == StreamTokenizer.TT_EOF) {
            break OUTER;
        } else if (nt == StreamTokenizer.TT_WORD) {
            if(st.sval.startsWith("MOTIF")) {
              st.nextToken();			// MOTIF x
              while(st.nextToken() != StreamTokenizer.TT_NUMBER) {} // width = w
              width = (int) st.nval;		// w
              break FINDMOTIF;
            }
        }
      }

     FINDWEIGHTS:
      while( true ) {
        int nt = st.nextToken();
        if (nt == StreamTokenizer.TT_EOF) {
            break OUTER;
        } else if (nt == StreamTokenizer.TT_WORD) {
            if(st.sval.startsWith("letter")) {
              while(st.nextToken() != StreamTokenizer.TT_EOL) {}
              break FINDWEIGHTS;
            }
        }
      }

      SimpleWeightMatrix matrix = new SimpleWeightMatrix(
        (FiniteAlphabet) symParser.getAlphabet(),
        width,
        DistributionFactory.DEFAULT
      );

      int r = 0;
      int c = 0;
     READMOTIF:
      while( true ) {
        int nt = st.nextToken();
        if (nt == StreamTokenizer.TT_EOF) {
            break OUTER;
        } else if (nt == StreamTokenizer.TT_EOL) {
            r = 0;
            c++;
            if(c == width)
              break READMOTIF;
        } else if (nt == StreamTokenizer.TT_NUMBER) {
          try {
            matrix.getColumn(c).setWeight(sym.symbolAt(r+1), st.nval);
            r++;
          } catch (ChangeVetoException cve) {
            throw new BioError("Couldn't set up the distribution ",cve);
          }
        }
      }

      motifs.add(matrix);
    }
  }
}
