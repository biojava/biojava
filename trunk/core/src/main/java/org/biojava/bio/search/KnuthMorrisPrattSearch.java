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


package org.biojava.bio.search;


import java.util.ArrayList;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.seq.CircularView;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;

/**
 * An object to find exact subsequences within a sequence.
 *
 * <p>
 * Reference: KNUTH D.E., MORRIS (Jr) J.H., PRATT V.R., 1977,
 * Fast pattern matching in strings, SIAM Journal on Computing 6(1):323-350.
 * </P.
 *
 * <p><b>USAGE:</b></p>
 * When the object is constructed the <code>findMatches()</code>
 * method would be called. This will return an int[] giving the offsets
 * (ie the location of the first symbol of each match in the text).
 * The <code>getKMPNextTable()</code> returns the table of border lengths
 * used in the algorithm. This method is protected as it is unlikely it
 * will be needed except for debugging.</p>
 *
 * <p>
 * The algorithm finds exact matches therefore ambiguity symbols will match
 * only themselves. The class cannot perform regular expressions. The class
 * operates on all alphabets thus if searching for a DNA pattern you should
 * compile both the pattern and its reverse complement.</p>
 *
 * <p><b>WARNING the behaviour of a pattern containing gaps is undefined.
 *  It's not recommended that you try it.</b></p>
 *
 * <p>Copyright:    Copyright (c) 2002</p>
 * <p>Company:      AgResearch</p>
 *
 * @author Mark Schreiber
 * @version 1.0
 */

public final class KnuthMorrisPrattSearch {
  private int[] kmpNext;// the table that defines the border lengths.
  private SymbolList pattern;

  /**
   * Constructs a KMP matcher to find exact occurances of
   * <code>pattern</code> in <code>text</code> using the
   * Knuth-Morris-Pratt algorithm.<p>
   *
   * A new class should be constructed for each occurance of
   * <code>pattern</code>.<p>
   * @param pattern the pattern to search for
   */
  public KnuthMorrisPrattSearch(SymbolList pattern) {
    Alphabet alpha = pattern.getAlphabet();
    ArrayList rList = new ArrayList(pattern.toList());

    /*
     *need to perform this hack to make the kmpNext capable of dealing with
     *overlapping patterns, unfortunately it means the behaviour of a pattern
     *containing a gap is unspecified.
     */
    rList.add(alpha.getGapSymbol());
    try{
      this.pattern = new SimpleSymbolList(alpha, rList);
    }catch(Exception e){
      //really shouldn't happen
      throw new BioError("Couldn't make KMP pattern", e);
    }


    kmpNext = new int[this.pattern.length()];
    compilePattern();
  }

  private void compilePattern(){
    int k = pattern.length()-1;
    //System.out.println("k = "+k);
    int i = 0;
    int j = kmpNext[0] = -1;

    while(i < k){
      while (j > -1 && pattern.symbolAt(i+1) != pattern.symbolAt(j+1)){
        j = kmpNext[j];
      }
      i++; j++;
      if(pattern.symbolAt(i+1) == pattern.symbolAt(j+1)){
        //System.out.println("i = "+i+" j = "+j);
        kmpNext[i] = kmpNext[j];
      }else{
        //System.out.println("i = "+i+" j = "+j);
        kmpNext[i] = j;
      }
    }
  }



  /**
   * This will return an int[] giving the offsets of the matches in <code>text</code>
   * (ie the location of the first symbol of each match in the <code>text</code>).
   *
   * @param text the text or sequence to search for the pattern
   *
   * @return an int[] giving the offsets to the matches or a zero length array if there are
   * no matches.
   */
  public int[] findMatches(SymbolList text){
    List matches = new ArrayList();
    int n; // the length of the text
    int m; //the length of the pattern
    int i = 0;
    int j = 0;

    m = this.pattern.length()-1; //-1 to remove the gap at the end hack
    if(text instanceof CircularView){
      n = text.length()+pattern.length() -1; //allow wrap around
    }else{
      n = text.length();
    }

    //find the matches
    while(j < n){
      Symbol sym = text.symbolAt(j+1);
      while( i > -1 && pattern.symbolAt(i+1) != sym)
        i = kmpNext[i];
      i++;
      j++;
      if(i >= m){
        //match found, add 1 for SymbolList coordinates.
        matches.add(new Integer(j - i +1));
        i = kmpNext[i];
      }
    }

    //turn matches into an int[]
    int[] mat = new int[matches.size()];
    for (int x = 0; x < mat.length; x++) {
      mat[x] = ((Integer)matches.get(x)).intValue();
    }
    return mat;
  }

  /**
   * Returns the table of border lengths
   * @return an int[] of border lenghts
   */
  protected int[] getKmpNextTable(){
    return kmpNext;
  }

  /**
   *
   * @return the pattern being searched for
   */
  public SymbolList getPattern() {
      return pattern;
  }

  /**
   * Demo and Test method
   * @param args no arguments required
   * @throws Exception if the test fails
   */
  public static void main(String[] args) throws Exception{
    KnuthMorrisPrattSearch kmp1;
    int[] table;
    int[] matches;


    SymbolList pattern = DNATools.createDNA("gcagagag");
    SymbolList pattern2 = DNATools.createDNA("agag");
    SymbolList text = DNATools.createDNA("gcatcgcagagagtatacagtacg");

    //check pattern
    kmp1 = new KnuthMorrisPrattSearch(pattern);

    table = kmp1.getKmpNextTable();
    System.out.println(pattern.seqString());
    for (int i = 0; i < table.length; i++) {
      System.out.print(table[i] +" ");
    }
    //table should be -1 0 0 -1 1 -1 1 -1
    System.out.println("");
    matches = kmp1.findMatches(text);
    System.out.print("Matches at: ");
    for (int i = 0; i < matches.length; i++) {
      System.out.print(matches[i]+" ");
    }
    //matches should be at 6.
    System.out.println("\n");

    //check pattern2
    kmp1 = new KnuthMorrisPrattSearch(pattern2);
    table = kmp1.getKmpNextTable();
    System.out.println(pattern2.seqString());
    for (int i = 0; i < table.length; i++) {
      System.out.print(table[i] +" ");
    }
    System.out.println("");
    //table should be    -1  0 -1  0  2

    matches = kmp1.findMatches(text);
    System.out.print("Matches at: ");
    for (int i = 0; i < matches.length; i++) {
      System.out.print(matches[i]+" ");
    }
    //matches should be at 8 and 10
    System.out.println("\n");
  }
}
