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

package org.biojava.bio.program.phred;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.OutputStream;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import org.biojava.bio.BioError;
import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.dist.Distribution;
import org.biojava.bio.dist.DistributionFactory;
import org.biojava.bio.dist.DistributionTools;
import org.biojava.bio.dist.DistributionTrainerContext;
import org.biojava.bio.dist.SimpleDistributionTrainerContext;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.Sequence;
import org.biojava.bio.seq.SequenceIterator;
import org.biojava.bio.seq.db.HashSequenceDB;
import org.biojava.bio.seq.db.IDMaker;
import org.biojava.bio.seq.db.SequenceDB;
import org.biojava.bio.seq.impl.SimpleSequence;
import org.biojava.bio.seq.io.FastaDescriptionLineParser;
import org.biojava.bio.seq.io.FastaFormat;
import org.biojava.bio.seq.io.SeqIOTools;
import org.biojava.bio.seq.io.SequenceBuilderFactory;
import org.biojava.bio.seq.io.SimpleSequenceBuilder;
import org.biojava.bio.seq.io.StreamReader;
import org.biojava.bio.seq.io.StreamWriter;
import org.biojava.bio.seq.io.SymbolTokenization;
import org.biojava.bio.symbol.Alphabet;
import org.biojava.bio.symbol.AlphabetManager;
import org.biojava.bio.symbol.BasisSymbol;
import org.biojava.bio.symbol.FiniteAlphabet;
import org.biojava.bio.symbol.IllegalAlphabetException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.IntegerAlphabet;
import org.biojava.bio.symbol.SimpleSymbolList;
import org.biojava.bio.symbol.Symbol;
import org.biojava.bio.symbol.SymbolList;
import org.biojava.utils.AssertionFailure;
import org.biojava.utils.ChangeVetoException;
import org.biojava.utils.ListTools;

/**
 * <p>PhredTools contains static methods for working with phred
 * quality data.</p>
 *
 * <p>Copyright (c) 2001</p>
 * <p>Company:      AgResearch</p>
 *
 * @author Mark Schreiber
 * @author Matthew Pocock
 * @since 1.1
 *
 * Note that Phred is a copyright of CodonCode Corporation.
 */

public final class PhredTools {

   static{
     try {
       List l = new ArrayList(2);
       l.add(DNATools.getDNA());
       l.add(IntegerAlphabet.getSubAlphabet(0,99));
       AlphabetManager.getCrossProductAlphabet(l, "PHRED");
     } catch (IllegalAlphabetException iae) {
       throw new BioError( "Could not create the Phred alphabet", iae);
     }
   }

  /**
   * Retrieves the PHRED alphabet from the AlphabetManager. The Phred alphabet
   * is a cross product of a subset of the IntegerAlphabet from 0...99 and the
   * DNA alphabet. The Phred alphabet is finite.
   *
   * The Phred Alphabet contains 400 BasisSymbols named, for example, (guanine 47).
   * The BasisSymbols can be fragmented into their component AtomicSymbols using
   * the <code>getSymbols()</code> method of BasisSymbol.
   */
  public static final FiniteAlphabet getPhredAlphabet(){
    return (FiniteAlphabet)AlphabetManager.alphabetForName("PHRED");
  }

  /**
   * Retrives the DNA symbol component of the Phred BasisSymbol from the
   * PHRED alphabet.
   * @throws IllegalSymbolException if the provided symbol is not from the
   * PHRED alphabet.
   */
  public static final Symbol dnaSymbolFromPhred(Symbol phredSym)
    throws IllegalSymbolException{
    //validate the symbol
    getPhredAlphabet().validate(phredSym);
    //get the DNA component of the Phred Symbol
    List l = ((BasisSymbol)phredSym).getSymbols();
    //the first symbol should be DNA
    return (Symbol) l.get(0);
  }

  /**
   * Retrives the IntegerSymbol component of the Phred BasisSymbol from the
   * PHRED alphabet.
   * @throws IllegalSymbolException if the provided symbol is not from the
   * PHRED alphabet.
   */
  public static final IntegerAlphabet.IntegerSymbol integerSymbolFromPhred(Symbol phredSym)
    throws IllegalSymbolException{
    //validate the symbol
    getPhredAlphabet().validate(phredSym);
    //get the IntegerSymbol component of the Phred Symbol
    List l = ((BasisSymbol)phredSym).getSymbols();
    //the second symbol should be the IntegerSymbol
    return (IntegerAlphabet.IntegerSymbol)(l.get(1));
  }

  /**
   * Merges a Symbol List from the DNA alphabet with a SymbolList from the
   * [0..99] subset of the IntegerAlphabet into a SymbolList from
   * the PHRED alphabet.
   * @throws IllegalAlphabetException if the alphabets are not of the required alphabets
   * @throws IllegalArgumentException if the two SymbolLists are not of equal length.
   * @throws IllegalSymbolException if a combination of Symbols cannot be represented by
   * the PHRED alphabet.
   */
  public static SymbolList createPhred(SymbolList dna, SymbolList quality)
    throws IllegalArgumentException, IllegalAlphabetException, IllegalSymbolException{
    //perform initial checks
    if(dna.length() != quality.length()){
      throw new IllegalArgumentException("SymbolLists must be of equal length "+
        dna.length()+" : "+quality.length());
    }
    if(dna.getAlphabet() != DNATools.getDNA()){
      throw new IllegalAlphabetException(
        "Expecting SymbolList 'dna' to use the DNA alphabet, uses "
        +dna.getAlphabet().getName());
    }
    Alphabet subint = IntegerAlphabet.getSubAlphabet(0,99);
    if(quality.getAlphabet() != subint && quality.getAlphabet() != IntegerAlphabet.getInstance()){
      throw new IllegalAlphabetException(
        "Expecting SymbolList quality to use the "+subint.getName()+" alphabet"+
        "or IntegerAlphabet instead uses "+
        quality.getAlphabet().getName());
    }

    //build the symbollist
    SimpleSymbolList sl = new SimpleSymbolList(getPhredAlphabet());

    for(int i = 1; i <= dna.length(); i++){
      Symbol d = dna.symbolAt(i);
      Symbol q = quality.symbolAt(i);
      try{
        sl.addSymbol(getPhredSymbol(d,q));
      }catch(ChangeVetoException e){
        throw new AssertionFailure(e);
      }
    }

    return sl;
  }

  /**
   * Creates a symbol from the PHRED alphabet by combining a Symbol from the
   * DNA alphabet and a Symbol from the IntegerAlphabet (or one of its subsets).
   * @throws IllegalSymbolException if there is no Symbol in the PHRED alphabet
   * that represents the two arguments.
   */
  public static final Symbol getPhredSymbol(Symbol dna, Symbol integer)
    throws IllegalSymbolException{
    return getPhredAlphabet().getSymbol(new ListTools.Doublet(dna, integer));
  }

  /**
   * Writes Phred quality data in a Fasta type format.
   * @param db a bunch of PhredSequence objects
   * @param qual the OutputStream to write the quality data to.
   * @param seq the OutputStream to write the sequence data to.
   * @since 1.2
   */
   public static void writePhredQuality(OutputStream qual, OutputStream seq, SequenceDB db)
    throws IOException, BioException{
      StreamWriter qualw = new StreamWriter(qual,new PhredFormat());
      StreamWriter seqw = new StreamWriter(seq, new FastaFormat());
      SequenceDB qualDB = new HashSequenceDB(IDMaker.byName);
      //Get the quality SymbolLists and add them to a SeqDB
      for(SequenceIterator i = db.sequenceIterator(); i.hasNext();){
        Sequence p = i.nextSequence();
        if(p instanceof PhredSequence){
          PhredSequence ps = (PhredSequence)p;
          SymbolList ql = ps.getQuality();
          try{
            qualDB.addSequence( new SimpleSequence(ql,p.getURN(),p.getName(),p.getAnnotation()));
          }catch(ChangeVetoException cve){
            throw new AssertionFailure("Cannot Add Quality Sequences to Database", cve);
          }
        }
        else{
          throw new BioException("Expecting PhredSequence, got " + p.getClass().getName());
        }
      }
      qualw.writeStream(qualDB.sequenceIterator());
      seqw.writeStream(db.sequenceIterator());//this works as sequence methods act on the underlying SimpleSequence
   }

  /**
   * Constructs a StreamReader to read in Phred quality data in FASTA format.
   * The data is converted into sequences consisting of Symbols from the IntegerAlphabet.
   */
  public static StreamReader readPhredQuality(BufferedReader br){
    return new StreamReader(br,
      new PhredFormat(),
      getQualityParser(),
      getFastaBuilderFactory());
  }



  /**
   * Calls SeqIOTools.readFastaDNA(br), added here for convinience.
   */
  public static StreamReader readPhredSequence(BufferedReader br){
    return (StreamReader)SeqIOTools.readFastaDNA(br);
  }


  private static SequenceBuilderFactory _fastaBuilderFactory;

    /**
     * Get a default SequenceBuilderFactory for handling FASTA
     * files.
     */
  private static SequenceBuilderFactory getFastaBuilderFactory() {
      if (_fastaBuilderFactory == null) {
          _fastaBuilderFactory = new FastaDescriptionLineParser.Factory(SimpleSequenceBuilder.FACTORY);
      }
      return _fastaBuilderFactory;
  }

  /**
   * returns the IntegerAlphabet parser
   */
  private static SymbolTokenization getQualityParser() {
    return IntegerAlphabet.getInstance().getTokenization("token");
  }

  /**
   * The quality value is related to the base call error probability
   * by the formula  QV = - 10 * log_10( P_e )
   * where P_e is the probability that the base call is an error.
   * @return a <code>double</code> value, note that for most Phred scores this will be rounded
   * to the nearest <code>int</code>
   */
   public static double qualityFromP(double probOfError){
     return (-10 * (Math.log(probOfError)/Math.log(10.0)));
   }

   /**
    * Calculates the probability of an error from the quality score via the formula
    *  P_e = 10**(QV/-10)
    */
    public static double pFromQuality(double quality){
      return Math.pow(10.0,(quality/-10.0));
    }

    /**
     * Calculates the probability of an error from the quality score via the formula
     *  P_e = 10**(QV/-10)
     */
    public static double pFromQuality(int quality){
      return pFromQuality((double)quality);
    }

    /**
     * Calculates the probability of an error from the quality score via the formula
     *  P_e = 10**(QV/-10)
     */
    public static double pFromQuality(IntegerAlphabet.IntegerSymbol quality){
      return pFromQuality(quality.intValue());
    }

    /**
     * Converts a Phred sequence to an array of distributions. Essentially a fuzzy sequence
     * Assumes that all of the non called bases are equiprobable
     */
    public static Distribution[] phredToDistArray(PhredSequence s){
      Distribution[] pos = new Distribution[s.length()];
      DistributionTrainerContext dtc = new SimpleDistributionTrainerContext();

      for (int i = 0; i < s.length(); i++) {// for each symbol in the phred sequence
        Symbol qual = s.getQualityAt(i);
        Symbol base = s.getDNAAt(i);
        double pBase = pFromQuality((IntegerAlphabet.IntegerSymbol)qual);
        double pOthers = (1.0 - pBase)/3;

        try{
          pos[i] = DistributionFactory.DEFAULT.createDistribution(DNATools.getDNA());
          dtc.registerDistribution(pos[i]);

          for(Iterator iter = (DNATools.getDNA().iterator()); iter.hasNext();){
            Symbol sym = (Symbol)iter.next();
            if(sym.equals(base)) pos[i].setWeight(sym,pBase);
            else pos[i].setWeight(sym,pOthers);
          }

          dtc.train();
        }catch(IllegalAlphabetException iae){
          throw new AssertionFailure("Sequence "+s.getName()+" contains an illegal alphabet", iae);
        }catch(ChangeVetoException cve){
          throw new AssertionFailure("The Distribution has become locked", cve);
        }catch(IllegalSymbolException ise){
          throw new AssertionFailure("Sequence "+s.getName()+" contains an illegal symbol", ise);
        }
      }
      return pos;
    }

    /**
     * converts an Alignment of PhredSequences to a Distribution[] where each position is the average
     * distribution of the underlying column of the alignment.
     * @throws ClassCastException if the sequences in the alignment are not instances of PhredSequence
     */
    public static Distribution[] phredAlignmentToDistArray(Alignment a){
      List<String> labels = a.getLabels();
      Distribution [] average = new Distribution[a.length()];

      Distribution[][] matrix = new Distribution[labels.size()][];
      for(int y = 0; y < a.length(); y++){// for eaxh position
        for(Iterator<String> i = labels.iterator(); i.hasNext();){
          SymbolList sl = a.symbolListForLabel(i.next());
          matrix[y] = phredToDistArray((PhredSequence)sl);
        }
        average[y] = DistributionTools.average(matrix[y]);
      }

      return average;
    }
}
