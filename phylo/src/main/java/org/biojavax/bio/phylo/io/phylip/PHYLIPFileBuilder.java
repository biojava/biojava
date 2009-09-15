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
package org.biojavax.bio.phylo.io.phylip;

import java.util.ArrayList;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map.Entry;

import org.biojava.bio.BioException;
import org.biojava.bio.alignment.Alignment;
import org.biojava.bio.alignment.FlexibleAlignment;
import org.biojava.bio.alignment.SimpleAlignmentElement;
import org.biojava.bio.seq.DNATools;
import org.biojava.bio.seq.ProteinTools;
import org.biojava.bio.seq.io.ParseException;
import org.biojava.bio.symbol.IllegalSymbolException;
import org.biojava.bio.symbol.Location;
import org.biojava.bio.symbol.LocationTools;
import org.biojava.bio.symbol.SymbolList;

/**
 * Builds a PHYLIP file by listening to events.
 * 
 * @author Richard Holland
 * @author Tobias Thierer
 * @author Jim Balhoff
 * @since 1.6
 */
public class PHYLIPFileBuilder implements PHYLIPFileListener {
  
  private LinkedHashMap sequences;
  private int sequenceCount;
  private int sitesCount;
  private String currentSequenceName;
  private Alignment alignment;
  
  public void startFile() {
    this.sequences = new LinkedHashMap();
  }

  public void endFile() throws ParseException {
    this.verifySequenceAndSitesCount();
    this.buildAlignment();
  }

  public void setSequenceCount(int count) {
    this.sequenceCount = count;
  }
  
  public void setSitesCount(int count) {
    this.sitesCount = count;
  }
  
  public void setCurrentSequenceName(String name) {
    if (!(this.sequences.containsKey(name))) {
      sequences.put(name, new StringBuffer());
    }
    this.currentSequenceName = name;
  }
  
  public void receiveSequence(String sequence) {
    StringBuffer buffer = (StringBuffer)(this.sequences.get(this.currentSequenceName));
    //System.out.println(sequence);
    buffer.append(sequence);
  }
  
  public Alignment getAlignment() {
    return this.alignment;
  }
  
  private void buildAlignment() throws ParseException {
    List importedSequences = null;
    try {
      importedSequences = this.createSequences();
    } catch (IllegalSymbolException e) {
      throw new ParseException("Illegal symbol in sequence: " + e);
    } catch (BioException e) {
      throw new ParseException("Could not create sequences: " + e);
    }
    FlexibleAlignment newAlignment;
    try {
      newAlignment = new FlexibleAlignment(importedSequences);
    } catch (BioException e) {
      throw new ParseException("Could not construct alignment object: " + e);
    }
    this.alignment = newAlignment;
  }
  
  private List createSequences() throws IllegalSymbolException, BioException {
    List importedSequences = new ArrayList();
    boolean checkedType = false;
    boolean isDNA = true;
    Location loc = LocationTools.makeLocation(0, this.sitesCount - 1);
    for (Iterator i = this.sequences.entrySet().iterator(); i.hasNext(); ) {
      Entry sequenceEntry = (Entry)i.next();
      String name = (String)sequenceEntry.getKey();
      String sequence = ((StringBuffer)sequenceEntry.getValue()).toString();
      SymbolList symbolList = null;
      if (!checkedType) {
        try {
          DNATools.createGappedDNASequence(sequence, name);
        } catch (IllegalSymbolException e) {
          isDNA = false;
        }
        checkedType = true;
      }
      if (isDNA) {
        // make DNA sequences
        symbolList = DNATools.createGappedDNASequence(sequence, name);
        
      } else {
        // make protein sequences
        symbolList = ProteinTools.createGappedProteinSequence(sequence, name);
      }
      importedSequences.add(new SimpleAlignmentElement(name, symbolList, loc));
    }
    return importedSequences;
  } 
  
  private void verifySequenceAndSitesCount() throws ParseException {
    if (this.sequences.size() != this.sequenceCount) {
      throw new ParseException("Number of sequences does not match header.");
    } else {
      for (Iterator i = this.sequences.values().iterator(); i.hasNext();) {
        String currentSequence = ((StringBuffer)i.next()).toString();
        if (currentSequence.length() != this.sitesCount) {
          throw new ParseException("Number of sites does not match header.");
        }
      }
    }
  }
 }
