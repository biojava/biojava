/* -*- c-basic-offset: 2; indent-tabs-mode: nil -*- */
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

package org.biojava.bio.program.gff;

import java.util.Iterator;
import java.util.List;
import java.util.Map;

import org.biojava.bio.seq.StrandedFeature;
import org.biojava.utils.SmallMap;

/**
 * A no-frills implementation of a <span class="type">GFFRecord</span>.
 *
 * @author Matthew Pocock
 * @author Greg Cox
 * @author Aroul Ramadass
 * @author Len Trigg
 * @author Richard Holland
 */
public class SimpleGFFRecord implements GFFRecord {
  /**
   * The sequence name.
   */
  private String seqName;
  /**
   * The source.
   */
  private String source;
  /**
   * The feature type.
   */
  private String feature;
  /**
   * The start coordinate.
   */
  private int start;
  /**
   * The end coordinate.
   */
  private int end;
  /**
   * The feature score.
   */
  private double score;
  /**
   * The feature strand.
   */
  private StrandedFeature.Strand strand;
  /**
   * The feature frame.
   */
  private int frame;
  /**
   * The group-name -> <span class="type">List</span> &lt;attribute&gt;
   * <span class="type">Map</span>
   */
  private Map groupAttributes;
  /**
   * The comment.
   */
  private String comment;

  /**
  * Create a new SimpleGFFRecord from GFFRecord object
  * 
  * @param rec - A GFFRecord object
  */

  public SimpleGFFRecord(GFFRecord rec) {
    this.seqName = rec.getSeqName();
    this.source = rec.getSource();
    this.feature = rec.getFeature();
    this.start = rec.getStart();
    this.end = rec.getEnd();
    this.score = rec.getScore();
    this.strand = rec.getStrand();
    this.frame = rec.getFrame();
    this.comment = rec.getComment();
    this.groupAttributes = new SmallMap(rec.getGroupAttributes());
  }

  public SimpleGFFRecord(
                         String seqName,
                         String source,
                         String feature,
                         int start,
                         int end,
                         double score,
                         StrandedFeature.Strand strand,
                         int frame,
                         String comment,
                         Map groupAttributes
                         ) {
    this.seqName = seqName;
    this.source = source;
    this.feature = feature;
    this.start = start;
    this.end = end;
    this.score = score;
    this.strand = strand;
    this.frame = frame;
    this.comment = comment;
    this.groupAttributes = new SmallMap(groupAttributes);
  }


   /**
   * Create a new SimpleGFFRecord with values set to null or zero
   */
  public SimpleGFFRecord() {
    this.seqName = null;
    this.source = null;
    this.feature = null;
    this.start = 0;
    this.end = 0;
    this.score = 0;
    this.strand = null;
    this.frame = 0;
    this.comment = null;
    this.groupAttributes = null;
  }

  /**
   * Set the sequence name to <span class="arg">seqName</span>.
   *
   * @param seqName  the new name
   */
  public void setSeqName(String seqName) {
    this.seqName = seqName;
  }

  public String getSeqName() {
    return seqName;
  }

  /**
   * Set the feature source to <span class="arg">source</source>.
   *
   * @param source  the new source
   */
  public void setSource(String source) {
    this.source = source;
  }

  public String getSource() {
    return source;
  }

  /**
   * Set the feature type to <span class="arg">type</source>.
   *
   * @param feature  the new feature type
   */
  public void setFeature(String feature) {
    this.feature = feature;
  }

  public String getFeature() {
    return feature;
  }

  /**
   * Set the start coordinate to <span class="arg">start</source>.
   *
   * @param start  the new start coordinate
   */
  public void setStart(int start) {
    this.start = start;
  }

  public int getStart() {
    return start;
  }

  /**
   * Set the end coordinate to <span class="arg">end</source>.
   *
   * @param end  the new end coordinate
   */
  public void setEnd(int end) {
    this.end = end;
  }

  public int getEnd() {
    return end;
  }

  /**
   * Set the score to <span class="arg">score</source>.
   * <p>
   * The score must be a double, inclusive of <code>0</code>.
   * If you wish to indicate that there is no score, then use
   * <span class="type">GFFRecord</span>.<span class="const">NO_SCORE</span>.
   *
   * @param score  the new score
   */
  public void setScore(double score) {
    this.score = score;
  }

  public double getScore() {
    return score;
  }

  /**
   * Set the strand to <span class="arg">strand</source>.
   *
   * @param strand the new Strand
   */
  public void setStrand(StrandedFeature.Strand strand) {
    this.strand = strand;
  }

  public StrandedFeature.Strand getStrand() {
    return strand;
  }

  /**
   * Set the frame to <span class="arg">frame</source>.
   * <p>
   * The score must be  one of <code>{0, 1, 2}</code> or
   * <span class="type">GFFRecord</span>.<span class="const">NO_FRAME</span>.
   *
   * @param frame the frame
   * @throws IllegalArgumentException if score is not valid.
   */
  public void setFrame(int frame) {
    if (frame != GFFTools.NO_FRAME &&
       (frame < 0 || frame > 2))
    {
      throw new IllegalArgumentException("Illegal frame: " + frame);
    }
    this.frame = frame;
  }

  public int getFrame() {
    return frame;
  }

  /**
   * Replace the group-attribute <span class="type">Map</span> with
   * <span class="arg">ga</span>.
   * <p>
   * To efficiently add a key, call <span class="method">getGroupAttributes()</span>
   * and modify the <span class="type">Map</span>.
   *
   * @param ga  the new group-attribute <span class="type">Map</span>
   */
  public void setGroupAttributes(Map ga) {
    this.groupAttributes = ga;
  }

  public Map getGroupAttributes() {
    if (groupAttributes == null) {
      groupAttributes = new SmallMap();
    }
    return groupAttributes;
  }

  /**
   * Set the comment to <span class="arg">comment</source>.
   * <p>
   * If you set it to null, then the comment for this line will be ignored.
   *
   * @param comment the new comment
   */
  public void setComment(String comment) {
    this.comment = comment;
  }

  public String getComment() {
    return comment;
  }

  /**
   * Create a <span class="type">String</span> representation of
   * <span class="arg">attMap</span>.
   *
   * <span class="arg">attMap</span> is assumed to contain
   * <span class="type">String</span> keys and
   * <span class="type">List</span> values.
   *
   * @param attMap  the <span class="type">Map</span> of attributes and value lists
   * @return  a GFF attribute/value <span class="type">String</span>
   */
  public static String stringifyAttributes(Map attMap) {
    StringBuffer sBuff = new StringBuffer();
    Iterator ki = attMap.keySet().iterator();
    while (ki.hasNext()) {
      String key = (String) ki.next();
      sBuff.append(key);
      List values = (List) attMap.get(key);
      for (Iterator vi = values.iterator(); vi.hasNext();) {
        String value = (String) vi.next();
        if (isText(value)) {
          sBuff.append(" \"" + value + "\"");
        } else {
          sBuff.append(" " + value);
        }
      }
      if (ki.hasNext()) {
        sBuff.append(" ;");
      }
    }
    return sBuff.substring(0);
  }

  /**
   * Returns true if a string is "textual". The GFF Spec says that
   * "textual" values must be quoted. This implementation just tests
   * if the string contains letters or whitespace.
   *
   * @param value a <code>String</code> value.
   * @return true if value is "textual".
   */
  private static boolean isText(String value) {
    for (int i = 0; i < value.length(); i++) {
      char c = value.charAt(i);
      if (Character.isLetter(c) || Character.isWhitespace(c)) {
        return true;
      }
    }
    return false;
  }
}

