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

package org.biojava.bio.program.tagvalue;


/**
 * <p>
 * A parser that splits a line into tag/value at a given column number. The
 * GENBANK and EMBL constants are parsers pre-configured for genbank and embl
 * style files respectively.
 * </p>
 *
 * <p>
 * There are many properties of the parser that can be set to change how lines
 * are split, and how the tag and value is produced from that split.
 * <ul>
 * <li>endOfRecord - string starting lines that mark record boundaries
 * e.g. "//"</li>
 * <li>splitOffset - column index of the first character of the value, and the
 * length of the raw tag e.g. 5 for EMBL files</li>
 * <li>trimTag - trim white-space from tags</li>
 * <li>trimValue - trim white-space from values</li>
 * <li>continueOnEmptyTag - if the tag is empty, use the previous tag e.g. this
 * is true for GENBANK files and false for EMBL files</li>
 * <li>mergeSameTag - if two consecutive tags have the same value, consider
 * their values to be a continuation of a single value so don't fire start/end
 * tag events e.g. true for EMBL</li>
 * </ul>
 *
 * @author Matthew Pocock
 * @author Keith James (enabled empty line EOR)
 * @since 1.2
 */
public class LineSplitParser
  implements
    TagValueParser,
    Cloneable
{
  /**
   * A LineSplitParser pre-configured to process EMBL-style flat files.
   */
  public static final LineSplitParser EMBL;

  /**
   * A LineSplitParser pre-configured to process GENBANK-style flat files.
   */
  public static final LineSplitParser GENBANK;
  
  static {
    EMBL = new LineSplitParser();
    EMBL.setEndOfRecord("//");
    EMBL.setSplitOffset(5);
    EMBL.setTrimTag(true);
    EMBL.setTrimValue(false);
    EMBL.setContinueOnEmptyTag(false);
    EMBL.setMergeSameTag(true);
    
    GENBANK = new LineSplitParser();
    GENBANK.setEndOfRecord("//");
    GENBANK.setSplitOffset(12);
    GENBANK.setTrimTag(true);
    GENBANK.setTrimValue(false);
    GENBANK.setContinueOnEmptyTag(true);
    GENBANK.setMergeSameTag(false);
  }

  // properties
  //
  
  private String endOfRecord = null;
  
  private int splitOffset;
  
  private boolean trimTag;
  
  private boolean trimValue;
  
  private boolean continueOnEmptyTag;
  
  private boolean mergeSameTag;

  // state
  //
  
  private String tag;
  
  public LineSplitParser() {}

  public LineSplitParser(LineSplitParser parser) {
    this.endOfRecord = parser.endOfRecord;
    this.splitOffset = parser.splitOffset;
    this.trimTag = parser.trimTag;
    this.trimValue = parser.trimValue;
    this.continueOnEmptyTag = parser.continueOnEmptyTag;
    this.mergeSameTag = parser.mergeSameTag;
  }
  
  /**
   * Set the string indicating that a record has ended.
   *
   * @param endOfRecord the new String delimiting records
   */
  public void setEndOfRecord(String endOfRecord) {
    this.endOfRecord = endOfRecord;
  }
  
  /**
   * Get the current string indicating that a record has ended.
   *
   * @return the current string delimiting records.
   */
  public String getEndOfRecord() {
    return endOfRecord;
  }
  
  /**
   * Set the offset to split lines at.
   *
   * @param splitOffset the new offset to split at
   */
  public void setSplitOffset(int splitOffset) {
    this.splitOffset = splitOffset;
  }
  
  /**
   * Get the current offset at which lines are split.
   *
   * @return the offset to split at
   */
  public int getSplitOffset() {
    return splitOffset;
  }
  
  /**
   * Enable or disable trimming of tags.
   *
   * @param trimTag  true if tags should be trimmed, otherwise false
   */
  public void setTrimTag(boolean trimTag) {
    this.trimTag = trimTag;
  }
  
  /**
   * See if tag trimming is enabled.
   *
   * @return true if tags are trimmed, otherwise false
   */
  public boolean getTrimTag() {
    return trimTag;
  }
  
  /**
   * Enable or disable trimming of values.
   *
   * @param trimValue  true if values should be trimmed, otherwise false
   */
  public void setTrimValue(boolean trimValue) {
    this.trimValue = trimValue;
  }
  
  /**
   * See if value trimming is enabled.
   *
   * @return true if values are trimmed, otherwise false
   */
  public boolean getTrimValue() {
    return trimValue;
  }
  
  /**
   * Choose whether to treat empty tags as a continuation of previous tags or as a
   * new tag with the value of the empty string.
   *
   * @param continueOnEmptyTag true to enable empty tags to be treated as a
   *        continuation of the previous tag, false otherwise
   */
  public void setContinueOnEmptyTag(boolean continueOnEmptyTag) {
    this.continueOnEmptyTag = continueOnEmptyTag;
  }
  
  /**
   * See if empty tags are treated as a continuation of previous tags or as a
   * new tag with the value of the empty string.
   *
   * @return true if continuation is enabled, false otherwise
   */
  public boolean getContinueOnEmptyTag() {
    return continueOnEmptyTag;
  }
  
  /**
   * Enable or disable treating runs of identical tags as a single tag start
   * event with multiple values or each as a separate tag start, value, and tag
   * end.
   *
   * @param mergeSameTag true if tags should be merged, false otherwise
   */
  public void setMergeSameTag(boolean mergeSameTag) {
    this.mergeSameTag = mergeSameTag;
  }
  
  /**
   * See if tags are being merged.
   *
   * @return true if merging is enabled, false otherwise
   */
  public boolean getMergeSameTag() {
    return mergeSameTag;
  }
  
  public TagValue parse(Object o) {
    String line = o.toString();

    // Use of the special value for the EOR marker allows a blank line
    // to be used to delimit records. Many file formats are like this.
    if (endOfRecord != null) {
        if (endOfRecord == TagValueParser.EMPTY_LINE_EOR) {
            if (line.equals(TagValueParser.EMPTY_LINE_EOR)) {
                return null;
            }
        }
        else
        {
            if (line.startsWith(endOfRecord)) {
                return null;
            }
        }
    }
    
    int length = line.length();
    
    String tag;
    if(length > splitOffset) {
      tag = line.substring(0, splitOffset);
    } else {
      tag = line;
    }
    if(trimTag) {
      tag = tag.trim();
    }
    
    String value;
    if(length > splitOffset) {
      value = line.substring(splitOffset);
    } else {
      value = "";
    }
    if(trimValue) {
      value = value.trim();
    }
    
    if(continueOnEmptyTag && (tag.length() == 0)) {
      return new TagValue(this.tag, value, false);
    } else if(mergeSameTag && tag.equals(this.tag)) {
      return new TagValue(tag, value, false);
    } else {
      return new TagValue(this.tag = tag, value, true);
    }
  }
  
  public Object clone()
  throws CloneNotSupportedException {
    return super.clone();
  }
}
